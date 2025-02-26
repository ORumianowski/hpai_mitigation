

setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")



library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)
library(lhs)
library(dplyr)
library(stringr)
library(rlang)

source("param_ranges.R")
source("scenarios.R")
scenarios = scenarios[c(2,4,5),]



load("simulation_dt/simulation_dt_1000it_cluster_1.RData")
simulation_dt1 = simulation_dt
load("simulation_dt/simulation_dt_1000it_ordi_1.RData")
simulation_dt2 = simulation_dt
load("simulation_dt/simulation_dt_10000it_cluster_1.RData")
simulation_dt3 = simulation_dt
load("simulation_dt/simulation_dt_2000it_cluster_1.RData")
simulation_dt4 = simulation_dt
simulation_dt = rbind(simulation_dt1,
                      simulation_dt2,
                      simulation_dt3,
                      simulation_dt4)


# Plot binned heatmap -----------------------------------------------------
# Function to create binned data with dynamic parameters and variable block sizes
create_binned_data = function(data,
                              params,
                              param_ranges_,
                              selected_output,
                              n_bins) {
  
  # Apply log transformation if specified
  if (param_ranges_[[params[1]]][[2]] == "logarithmic") {
    data$var1 = log(data[[params[1]]])
  } else {
    data$var1 = data[[params[1]]]
  }
  
  if (param_ranges_[[params[2]]][[2]] == "logarithmic") {
    data$var2 = log(data[[params[2]]])
  } else {
    data$var2 = data[[params[2]]]
  }
  
  # Calculate min and max for both parameters
  min_x = min(data$var1, na.rm = TRUE)
  max_x = max(data$var1, na.rm = TRUE)
  min_y = min(data$var2, na.rm = TRUE)
  max_y = max(data$var2, na.rm = TRUE)
  
  # Define dynamic block sizes
  block_size_x = (max_x - min_x) / n_bins
  block_size_y = (max_y - min_y) / n_bins
  
  # Convert selected_output to a symbol
  selected_output_sym = sym(selected_output)
  
  res = data %>%
    mutate(
      V1 = cut(data$var1, breaks = seq(min_x, max_x, by = block_size_x), include.lowest = TRUE),
      V2 = cut(data$var2, breaks = seq(min_y, max_y, by = block_size_y), include.lowest = TRUE)
    ) %>%
    group_by(V1, V2) %>%
    summarise(output_avg = mean(!!selected_output_sym, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    mutate(
      # Use str_extract to extract the first numeric value from V1 and V2
      x_mid = as.numeric(str_extract(V1, "[-+]?[0-9]*\\.?[0-9]+")) + block_size_x / 2,
      y_mid = as.numeric(str_extract(V2, "[-+]?[0-9]*\\.?[0-9]+")) + block_size_y / 2
    )
  
  # Apply exp transformation if specified
  if (param_ranges_[[params[1]]][[2]] == "logarithmic") {
    res$x_mid = exp(res$x_mid)
  } else {
    res$x_mid = res$x_mid
  }
  
  if (param_ranges_[[params[2]]][[2]] == "logarithmic") {
    res$y_mid = exp(res$y_mid)
  } else {
    res$y_mid = res$y_mid
  }
  
  return(res)
}

# Define a function to create the binned data and calculate the difference
create_diff_between_scenarios <- function(scenario_1,
                                          scenario_2,
                                          simulation_dt,
                                          evaluated_parameter,
                                          param_ranges,
                                          selected_output,
                                          n_bins) {
  
  # Create binned data for scenario 1
  data_binned_1 <- create_binned_data(data = simulation_dt %>% subset(., scenario == scenario_1),
                                      params = evaluated_parameter,
                                      param_ranges,
                                      selected_output,
                                      n_bins) %>%
    dplyr::select(output_avg, x_mid, y_mid)
  
  # Create binned data for scenario 2
  data_binned_2 <- create_binned_data(simulation_dt %>% subset(., scenario == scenario_2),
                                      evaluated_parameter,
                                      param_ranges,
                                      selected_output,
                                      n_bins = n_bins) %>%
    dplyr::select(output_avg, x_mid, y_mid)
  
  # Merge the two datasets and calculate the difference
  diff_data <- merge(data_binned_1, data_binned_2, by = c("x_mid", "y_mid")) %>%
    mutate(diff = output_avg.x - output_avg.y) %>%
    dplyr::select(x_mid, y_mid, diff)
  
  return(diff_data)
}

# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned = function(data, params, param_ranges) {
  
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = output_avg)
    ) +
    scale_fill_gradient(low = "blue", high = "yellow") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "lightgray", color = NA))+
    labs(x = params[1], y = params[2], fill = "ENLA")
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  
  return(p)
}

# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned_diff = function(data, params, param_ranges) {
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = diff)
    ) +
    scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0)+
    theme_minimal() +
    theme(panel.background = element_rect(fill = "lightgray", color = NA)) +
    labs(x = params[1], y = params[2], fill = "ENLA")
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  
  return(p)
}


#   -----------------------------------------------------------------------

################################################################
N_BINS = 5
evaluated_parameter = c("theta", "beta_I_colony")
SELECTED_OUTPUT = "nb_adults_equi"

################################################################

# List of scenarios
scenarios_to_compare <- rownames(scenarios) 

# Empty list to store results
all_diff_results <- list()

# Loop through scenarios and calculate differences
for (scenario_name in scenarios_to_compare) {
  if (scenario_name != "BO") {  
    diff_result <- create_diff_between_scenarios("BO",
                                                 scenario_name, 
                                                 simulation_dt, 
                                                 evaluated_parameter,
                                                 param_ranges, 
                                                 selected_output = SELECTED_OUTPUT, 
                                                 n_bins = N_BINS)
    all_diff_results[[paste("BO_vs", scenario_name, sep = "_")]] <- diff_result
  }
}


# Create binned data with dynamic parameters and dynamic block sizes
plot_data_binned = create_binned_data(data = simulation_dt %>% subset(., scenario == "BO"),
                                      params = evaluated_parameter,
                                      param_ranges, 
                                      selected_output = SELECTED_OUTPUT,
                                      n_bins = N_BINS)



# Visualize the heatmap with dynamic block sizes
p0 = plot_heatmap_binned(plot_data_binned, evaluated_parameter, param_ranges)+
  ggtitle("BO")
  #+labs(x = "Transmission rate", y = "P(Fledging->Breeder)")
p2 = plot_heatmap_binned_diff(all_diff_results$BO_vs_RS, evaluated_parameter, param_ranges)+
  ggtitle("BO - RS")
  #+labs(x = "Transmission rate", y = "P(Fledging->Breeder)")
p4 = plot_heatmap_binned_diff(all_diff_results$BO_vs_P2, evaluated_parameter, param_ranges)+
  ggtitle("BO - P2")
  #+labs(x = "Transmission rate", y = "P(Fledging->Breeder)")



gridExtra::grid.arrange(p0, p2,
                        p4,
                        ncol = 2, nrow = 2,
                        top = " ")


# -------------------------------------------------------------------------

# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned_diff = function(data, params, param_ranges) {
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = diff)
    ) +
    scale_fill_gradient2(low = "red", mid = "white", high = "green",
                         midpoint = 0,
                         limits = c(-10, 35))+
    theme_minimal() +
    theme(panel.background = element_rect(fill = "lightgray", color = NA)) +
    labs(
      #x = params[1], y = params[2],
      x = NULL, y = NULL,
         fill = "ENLA")+
    guides(fill = "none") 
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  
  return(p)
}



N_BINS = 10
evaluated_parameter = c("beta_I_colony", 
                        "initial_number_infected_breeders_A", 
                        "theta",
                        "avrg_stay_NB_sea",
                        "infectious_period")

SELECTED_OUTPUT = "nb_adults_equi"


library(gridExtra)
library(grid)

# List of parameters

# Empty list to store results
all_diff_results <- list()


# Loop through scenarios and calculate differences
for (index_param_1 in 1:(length(evaluated_parameter))) {
  for (index_param_2 in 1:(length(evaluated_parameter))) {  
    
    
    diff_result <- create_diff_between_scenarios("BO",
                                                 "P2", 
                                                 simulation_dt, 
                                                 c(evaluated_parameter[index_param_1],
                                                   evaluated_parameter[index_param_2]),
                                                 param_ranges, 
                                                 selected_output = SELECTED_OUTPUT, 
                                                 n_bins = N_BINS)
    all_diff_results[[paste(evaluated_parameter[index_param_1], evaluated_parameter[index_param_2], sep = "__")]] <- list(value = diff_result,
                                                                                                                          param = c(evaluated_parameter[index_param_1],
                                                                                                                                    evaluated_parameter[index_param_2])                                                                                                                          )
  }
}


plot_heatmap_binned_diff(all_diff_results[[1]]$value, all_diff_results[[1]]$param, param_ranges)+
  ggtitle("BO - P2")


# Create a list of plots with lapply
plots <- lapply(1:length(all_diff_results), function(i) {
  plot_heatmap_binned_diff(all_diff_results[[i]]$value,
                           all_diff_results[[i]]$param,
                           param_ranges)
})

library(gridExtra)
library(grid)

# List of column labels
col_labels_ <- c("Transmission rate", 
                 "Inital infected", 
                 "Connectivity",
                 "Avrg stay NB sea",
                 "Infectious period")

# List of row labels
row_labels_ <- col_labels_

# Create text grobs for column labels (rotated 45 degrees)
col_labels <- lapply(1:(length(evaluated_parameter)), function(i) {
  textGrob(label = col_labels_[i], rot = 45, gp = gpar(fontsize = 10, fontface = "bold"))
})

# Create text grobs for row labels (rotated 45 degrees for the left side)
row_labels <- lapply(1:(length(evaluated_parameter)), function(i) {
  textGrob(label = row_labels_[i], rot = 45, gp = gpar(fontsize = 10, fontface = "bold"))
})

# Create layout matrix
n <- length(evaluated_parameter)

# Create an empty matrix for layout with n+1 rows and columns for the labels and plots
layout_matrix <- matrix(0, nrow = n+1, ncol = n+1)

# Fill the layout matrix
layout_matrix[1, 2:(n+1)] <- 1:n   # Top row for column labels
layout_matrix[2:(n+1), 1] <- (n+1):(2*n)  # First column for row labels
layout_matrix[2:(n+1), 2:(n+1)] <- (2*n+1):(n*n+2*n)  # The rest for the plots

# Combine nullGrob, row labels, col labels, and plots into a single list of grobs
all_grobs <- c(list(nullGrob()), col_labels, row_labels, plots)






# empty_plot <- ggplot() + theme_void()
# 
# gridExtra::grid.arrange(grobs = list(plots[[1]], plots[[2]],plots[[3]],plots[[4]],plots[[5]],
#                                      empty_plot,
#                                      plots[[6]], plots[[7]],plots[[8]],plots[[9]],
#                                      empty_plot, empty_plot,
#                                      plots[[10]], plots[[11]],plots[[12]],
#                                      empty_plot, empty_plot, empty_plot,
#                                      plots[[13]], plots[[14]],plots[[15]]
#                                      ),
#                         ncol = (length(evaluated_parameter)-1)
#                         , nrow = 5, top = " ")





# -------------------------------------------------------------------------





# 
# 
# best_strat_dt =
#   data.frame(
#   x_mid = all_diff_results$BO_vs_P2$x_mid,
#   y_mid = all_diff_results$BO_vs_P2$y_mid,
#   BO_vs_RS = all_diff_results$BO_vs_RS$diff,
#   BO_vs_P2 = all_diff_results$BO_vs_P2$diff) %>%
#   mutate(
#     best_strat = apply(.[, 3:ncol(.)], 1, function(row) {
#       if (all(row < 0)) {
#         return("BO")  # Retourne "BO" si toutes les valeurs sont négatives
#       } else {
#         return(colnames(.)[which.max(row) + 2])  # Retourne le nom de la colonne avec la valeur max
#       }
#     })
#     )
# 
# best_strat_dt
# 
# 
# plot_heatmap_best_strat = function(data, params, param_ranges) {
# 
#   p = ggplot(data) +
#     geom_tile(aes(x = x_mid, y = y_mid, fill = best_strat)
#     ) +
#     theme_minimal() +
#     theme(panel.background = element_rect(fill = "lightgray", color = NA))+
#     labs(x = params[1], y = params[2], fill = "Best \n Strategy")
# 
#   if (param_ranges[[params[1]]][[2]] == "logarithmic"){
#     p = p + scale_x_log10()
#   }
#   if (param_ranges[[params[2]]][[2]] == "logarithmic"){
#     p = p + scale_y_log10()
#   }
# 
#   return(p)
# }
# 
# plot_heatmap_best_strat(best_strat_dt,evaluated_parameter,param_ranges)
# 

# Plot boxplots -----------------------------------------------------------


plot_one_param = function(df, evaluated_scenario, evaluated_param){
  
  # Filter and process the data
  mono_dt = df  %>%
    subset(., scenario == evaluated_scenario) 
  
  # Dynamic plotting based on the evaluated parameter
  p = ggplot() +
    geom_point(data=mono_dt,
               aes_string(x = evaluated_param, y = "nb_adults_equi"),  # Couleur dans aes()
               size = 0.4, alpha = 0.14) +  # Transparence appliquée
    #geom_smooth(method = "loess") +
    theme_minimal() +
    labs(x = NULL, y = NULL)
  
  return(p)
}

plot_one_param(df = simulation_dt,
               evaluated_scenario = "P2",
               evaluated_param = evaluated_parameter[1])



plot_one_param_bank = list()

for (i in 1:length(evaluated_parameter)){
  
  plot_one_param_bank[[i]] = plot_one_param(df = simulation_dt,
                           evaluated_scenario = "P2",
                           evaluated_param = evaluated_parameter[i])
}



#   -----------------------------------------------------------------------

for (i in 1:length(evaluated_parameter)){
  all_grobs[[6+(6*i)]] = plot_one_param_bank[[i]]
}


grid.arrange(
  grobs = all_grobs,  
  layout_matrix = layout_matrix, 
  top = "BO - P2"
)



