

setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")



# Downloading model -------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)

source("model.R")

# LHS ---------------------------------------------------------------------

library(lhs)
library(ggplot2)
library(dplyr)
library(stringr)
library(rlang)

source("param_ranges.R")
source("scenarios.R")

# source("sampling.R")
# load("simulation_dt/simulation_dt_50_2.RData")
# simulation_dt1 = simulation_dt
# load("simulation_dt/simulation_dt_150_2.RData")
# simulation_dt2 = simulation_dt
# load("simulation_dt/simulation_dt_300_2.RData")
# simulation_dt3 = simulation_dt
# load("simulation_dt/simulation_dt_500_2.RData")
# simulation_dt4 = simulation_dt
# 
# simulation_dt = rbind(simulation_dt1,
#                       simulation_dt2,
#                       simulation_dt3,
#                       simulation_dt4)

load("simulation_dt/simulation_dt_50_3.RData")


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
N_BINS = 6
evaluated_parameter = c("beta_I_colony", "rho_to_colony")
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
p1 = plot_heatmap_binned_diff(all_diff_results$BO_vs_HS, evaluated_parameter, param_ranges)+
  ggtitle("BO - HS")
 #+labs(x = "Transmission rate", y = "P(Fledging->Breeder)")
p2 = plot_heatmap_binned_diff(all_diff_results$BO_vs_RS, evaluated_parameter, param_ranges)+
  ggtitle("BO - RS")
  #+labs(x = "Transmission rate", y = "P(Fledging->Breeder)")
p3 = plot_heatmap_binned_diff(all_diff_results$BO_vs_PS, evaluated_parameter, param_ranges)+
  ggtitle("BO - PS")
 #+labs(x = "Transmission rate", y = "P(Fledging->Breeder)")
p4 = plot_heatmap_binned_diff(all_diff_results$BO_vs_P2, evaluated_parameter, param_ranges)+
  ggtitle("BO - P2")
  #+labs(x = "Transmission rate", y = "P(Fledging->Breeder)")


# gridExtra::grid.arrange(p0, grid::rectGrob(gp = grid::gpar(col = NA)),
#                         p1, p2,
#                         p3, p4,
#                         ncol = 2, nrow = 3,
#                         top = " ")

gridExtra::grid.arrange(p0, p2,
                        p3, p4,
                        ncol = 2, nrow = 2,
                        top = " ")


best_strat_dt = 
  data.frame(
  x_mid = all_diff_results$BO_vs_HS$x_mid,
  y_mid = all_diff_results$BO_vs_HS$y_mid,
  BO_vs_RS = all_diff_results$BO_vs_RS$diff,
  #BO_vs_PS = all_diff_results$BO_vs_PS$diff,
  BO_vs_P2 = all_diff_results$BO_vs_P2$diff) %>% 
  mutate(
    best_strat = apply(.[, 3:ncol(.)], 1, function(row) {
      if (all(row < 0)) {
        return("BO")  # Retourne "BO" si toutes les valeurs sont négatives
      } else {
        return(colnames(.)[which.max(row) + 2])  # Retourne le nom de la colonne avec la valeur max
      }
    })
    )

best_strat_dt


plot_heatmap_best_strat = function(data, params, param_ranges) {
  
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = best_strat)
    ) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "lightgray", color = NA))+
    labs(x = params[1], y = params[2], fill = "Best \n Strategy")
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  
  return(p)
}

plot_heatmap_best_strat(best_strat_dt,evaluated_parameter,param_ranges)


# Plot boxplots -----------------------------------------------------------




interval_size = 0.5
evaluated_parameter = "initial_number_infected_breeders_A"


boxplot_dt = simulation_dt  %>%
  subset(., scenario == "P2")%>%
  mutate(parameter = cut(
    get(evaluated_parameter),
    breaks = seq(
      param_ranges[[evaluated_parameter]][[1]][[1]],
      param_ranges[[evaluated_parameter]][[1]][[2]],
      by = interval_size
    ),
    include.lowest = TRUE
  )) %>%
  na.omit()


ggplot(boxplot_dt, aes(x = parameter, y = nb_adults_equi)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = evaluated_parameter, y = "ENLA")





