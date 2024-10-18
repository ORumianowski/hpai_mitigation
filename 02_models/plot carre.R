

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
library(gridExtra)
library(grid)

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

param = c("initial_number_infected_breeders_A", "initial_number_infected_breeders_B", "initial_number_infected_breeders_C", "initial_number_breeders_A",         
          "initial_number_breeders_B",          "initial_number_breeders_C",          "dispersal_reaction_time",            "dispersal_date",                    
          "hatching_date",                      "tau",                                "total_time",                         "prop_dispersal",                    
          "beta_E_colony",                      "beta_I_colony",                      "incubation_period",                  "eta",                               
          "infectious_period",                  "adult_mortality",                    "nestling_mortality",                 "avrg_stay_B_colony",                
          "avrg_stay_B_sea",                    "avrg_stay_NB_colony",                "avrg_stay_NB_sea",                   "theta",                             
          "psi",                                "hatching_sd",                        "reaching_repro_prob",                "prob_detection"  )

param_evaluated = c("BO_P2_nb_adults_equi",
                    "initial_number_infected_breeders_A",        
                    "hatching_date",                      "prop_dispersal",                    
                    "beta_I_colony",                      "incubation_period",                                              
                    "infectious_period",                  "adult_mortality",   
                    "nestling_mortality",                         
                    "avrg_stay_NB_sea",                   "theta",                             
                    "reaching_repro_prob",                "prob_detection"  )

dt = simulation_dt[,1:32]


dt2 <- dt %>%
  pivot_wider(
    names_from = scenario,  
    values_from = c(nb_adults_equi, nb_infected_colonies, infected_X_time), 
    names_glue = "{scenario}_{.value}" 
  ) %>% 
  unnest(cols = everything()) %>% 
  mutate(BO_P2_nb_adults_equi = BO_nb_adults_equi - P2_nb_adults_equi) %>% 
  mutate(BO_RS_nb_adults_equi = BO_nb_adults_equi - RS_nb_adults_equi)

# -------------------------------------------------------------------------

SELECTED_OUTPUT = "BO_P2_nb_adults_equi"

SCENARIO = "RS"

N_BINS = 10

evaluated_parameter = c(
  "beta_I_colony", 
  "initial_number_infected_breeders_A", 
  "theta",
  "avrg_stay_NB_sea",
  "hatching_date",
  "reaching_repro_prob"
)



graph_param_name = c(
  "Transmission rate", 
  "Inital infected", 
  "Connectivity",
  "hatching_date",
  "Infectious period")

graph_param_name = evaluated_parameter

dt2 = dt2 %>% subset(., initial_number_infected_breeders_A!=0)



# -------------------------------------------------------------------------

simulation_dt = dt2


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
  
  res = res %>% 
    dplyr::select(x_mid, y_mid, output_avg)
  
  return(res)
}



# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned_diff = function(data, params, param_ranges) {
  
  
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = output_avg)
    ) +
    scale_fill_gradient2(low = "brown3", mid = "white", high = "chartreuse4",
                         midpoint = 0,
                         limits = c(min(data$output_avg), max(data$output_avg)))+
    labs(
      #x = params[1], y = params[2],
      x = NULL, y = NULL,
      fill = "ENLA")+
    guides(fill = "none") +
    theme_minimal()
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  
  return(p)
}


# Empty list to store results
all_diff_results <- list()

# Loop through scenarios and calculate differences
for (index_param_1 in 1:(length(evaluated_parameter))) {
  for (index_param_2 in 1:(length(evaluated_parameter))) {  
    

    
    diff_result = create_binned_data(data = simulation_dt,
                       params = c(evaluated_parameter[index_param_1],
                                  evaluated_parameter[index_param_2]),
                       param_ranges_ = param_ranges ,
                       selected_output = SELECTED_OUTPUT,
                       n_bins = N_BINS)
    
    all_diff_results[[paste(evaluated_parameter[index_param_1], evaluated_parameter[index_param_2], sep = "__")]] <- list(value = diff_result,
                                                                                                                          param = c(evaluated_parameter[index_param_1],
                                                                                                                                    evaluated_parameter[index_param_2])                                                                                                                          )
  }
}


# Create a list of plots with lapply
plots <- lapply(1:length(all_diff_results), function(i) {
  plot_heatmap_binned_diff(all_diff_results[[i]]$value,
                           all_diff_results[[i]]$param,
                           param_ranges)
})

# List of column labels
col_labels_ <- graph_param_name

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


plot_one_param = function(df, evaluated_param){

  
  # Dynamic plotting based on the evaluated parameter
  p = ggplot() +
    geom_point(data=df,
               aes_string(x = evaluated_param, y = SELECTED_OUTPUT),  # Couleur dans aes()
               size = 0.4, alpha = 0.14) +  # Transparence appliquée
    #geom_smooth(method = "loess") +
    theme_minimal() +
    labs(x = NULL, y = NULL)
  
  return(p)
}


plot_one_param_bank = list()

for (i in 1:length(evaluated_parameter)){
  
  plot_one_param_bank[[i]] = plot_one_param(df = simulation_dt,
                                            evaluated_param = evaluated_parameter[i])
}



#   -----------------------------------------------------------------------

nb_param = length(evaluated_parameter)

for (i in 1:length(evaluated_parameter)){
  all_grobs[[(nb_param+1)+((nb_param+1)*i)]] = plot_one_param_bank[[i]]
}


grid.arrange(
  grobs = all_grobs,  
  layout_matrix = layout_matrix, 
  top = paste0("BO - ", SCENARIO)
)



