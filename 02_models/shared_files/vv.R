

setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models/shared_files")



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




load("simulation_dt/dt_cluster_3_1E4_for_v2.RData")

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


dt %>% subset(., prop_dispersal>=0.99)

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

SCENARIO = "P2"

N_BINS = 10


evaluated_parameter = c(
  "beta_I_colony", 
  "initial_number_infected_breeders_A", 
  # "theta",
  # "prop_dispersal",
  # "avrg_stay_NB_sea",
  # "hatching_date",
  "reaching_repro_prob")



graph_param_name = c(
  "Transmission rate", 
  "Inital infected", 
  # "Connectivity",
  # "Prop. dispersed",
  # "Avrg. stay. NB sea",
  # "Hatching date",
  "Reaching Repro. Prob.")

# graph_param_name = evaluated_parameter

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
    summarise(output_q = quantile(!!selected_output_sym, na.rm = TRUE, probs = 0.5),
              output_avg = mean(!!selected_output_sym, na.rm = TRUE), .groups = 'drop') %>%
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
    dplyr::select(x_mid, y_mid, output_q, output_avg)
  
  return(res)
}




# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned_diff_mean = function(data_, params, param_ranges) {
  
  
  p = ggplot() +
    geom_tile(data = data_, aes(x = x_mid, y = y_mid, fill = output_avg)
    ) +
    scale_fill_gradient2(low = "brown3", mid = "white", high = "chartreuse4",
                         midpoint = 0,
                         limits = c(min(data_$output_avg), max(data_$output_avg)))+
    labs(
      x = graph_param_name[1], y = graph_param_name[2],
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

plot_heatmap_binned_diff_q = function(data_, params, param_ranges) {
  
  
  p = ggplot() +
    geom_tile(data = data_, aes(x = x_mid, y = y_mid, fill = output_q)
    ) +
    scale_fill_gradient2(low = "deepskyblue4", mid = "white", high = "chartreuse4",
                         midpoint = 0,
                         limits = c(min(data_$output_q), max(data_$output_q)))+
    labs(
      x = graph_param_name[1], y = graph_param_name[2],
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






# Empty list to store results
all_diff_results <- list()

for (index_param_1 in 1:(length(evaluated_parameter))) {
  for (index_param_2 in 1:(length(evaluated_parameter))) {  
    if (index_param_1 <= index_param_2) {
      diff_result = create_binned_data(data = simulation_dt,
                                         params = c(evaluated_parameter[index_param_1], evaluated_parameter[index_param_2]),
                                         param_ranges_ = param_ranges,
                                         selected_output = SELECTED_OUTPUT,
                                         n_bins = N_BINS)
    } else {
      diff_result = create_binned_data(data = simulation_dt,
                                            params = c(evaluated_parameter[index_param_2], evaluated_parameter[index_param_1]),
                                            param_ranges_ = param_ranges,
                                            selected_output = SELECTED_OUTPUT,
                                            n_bins = N_BINS)
    }
    
    all_diff_results[[paste(evaluated_parameter[index_param_1], evaluated_parameter[index_param_2], sep = "__")]] <- list(value = diff_result,
                                                                                                                          param = c(evaluated_parameter[index_param_1], evaluated_parameter[index_param_2]))
  }
}

# Create the heatmap matrix layout with correct functions
plots <- list()
nb_param <- length(evaluated_parameter)
for (i in 1:nb_param) {
  for (j in 1:nb_param) {
    if (i == j) {
      # Diagonal: Scatterplots
      plots[[length(plots) + 1]] <- plot_one_param_bank[[i]]
    } else if (i < j) {
      # Upper triangle: Mean heatmaps
      plots[[length(plots) + 1]] <- plot_heatmap_binned_diff_q(
        all_diff_results[[paste(evaluated_parameter[i], evaluated_parameter[j], sep = "__")]]$value,
        c(evaluated_parameter[i], evaluated_parameter[j]),
        param_ranges
      )
    } else {
      # Lower triangle: 5th percentile quantile heatmaps
      plots[[length(plots) + 1]] <- plot_heatmap_binned_diff_mean(
        all_diff_results[[paste(evaluated_parameter[i], evaluated_parameter[j], sep = "__")]]$value,
        c(evaluated_parameter[j], evaluated_parameter[i]),
        param_ranges
      )
    }
  }
}

# # Define the layout matrix
layout_matrix <- matrix(1:(nb_param*nb_param), nrow = nb_param, ncol = nb_param )
# layout_matrix[1, 2:(nb_param)] <- 1:nb_param   # Column labels
# layout_matrix[2:(nb_param), 1] <- (nb_param + 1):(2 * nb_param)  # Row labels
# layout_matrix[2:(nb_param), 2:(nb_param + 1)] <- (2 * nb_param + 1):(nb_param * nb_param + 2 * nb_param)
# 



# Generate the final plot
p = grid.arrange(
  grobs = plots, #all_grobs,  
  layout_matrix = layout_matrix, 
  top = paste0("", SCENARIO)
)

plot(p)
