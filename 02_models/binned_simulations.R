

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
load("simulation_dt_4.RData")


# Plot binned heatmap -----------------------------------------------------

evaluated_parameter = c("beta_I_colony", "reaching_repro_prob")
plotted_scenario = "BO"

# Number of bins (same for both dimensions)
n_bins = 10  


data = simulation_dt %>% subset(., scenario == "BO")

# Function to create binned data with dynamic parameters and variable block sizes
create_binned_data = function(data,
                              params, param_ranges,
                              selected_output,
                              n_bins) {
  
  # Apply log transformation if specified
  if (param_ranges[[params[1]]][[2]] == "logarithmic") {
    data$var1 = log(data[[params[1]]])
  } else {
    data$var1 = data[[params[1]]]
  }
  
  if (param_ranges[[params[2]]][[2]] == "logarithmic") {
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
  selected_output_sym = ensym(selected_output)
  
  
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
  if (param_ranges[[params[1]]][[2]] == "logarithmic") {
    res$x_mid = exp(res$x_mid)
  } else {
    res$x_mid = res$x_mid
  }
  
  if (param_ranges[[params[2]]][[2]] == "logarithmic") {
    res$y_mid = exp(res$y_mid)
  } else {
    res$y_mid = res$y_mid
  }
  
  return(res)
}

# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned = function(data, params, param_ranges) {
  
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = output_avg)
    ) +
    scale_fill_gradient(low = "blue", high = "yellow") +
    theme_minimal() +
    labs(x = params[1], y = params[2], fill = "ENLA")
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  
  return(p)
}

# Create binned data with dynamic parameters and dynamic block sizes
plot_data_binned = create_binned_data(simulation_dt %>% subset(., scenario == "BO"),
                                      evaluated_parameter, param_ranges,
                                      selected_output = "nb_infected_colonies",
                                      n_bins)



# Visualize the heatmap with dynamic block sizes
p0 = plot_heatmap_binned(plot_data_binned, 
                         evaluated_parameter, param_ranges)


# -------------------------------------------------------------------------

# Define a function to create the binned data and calculate the difference
create_diff_between_scenarios <- function(scenario_1,
                                          scenario_2,
                                          simulation_dt,
                                          evaluated_parameter,
                                          param_ranges,
                                          n_bins) {
  
  # Create binned data for scenario 1
  data_binned_1 <- create_binned_data(simulation_dt %>% subset(., scenario == scenario_1),
                                      evaluated_parameter,
                                      param_ranges,
                                      n_bins) %>%
    dplyr::select(nb_adults_equi_avg, x_mid, y_mid)
  
  # Create binned data for scenario 2
  data_binned_2 <- create_binned_data(simulation_dt %>% subset(., scenario == scenario_2),
                                      evaluated_parameter,
                                      param_ranges,
                                      n_bins) %>%
    dplyr::select(nb_adults_equi_avg, x_mid, y_mid)
  
  # Merge the two datasets and calculate the difference
  diff_data <- merge(data_binned_1, data_binned_2, by = c("x_mid", "y_mid")) %>%
    mutate(diff = nb_adults_equi_avg.x - nb_adults_equi_avg.y) %>%
    dplyr::select(x_mid, y_mid, diff)
  
  return(diff_data)
}

scenarios_to_compare <- rownames(scenarios) # List of scenarios

# Empty list to store results
all_diff_results <- list()

# Loop through scenarios and calculate differences
for (scenario_name in scenarios_to_compare) {
  if (scenario_name != "BO") {  # We already calculated BO vs RS, skip it to avoid redundancy
    diff_result <- create_diff_between_scenarios("BO", scenario_name, simulation_dt, evaluated_parameter, param_ranges, n_bins)
    all_diff_results[[paste("BO_vs", scenario_name, sep = "_")]] <- diff_result
  }
}



# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned_diff = function(data, params, param_ranges) {
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = diff)
    ) +
    scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0)+
    theme_minimal() +
    labs(x = params[1], y = params[2], fill = "ENLA")
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  
  return(p)
}


p1 = plot_heatmap_binned_diff(all_diff_results$BO_vs_RS, evaluated_parameter, param_ranges)
p2 = plot_heatmap_binned_diff(all_diff_results$BO_vs_PS, evaluated_parameter, param_ranges)
p3 = plot_heatmap_binned_diff(all_diff_results$BO_vs_P2, evaluated_parameter, param_ranges)

gridExtra::grid.arrange(p0, p1, p2, p3, ncol = 2, nrow = 2)



# -------------------------------------------------------------------------






# Plot boxplots -----------------------------------------------------------




interval_size = 0.05
evaluated_parameter = "beta_I_colony"


boxplot_dt = simulation_dt %>%
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





