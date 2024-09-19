

setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")


# Library -----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)


# Downloading model -------------------------------------------------------

source("model.R")

# LHS ---------------------------------------------------------------------

library(lhs)
library(ggplot2)
library(dplyr)
library(stringr)

# Parameter ranges --------------------------------------------------------

param_ranges = list(
  
  # Integer parameters
  initial_number_infected_breeders_A = list(range = c(1, 5), distribution = "integer"),
  initial_number_infected_breeders_B = list(range = c(0, 0), distribution = "integer"),
  initial_number_infected_breeders_C = list(range = c(0, 0), distribution = "integer"),
  initial_number_breeders_A = list(range = c(50, 50), distribution = "integer"),
  initial_number_breeders_B = list(range = c(80, 80), distribution = "integer"),
  initial_number_breeders_C = list(range = c(20, 20), distribution = "integer"),
  dispersal_reaction_time =  list(range = c(1, 10), distribution = "integer"),
  dispersal_date = list(range = c(0, 0), distribution = "integer"),
  hatching_date = list(range = c(10, 10), distribution = "integer"),
  
  # Continuous parameters
  tau = list(range = c(0.15, 0.15), distribution = "simple_uniform"),
  total_time = list(range = c(70, 70), distribution = "simple_uniform"),
  prop_dispersal = list(range = c(1, 1), distribution = "simple_uniform"),
  beta_E_colony = list(range = c(0, 0), distribution = "simple_uniform"),
  beta_I_colony = list(range = c(0.01, 0.50), distribution = "logarithmic"),
  sigma = list(range = c(1/1, 1/1), distribution = "logarithmic"),
  eta = list(range = c(0, 0), distribution = "simple_uniform"),
  gamma = list(range = c(1/6, 1/6), distribution = "logarithmic"),
  mu_adult = list(range = c(1/6 * (0.5 / (1 - 0.5)), 1/6 * (0.5 / (1 - 0.5))), distribution = "logarithmic"),
  mu_nestling = list(range = c(1/6 * (0.8 / (1 - 0.8)), 1/6 * (0.8 / (1 - 0.8))), distribution = "logarithmic"),
  zeta_to_sea = list(range = c(1/2, 1/2), distribution = "logarithmic"),
  zeta_to_colony = list(range = c(1/2, 1/2), distribution = "logarithmic"),
  rho_to_sea = list(range = c(1/2, 1/2), distribution = "logarithmic"),
  rho_to_colony = list(range = c(1/40, 1/2), distribution = "logarithmic"),
  psi = list(range = c(1/500, 1/500), distribution = "logarithmic"),
  hatching_sd = list(range = c(3, 3), distribution = "simple_uniform"),
  reaching_repro_prob = list(range = c(0.3, 0.7), distribution = "simple_uniform")
)


# Convert LHS samples  ----------------------------------------------------

convert_samples = function(lhs_samples, param_ranges, nb_params) {
  # Convert integer parameters
  for (i in 1:nb_params) {
    if (param_ranges[[i]]$distribution == "integer"){
      lhs_samples[, i] = qinteger(lhs_samples[, i], param_ranges[[i]]$range[1], param_ranges[[i]]$range[2])
    }
    else if (param_ranges[[i]]$distribution == "logarithmic"){
      log_min <- log(param_ranges[[i]]$range[1])
      log_max <- log(param_ranges[[i]]$range[2])
      lhs_samples[, i] <- exp(qunif(lhs_samples[, i], log_min, log_max))
    }
    else if (param_ranges[[i]]$distribution == "simple_uniform"){
      lhs_samples[, i] = qunif(lhs_samples[, i], param_ranges[[i]]$range[1], param_ranges[[i]]$range[2])
    }
  }
  return(lhs_samples)
}

# Run simulations and store results ---------------------------------------

run_simulations = function(samples,
                            induced_dispersal_, 
                            dispersal_stochastic_,
                             initially_infected_) {
  
  nb_samples = nrow(samples)
  
  output_bank = data.frame(
    nb_adults_equi = numeric(nb_samples),
    nb_infected_colonies = numeric(nb_samples),
    infected_X_time = numeric(nb_samples)
  )
  
  for (i in 1:nb_samples) {
    output = gillespie_seir(
      
      # Scenario
      induced_dispersal = induced_dispersal_,
      dispersal_stochastic  = dispersal_stochastic_,
      initially_infected = initially_infected_,
      
      # Parameter samples
      initial_number_infected_breeders_A = samples[i, 1],
      initial_number_infected_breeders_B = samples[i, 2],
      initial_number_infected_breeders_C = samples[i, 3],
      initial_number_breeders_A = samples[i, 4],
      initial_number_breeders_B = samples[i, 5],
      initial_number_breeders_C = samples[i, 6],
      dispersal_reaction_time = samples[i, 7],
      dispersal_date = samples[i, 8],
      hatching_date = samples[i, 9],
      tau = samples[i, 10],
      total_time = samples[i, 11],
      prop_dispersal = samples[i, 12],
      beta_E_colony = samples[i, 13],
      beta_I_colony = samples[i, 14],
      sigma = samples[i, 15],
      eta = samples[i, 16],
      gamma = samples[i, 17],
      mu_adult = samples[i, 18],
      mu_nestling = samples[i, 19],
      zeta_to_sea = samples[i, 20],
      zeta_to_colony = samples[i, 21],
      rho_to_sea = samples[i, 22],
      rho_to_colony = samples[i, 23],
      psi = samples[i, 24],
      hatching_sd = samples[i, 25],
      reaching_repro_prob = samples[i, 26]
    )
    
    output_bank[i, ] = c(output$nb_adults_equi, output$nb_infected_colonies, output$infected_X_time)
  }
  
  return(output_bank)
}

# Number of samples
nb_samples = 500

# Total number of parameters 
nb_params = length(param_ranges)

# Generate samples
lhs_samples = randomLHS(nb_samples, nb_params)
samples = convert_samples(lhs_samples, param_ranges, nb_params)
colnames(samples) = names(param_ranges)

# Describe the 5 scenario
scenarios = data.frame(
  induced_dispersal = c(F,F,T,T,T),
  dispersal_stochastic = c(F,F,F,F,T),
  initially_infected = c(F,T,F,T,T))

rownames(scenarios) = c("HS","BO","PS","P2","RS")



# Get outputs for all scenarios -------------------------------------------

simulation_dt = data.frame(
  
  scenario = numeric(),
  nb_adults_equi = numeric(),
  nb_infected_colonies = numeric(),
  infected_X_time = numeric()
)

for (i in 1:nrow(scenarios)){
  
  scenario_output = run_simulations(samples, 
                                       induced_dispersal_= scenarios[i,"induced_dispersal"],
                                       dispersal_stochastic_ = scenarios[i,"dispersal_stochastic"],
                                       initially_infected_ = scenarios[i,"initially_infected"])
  
  scenario =  data.frame(scenario = rep(rownames(scenarios)[i], nrow(samples)))
  
  scenario_output = cbind(scenario, scenario_output)
  
  simulation_dt = rbind(simulation_dt, scenario_output)
  
}


simulation_dt = cbind(samples, simulation_dt)

save(simulation_dt, file = "simulation_dt_2.RData")
load("simulation_dt.RData")


# Plots -------------------------------------------------------------------

evaluated_parameter = c("beta_I_colony", "reaching_repro_prob")
plotted_scenario = "BO"

# Plot heatmap ------------------------------------------------------------

plot_heatmap = function(data, params, param_ranges) {
  
  p = ggplot(data) +
    geom_tile(aes(x = get(params[1]), 
                  y = get(params[2]), 
                  fill = nb_adults_equi),
              width = 0.03, height = 0.02) +
    scale_fill_gradient(low = "yellow", high = "blue") +
    theme_minimal() +
    labs(x = params[1],  
         y = params[2],
         fill = "ENLA")
  
  if (param_ranges[[params[1]]][[2]] == "logarithmic"){
    p = p + scale_x_log10()
  }
  if (param_ranges[[params[2]]][[2]] == "logarithmic"){
    p = p + scale_y_log10()
  }
  return(p)
}

plot_heatmap(simulation_dt %>% subset(., scenario == plotted_scenario),
             evaluated_parameter,
             param_ranges)





# Plot binned heatmap -----------------------------------------------------

# Number of bins (same for both dimensions)
n_bins = 6  # You can adjust this number

data = simulation_dt %>% subset(., scenario == plotted_scenario)
params = evaluated_parameter



# Function to create binned data with dynamic parameters and variable block sizes
create_binned_data = function(data, params, param_ranges, n_bins) {
  
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

  
  data = data %>%
    mutate(
      V1 = cut(data$var1, breaks = seq(min_x, max_x, by = block_size_x), include.lowest = TRUE),
      V2 = cut(data$var2, breaks = seq(min_y, max_y, by = block_size_y), include.lowest = TRUE)
    ) %>%
    group_by(V1, V2) %>%
    summarise(nb_adults_equi_avg = mean(nb_adults_equi, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    mutate(
      # Use str_extract to extract the first numeric value from V1 and V2
      x_mid = as.numeric(str_extract(V1, "[-+]?[0-9]*\\.?[0-9]+")) + block_size_x / 2,
      y_mid = as.numeric(str_extract(V2, "[-+]?[0-9]*\\.?[0-9]+")) + block_size_y / 2
    )
  
  # Apply exp transformation if specified
  if (param_ranges[[params[1]]][[2]] == "logarithmic") {
    data$x_mid = exp(data$x_mid)
  } else {
    data$x_mid = data$x_mid
  }
  
  if (param_ranges[[params[2]]][[2]] == "logarithmic") {
    data$y_mid = exp(data$y_mid)
  } else {
    data$y_mid = data$y_mid
  }
  
  return(data)
}

# Function to create a heatmap with dynamic block sizes
plot_heatmap_binned = function(data, params, param_ranges) {
  p = ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = nb_adults_equi_avg)
              ) +
    scale_fill_gradient(low = "yellow", high = "blue") +
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
plot_data_binned = create_binned_data(simulation_dt %>% subset(., scenario == plotted_scenario),
                                      evaluated_parameter,
                                      param_ranges,
                                      n_bins)



# Visualize the heatmap with dynamic block sizes
plot_heatmap_binned(plot_data_binned, evaluated_parameter, param_ranges)



# # Plot boxplots
# 
# 
# interval_size = 0.05
# evaluated_parameter = "beta_I_colony"
# 
# 
# boxplot_dt = simulation_dt %>%
#   mutate(parameter = cut(
#     get(evaluated_parameter),  
#     breaks = seq(
#       param_ranges[[evaluated_parameter]][[1]][[1]],  
#       param_ranges[[evaluated_parameter]][[1]][[2]],
#       by = interval_size
#     ),
#     include.lowest = TRUE
#   )) %>%
#   na.omit()
# 
# 
# ggplot(boxplot_dt, aes(x = parameter, y = nb_adults_equi)) +
#   geom_boxplot(fill = "lightblue", color = "darkblue") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x = evaluated_parameter, y = "ENLA")

