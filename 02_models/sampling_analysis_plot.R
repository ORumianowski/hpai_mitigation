

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

# Parameter ranges --------------------------------------------------------

param_ranges = list(
  
  # Integer parameters
  initial_number_infected_breeders_A = list(range = c(1, 1), distribution = "integer"),
  initial_number_infected_breeders_B = list(range = c(0, 0), distribution = "integer"),
  initial_number_infected_breeders_C = list(range = c(0, 0), distribution = "integer"),
  initial_number_breeders_A = list(range = c(50, 50), distribution = "integer"),
  initial_number_breeders_B = list(range = c(80, 80), distribution = "integer"),
  initial_number_breeders_C = list(range = c(20, 20), distribution = "integer"),
  dispersal_reaction_time =  list(range = c(5, 5), distribution = "integer"),
  dispersal_date = list(range = c(0, 0), distribution = "integer"),
  hatching_date = list(range = c(10, 10), distribution = "integer"),
  
  # Continuous parameters
  tau = list(range = c(0.25, 0.25), distribution = "simple_uniform"),
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
  reaching_repro_prob = list(range = c(0.2, 0.8), distribution = "simple_uniform")
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

run_simulations = function(samples, nb_samples,
                            induced_dispersal_, 
                            dispersal_stochastic_) {
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
nb_samples = 100

# Total number of parameters 
nb_params = length(param_ranges)

# Generate samples
lhs_samples = randomLHS(nb_samples, nb_params)
samples = convert_samples(lhs_samples, param_ranges, nb_params)
colnames(samples) = names(param_ranges)


scenarios = list(
  
  HS = c(),
  BO = ,
  RS = ,
  PA = ,
  P2 = ,
)


simulation_results = run_simulations(samples, nb_samples,
                                      induced_dispersal_= T,
                                      dispersal_stochastic_ = T)


simulation_dt = cbind(samples, simulation_results)


# Plot heatmap ------------------------------------------------------------

evaluated_parameter = c("rho_to_colony", "reaching_repro_prob")

plot_heatmap = function(data, params) {
  ggplot(data) +
    geom_tile(aes(x = get(params[1]), 
                  y = get(params[2]), 
                  fill = nb_adults_equi), width = 0.01, height = 0.01) +
    scale_fill_gradient(low = "yellow", high = "blue") +
    theme_minimal() +
    labs(x = params[1],  
         y = params[2],
         fill = "ENLA")
}

plot_heatmap(simulation_dt, evaluated_parameter)


# Plot binned heatmap -----------------------------------------------------

evaluated_parameter = evaluated_parameter #c("rho_to_colony", "reaching_repro_prob")

# Block size
block_size = 0.05

# Function to create binned data with dynamic parameters
create_binned_data = function(data, params, block_size) {
  data %>%
    mutate(
      V1 = cut(get(params[1]), breaks = seq(min(get(params[1])), max(get(params[1])), by = block_size), include.lowest = TRUE),
      V2 = cut(get(params[2]), breaks = seq(min(get(params[2])), max(get(params[2])), by = block_size), include.lowest = TRUE)
    ) %>%
    group_by(V1, V2) %>%
    summarise(nb_adults_equi_avg = mean(nb_adults_equi, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      x_mid = as.numeric(sub("\\((.+),.*", "\\1", V1)) + block_size / 2,
      y_mid = as.numeric(sub("\\((.+),.*", "\\1", V2)) + block_size / 2
    )
}

# Function to create a heatmap with dynamic parameters
plot_heatmap_binned = function(data, params) {
  ggplot(data) +
    geom_tile(aes(x = x_mid, y = y_mid, fill = nb_adults_equi_avg), width = block_size, height = block_size) +
    scale_fill_gradient(low = "yellow", high = "blue") +
    theme_minimal() +
    labs(x = params[1], y = params[2], fill = "ENLA")
}

# Create binned data with dynamic parameters
plot_data_binned = create_binned_data(simulation_dt, evaluated_parameter, block_size)

# Visualize the heatmap with dynamic parameters
plot_heatmap_binned(plot_data_binned, evaluated_parameter)

# Plot boxplots

# Définir la taille des intervalles et le nom du paramètre
interval_size = 0.05
evaluated_parameter = "reaching_repro_prob"

# Préparer les données pour le boxplot avec des noms de variables dynamiques
boxplot_dt = simulation_dt %>%
  mutate(parameter = cut(
    get(evaluated_parameter),  # Utiliser get() pour accéder dynamiquement à la variable
    breaks = seq(
      param_ranges[[evaluated_parameter]][1],  # Accéder dynamiquement à la plage de paramètres
      param_ranges[[evaluated_parameter]][2],
      by = interval_size
    ),
    include.lowest = TRUE
  )) %>%
  na.omit()

# Créer le boxplot
ggplot(boxplot_dt, aes(x = parameter, y = nb_adults_equi)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = evaluated_parameter, y = "ENLA")

