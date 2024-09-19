

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

source("param_ranges.R")
source("scenarios.R")

# Convert LHS samples  ----------------------------------------------------

convert_samples = function(lhs_samples, param_ranges, nb_params) {
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

induced_dispersal_ = F
dispersal_stochastic_ = F
initially_infected_ = T

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
      
      # # Parameter samples
      initial_number_infected_breeders_A = samples[[i, 1]],
      initial_number_infected_breeders_B = samples[[i, 2]],
      initial_number_infected_breeders_C = samples[[i, 3]],
      initial_number_breeders_A = samples[[i, 4]],
      initial_number_breeders_B = samples[[i, 5]],
      initial_number_breeders_C = samples[[i, 6]],
      dispersal_reaction_time = samples[[i, 7]],
      dispersal_date = samples[[i, 8]],
      hatching_date = samples[[i, 9]],
      tau = samples[[i, 10]],
      total_time = samples[[i, 11]],
      prop_dispersal = samples[[i, 12]],
      beta_E_colony = samples[[i, 13]],
      beta_I_colony = samples[[i, 14]],
      sigma = samples[[i, 15]],
      eta = samples[[i, 16]],
      gamma = samples[[i, 17]],
      mu_adult = samples[[i, 18]],
      mu_nestling = samples[[i, 19]],
      zeta_to_sea = samples[[i, 20]],
      zeta_to_colony = samples[[i, 21]],
      rho_to_sea = samples[[i, 22]],
      rho_to_colony = samples[[i, 23]],
      psi = samples[[i, 24]],
      hatching_sd = samples[[i, 25]],
      reaching_repro_prob = samples[[i, 26]]
    )
  
    output_bank[i, ] = c(output$nb_adults_equi, output$nb_infected_colonies, output$infected_X_time)
  }
  
  return(output_bank)
}

# Number of samples
nb_samples = 1500

# Total number of parameters 
nb_params = length(param_ranges)

# Generate samples
lhs_samples = randomLHS(nb_samples, nb_params)
samples = convert_samples(lhs_samples, param_ranges, nb_params)
colnames(samples) = names(param_ranges)



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

save(simulation_dt, file = "simulation_dt_1500_1.RData")
#load("simulation_dt.RData")


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

# plot_heatmap(simulation_dt %>% subset(., scenario == plotted_scenario),
#              evaluated_parameter,
#              param_ranges)




