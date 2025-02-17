

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
scenarios = scenarios[c(2,4,5),]

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
  
  output_bank = data.frame()
  
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
      incubation_period = samples[[i, 15]],
      eta = samples[[i, 16]],
      infectious_period = samples[[i, 17]],
      adult_mortality = samples[[i, 18]],
      nestling_mortality = samples[[i, 19]],
      avrg_stay_B_colony = samples[[i, 20]],
      avrg_stay_B_sea = samples[[i, 21]],
      avrg_stay_NB_colony = samples[[i, 22]],
      avrg_stay_NB_sea = samples[[i, 23]],
      theta = samples[[i, 24]],
      psi = samples[[i, 25]],
      hatching_sd = samples[[i, 26]],
      reaching_repro_prob = samples[[i, 27]],
      prob_detection = samples[[i, 28]]
    )
  
    output_bank = rbind(output_bank,
                        cbind(
                         data.frame(
                         nb_adults_equi = output$nb_adults_equi,
                         nb_infected_colonies = output$nb_infected_colonies,
                         infected_X_time = output$infected_X_time),
                         
                         data.frame(output$final_states[1,])
                         )  
                        )
  }
  
  return(output_bank)
}

# Number of samples
nb_samples = 4

# Total number of parameters 
nb_params = length(param_ranges)

# Generate samples
lhs_samples = randomLHS(nb_samples, nb_params)
samples = convert_samples(lhs_samples, param_ranges, nb_params)
colnames(samples) = names(param_ranges)



# Get outputs for all scenarios -------------------------------------------



simulation_dt = data.frame()

for (i in 1:nrow(scenarios)){

  scenario_output = run_simulations(samples,
                                       induced_dispersal_= scenarios[i,"induced_dispersal"],
                                       dispersal_stochastic_ = scenarios[i,"dispersal_stochastic"],
                                       initially_infected_ = scenarios[i,"initially_infected"])

  
  scenario_output = cbind(data.frame(scenario = rep(rownames(scenarios)[i], nrow(samples))),
                          scenario_output)
  
  simulation_dt = rbind(simulation_dt, cbind(samples, scenario_output))
  

}


#save(simulation_dt, file = "simulation_dt_modelv21_100_1.RData")
#save(simulation_dt, file = "simulation_dt_2000_cluster_1.RData")


