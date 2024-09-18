
setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")

#
library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)
library(gridExtra)

source("model.R")


# Run simulation ----------------------------------------------------------

time1 <- Sys.time()
output = gillespie_seir(
  
  # Parameter of the taul-leap agorithm
  tau = 0.05,
  # Number of simu_adultlation days
  total_time = 50,
  # Do we induce dispersion ?
  induced_dispersal = F,
  # Induced dispersion mode (deterministic or stochastic)
  dispersal_stochactic = T,
  # Reaction time between 1rst death and induced dispersal
  dispersal_reaction_time = 4,
  # Initial conditions
  initial_number_infected_breeders_A = 3,
  initial_number_breeders_A = 100,
  initial_number_breeders_B = 80,
  initial_number_breeders_C = 20,
  # Transmission rate from exposed individuals and from infectious individuals in a colony
  BETA = 0.02,
  # Time at sea before returning to a colony (non-breeders)
  TIME_AT_SEA_NB = 40
  
)

time2 <- Sys.time()
time2 - time1

plot_seir(output_ = output)