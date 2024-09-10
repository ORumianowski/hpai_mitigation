

##

# Meilleur dicernement de la prob recrutement

# sortir les paramters

# regroupe le parametre de probabilite de recrutement

# Should I add failures unrelated to the disease?


# color in relation the position compared to baseline outbreak

##






# Library -----------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)


# Model description -------------------------------------------------------

# Each compartment in the model is described by

# the epidemiological status : Susceptible (S), Exposed (E), Infectious (I), Recovered (R)
# the reproductive status : Breeders (B), Non-Breeders (NB), Nestlings (N)
# the localisation : 
# colony A (A), colony (B),
# "sea associted to colony A (sea_a), sea associted to colony B (sea_B),
# sea Non-Breeders (sea_NB)


# Event rates -------------------------------------------------------------


calculate_rates = function(  beta_E_colony, beta_I_colony,
                             sigma,eta, gamma, mu_adult, mu_nestling,
                             zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                             # A
                             S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                             S_a_B, E_a_B, I_a_B, R_a_B,  D_a_B, 
                             S_sea_a_B, E_sea_a_B, I_sea_a_B, R_sea_a_B, D_sea_a_B,
                             S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                             # B
                             S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                             S_b_B, E_b_B, I_b_B, R_b_B, D_b_B,
                             S_sea_b_B, E_sea_b_B, I_sea_b_B, R_sea_b_B, D_sea_b_B,
                             S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                             # C
                             S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                             S_c_B, E_c_B, I_c_B, R_c_B, D_c_B,
                             S_sea_c_B, E_sea_c_B, I_sea_c_B, R_sea_c_B, D_sea_c_B,
                             S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                             # Non-breeders at sea
                             S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB
                             
                             )
                            
                             {
  
  rates = c(
    
  # SEIR transitions
  
  ## Nestlings
  ### Colony A
  "S_a_N_to_E_a_N" = beta_E_colony * S_a_N * (E_a_B+E_a_NB+E_a_N) + 
                     beta_I_colony * S_a_N * (I_a_B+I_a_NB+I_a_N),
  "E_a_N_to_S_a_N" = eta * E_a_N,
  "E_a_N_to_I_a_N" = sigma * E_a_N,
  "I_a_N_to_R_a_N" = gamma * I_a_N,
  "I_a_N_to_D_a_N" = mu_nestling * I_a_N,
  ### Colony B
  "S_b_N_to_E_b_N" = beta_E_colony * S_b_N * (E_b_B+E_b_NB+E_b_N) +
                     beta_I_colony * S_b_N * (I_b_B+I_b_NB+I_b_N),
  "E_b_N_to_S_b_N" = eta * E_b_N,
  "E_b_N_to_I_b_N" = sigma * E_b_N,
  "I_b_N_to_R_b_N" = gamma * I_b_N,
  "I_b_N_to_D_b_N" = mu_nestling * I_b_N,
  ### Colony C
  "S_c_N_to_E_c_N" = beta_E_colony * S_c_N * (E_c_B+E_c_NB+E_c_N) +
                     beta_I_colony * S_c_N * (I_c_B+I_c_NB+I_c_N),
  "E_c_N_to_S_c_N" = eta * E_c_N,
  "E_c_N_to_I_c_N" = sigma * E_c_N,
  "I_c_N_to_R_c_N" = gamma * I_c_N,
  "I_c_N_to_D_c_N" = mu_nestling * I_c_N,
  
  ## Non-Breeders
  ### Colony A
  "S_a_NB_to_E_a_NB" = beta_E_colony * S_a_NB * (E_a_B+E_a_NB+E_a_N) +
                       beta_I_colony * S_a_NB * (I_a_B+I_a_NB+I_a_N),
  "E_a_NB_to_S_a_NB" = eta * E_a_NB,
  "E_a_NB_to_I_a_NB" = sigma * E_a_NB,
  "I_a_NB_to_R_a_NB" = gamma * I_a_NB,
  "I_a_NB_to_D_a_NB" = mu_adult * I_a_NB,
  ### Colony B
  "S_b_NB_to_E_b_NB" = beta_E_colony * S_b_NB * (E_b_B+E_b_NB+E_b_N) + 
                       beta_I_colony * S_b_NB * (I_b_B+I_b_NB+I_b_N),
  "E_b_NB_to_S_b_NB" = eta * E_b_NB,
  "E_b_NB_to_I_b_NB" = sigma * E_b_NB,
  "I_b_NB_to_R_b_NB" = gamma * I_b_NB,
  "I_b_NB_to_D_b_NB" = mu_adult * I_b_NB,
  ### Colony C
  "S_c_NB_to_E_c_NB" = beta_E_colony * S_c_NB * (E_c_B+E_c_NB+E_c_N) + 
                       beta_I_colony * S_c_NB * (I_c_B+I_c_NB+I_c_N),
  "E_c_NB_to_S_c_NB" = eta * E_c_NB,
  "E_c_NB_to_I_c_NB" = sigma * E_c_NB,
  "I_c_NB_to_R_c_NB" = gamma * I_c_NB,
  "I_c_NB_to_D_c_NB" = mu_adult * I_c_NB,
  ### Sea 
  "S_sea_NB_to_E_sea_NB" = 0,
  "E_sea_NB_to_S_sea_NB" = eta * E_sea_NB,
  "E_sea_NB_to_I_sea_NB" = sigma * E_sea_NB,
  "I_sea_NB_to_R_sea_NB" = gamma * I_sea_NB,
  "I_sea_NB_to_D_sea_NB" = mu_adult * I_sea_NB,

  
  ## Breeders
  ### Colony A
  "S_a_B_to_E_a_B" = beta_E_colony * S_a_B * (E_a_B+E_a_NB+E_a_N) + 
                 beta_I_colony * S_a_B * (I_a_B+I_a_NB+I_a_N),
  "E_a_B_to_S_a_B" = eta * E_a_B,
  "E_a_B_to_I_a_B" = sigma * E_a_B,
  "I_a_B_to_R_a_B" = gamma * I_a_B,
  "I_a_B_to_D_a_B" = mu_adult * I_a_B,
  ### Sea A
  "S_sea_a_B_to_E_sea_a_B" = 0,
  "E_sea_a_B_to_S_sea_a_B" = eta * E_sea_a_B,
  "E_sea_a_B_to_I_sea_a_B" = sigma * E_sea_a_B,
  "I_sea_a_B_to_R_sea_a_B" = gamma * I_sea_a_B,
  "I_sea_a_B_to_D_sea_a_B" = mu_adult * I_sea_a_B,
  ### Colony B
  "S_b_B_to_E_b_B" = beta_E_colony * S_b_B * (E_b_B+E_b_NB+E_b_N) + 
                 beta_I_colony * S_b_B * (I_b_B+I_b_NB+I_b_N),
  "E_b_B_to_S_b_B" = eta * E_b_B,
  "E_b_B_to_I_b_B" = sigma * E_b_B,
  "I_b_B_to_R_b_B" = gamma * I_b_B,
  "I_b_B_to_D_b_B" = mu_adult * I_b_B,
  ### Sea B
  "S_sea_b_B_to_E_sea_b_B" = 0,
  "E_sea_b_B_to_S_sea_b_B" = eta * E_sea_b_B,
  "E_sea_b_B_to_I_sea_b_B" = sigma * E_sea_b_B,
  "I_sea_b_B_to_R_sea_b_B" = gamma * I_sea_b_B,
  "I_sea_b_B_to_D_sea_b_B" = mu_adult * I_sea_b_B,
  ### Colony C
  "S_c_B_to_E_c_B" = beta_E_colony * S_c_B * (E_c_B+E_c_NB+E_c_N) + 
                 beta_I_colony * S_c_B * (I_c_B+I_c_NB+I_c_N),
  "E_c_B_to_S_c_B" = eta * E_c_B,
  "E_c_B_to_I_c_B" = sigma * E_c_B,
  "I_c_B_to_R_c_B" = gamma * I_c_B,
  "I_c_B_to_D_c_B" = mu_adult * I_c_B,
  ### Sea C
  "S_sea_c_B_to_E_sea_c_B" = 0,
  "E_sea_c_B_to_S_sea_c_B" = eta * E_sea_c_B,
  "E_sea_c_B_to_I_sea_c_B" = sigma * E_sea_c_B,
  "I_sea_c_B_to_R_sea_c_B" = gamma * I_sea_c_B,
  "I_sea_c_B_to_D_sea_c_B" = mu_adult * I_sea_c_B,

  
  
  # Mobility
  ## Non-Breeders
  ### In A
  ### From colony A to sea 
  "S_a_NB_to_S_sea_NB" = rho_to_sea * S_a_NB,
  "E_a_NB_to_E_sea_NB" = rho_to_sea * E_a_NB,
  "I_a_NB_to_I_sea_NB" = rho_to_sea * I_a_NB,
  "R_a_NB_to_R_sea_NB" = rho_to_sea * R_a_NB,
  ### From sea to colony A (conspecific attraction)
  "S_sea_NB_to_S_a_NB" = rho_to_colony * S_sea_NB * ( (S_a_B+E_a_B+I_a_B+R_a_B)
                                                     /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                       S_b_B+E_b_B+I_b_B+R_b_B+
                                                       S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                     ),
  
  "E_sea_NB_to_E_a_NB" = rho_to_colony * E_sea_NB * ( (S_a_B+E_a_B+I_a_B+R_a_B)
                                                     /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                       S_b_B+E_b_B+I_b_B+R_b_B+
                                                       S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                     ),
  "I_sea_NB_to_I_a_NB" = rho_to_colony * I_sea_NB * ( (S_a_B+E_a_B+I_a_B+R_a_B)
                                                     /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                       S_b_B+E_b_B+I_b_B+R_b_B+
                                                       S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                     ),
  "R_sea_NB_to_R_a_NB" = rho_to_colony * R_sea_NB * ( (S_a_B+E_a_B+I_a_B+R_a_B)
                                                      /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                        S_b_B+E_b_B+I_b_B+R_b_B+
                                                        S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                      ),
  ### In B
  ### From colony B to sea
  "S_b_NB_to_S_sea_NB" = rho_to_sea * S_b_NB,
  "E_b_NB_to_E_sea_NB" = rho_to_sea * E_b_NB,
  "I_b_NB_to_I_sea_NB" = rho_to_sea * I_b_NB,
  "R_b_NB_to_R_sea_NB" = rho_to_sea * R_b_NB,
  ### From sea to colony B
  "S_sea_NB_to_S_b_NB" = rho_to_colony * S_sea_NB * ( (S_b_B+E_b_B+I_b_B+R_b_B)
                                                     /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                       S_b_B+E_b_B+I_b_B+R_b_B+
                                                       S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                     ),
  "E_sea_NB_to_E_b_NB" = rho_to_colony * E_sea_NB * ( (S_b_B+E_b_B+I_b_B+R_b_B)
                                                     /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                       S_b_B+E_b_B+I_b_B+R_b_B+
                                                       S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                     ),
  "I_sea_NB_to_I_b_NB" = rho_to_colony * I_sea_NB * ( (S_b_B+E_b_B+I_b_B+R_b_B)
                                                     /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                       S_b_B+E_b_B+I_b_B+R_b_B+
                                                       S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                     ),
  "R_sea_NB_to_R_b_NB" = rho_to_colony * R_sea_NB * ( (S_b_B+E_b_B+I_b_B+R_b_B)
                                                     /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                       S_b_B+E_b_B+I_b_B+R_b_B+
                                                       S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                     ),
  ### In C
  ### From colony C to sea
  "S_c_NB_to_S_sea_NB" = rho_to_sea * S_c_NB,
  "E_c_NB_to_E_sea_NB" = rho_to_sea * E_c_NB,
  "I_c_NB_to_I_sea_NB" = rho_to_sea * I_c_NB,
  "R_c_NB_to_R_sea_NB" = rho_to_sea * R_c_NB,
  ### From sea to colony C
  "S_sea_NB_to_S_c_NB" = rho_to_colony * S_sea_NB * ( (S_c_B+E_c_B+I_c_B+R_c_B)
                                                      /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                        S_b_B+E_b_B+I_b_B+R_b_B+
                                                        S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                      ),
  "E_sea_NB_to_E_c_NB" = rho_to_colony * E_sea_NB * ( (S_c_B+E_c_B+I_c_B+R_c_B)
                                                      /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                        S_b_B+E_b_B+I_b_B+R_b_B+
                                                        S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                      ),
  "I_sea_NB_to_I_c_NB" = rho_to_colony * I_sea_NB * ( (S_c_B+E_c_B+I_c_B+R_c_B)
                                                      /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                        S_b_B+E_b_B+I_b_B+R_b_B+
                                                        S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                      ),
  "R_sea_NB_to_R_c_NB" = rho_to_colony * R_sea_NB * ( (S_c_B+E_c_B+I_c_B+R_c_B)
                                                      /(S_a_B+E_a_B+I_a_B+R_a_B+
                                                        S_b_B+E_b_B+I_b_B+R_b_B+
                                                        S_c_B+E_c_B+I_c_B+R_c_B+1)
                                                      ),
  
  ## Breeders
  ### In A
  ### From colony A to sea A
  "S_a_B_to_S_sea_a_B" = zeta_to_sea * S_a_B,
  "E_a_B_to_E_sea_a_B" = zeta_to_sea * E_a_B,
  "I_a_B_to_I_sea_a_B" = zeta_to_sea * I_a_B,
  "R_a_B_to_R_sea_a_B" = zeta_to_sea * R_a_B, 
  ### From sea A to colony A
  "S_sea_a_B_to_S_a_B" = zeta_to_colony * S_sea_a_B,
  "E_sea_a_B_to_E_a_B" = zeta_to_colony * E_sea_a_B,
  "I_sea_a_B_to_I_a_B" = zeta_to_colony * I_sea_a_B,
  "R_sea_a_B_to_R_a_B" = zeta_to_colony * R_sea_a_B,
  ### In B
  ### From colony B to sea B
  "S_b_B_to_S_sea_b_B" = zeta_to_sea * S_b_B,
  "E_b_B_to_E_sea_b_B" = zeta_to_sea * E_b_B,
  "I_b_B_to_I_sea_b_B" = zeta_to_sea * I_b_B,
  "R_b_B_to_R_sea_b_B" = zeta_to_sea * R_b_B,
  ### From sea B to colony B
  "S_sea_b_B_to_S_b_B" = zeta_to_colony * S_sea_b_B,
  "E_sea_b_B_to_E_b_B" = zeta_to_colony * E_sea_b_B,
  "I_sea_b_B_to_I_b_B" = zeta_to_colony * I_sea_b_B,
  "R_sea_b_B_to_R_b_B" = zeta_to_colony * R_sea_b_B,
  ### In C
  ### From colony B to sea C
  "S_c_B_to_S_sea_c_B" = zeta_to_sea * S_c_B,
  "E_c_B_to_E_sea_c_B" = zeta_to_sea * E_c_B,
  "I_c_B_to_I_sea_c_B" = zeta_to_sea * I_c_B,
  "R_c_B_to_R_sea_c_B" = zeta_to_sea * R_c_B,
  ### From sea B to colony C
  "S_sea_c_B_to_S_c_B" = zeta_to_colony * S_sea_c_B,
  "E_sea_c_B_to_E_c_B" = zeta_to_colony * E_sea_c_B,
  "I_sea_c_B_to_I_c_B" = zeta_to_colony * I_sea_c_B,
  "R_sea_c_B_to_R_c_B" = zeta_to_colony * R_sea_c_B
  
  # Breeders become Non-Breeders
  
  # "S_a_B_to_S_a_NB" = psi * S_a_B,
  # "E_a_B_to_E_a_NB" = psi * E_a_B,
  # "I_a_B_to_I_a_NB" = psi * I_a_B,
  # "R_a_B_to_R_a_NB" = psi * R_a
  
  )
  return(rates)
}


# Gillespie SEIR model function -------------------------------------------

gillespie_seir = function(# Parameter of the taul-leap agorithm
                          tau = 0.25,
                          # Number of simu_adultlation days
                          total_time = 70,
                          # Initial conditions
                          initial_number_infected_breeders_A = 1,
                          initial_number_breeders_A = 50,
                          initial_number_breeders_B = 50,
                          initial_number_breeders_C = 50,
                          # Induced dispersion parameters
                          # Do we induce dispersion ?
                          induced_dispersal = T,
                          # Induced dispersion mode (deterministic or stochastic)
                          dispersal_stochastic = T,
                          # Reaction time between 1rst death and induced dispersal
                          dispersal_reaction_time = 5,
                          ## Proportion of dispersed adults
                          prop_dispersal = 1,
                          ## Date of induced dispersion (if deterministic)
                          dispersal_date = 0,
                          # Epidemiological parameters
                          ## Transmission rate from exposed and infectious individuals in a colony
                          beta_E_colony = 0,
                          beta_I_colony = 0.05,
                          ## Rate of progression from exposed to infectious (inverse of incubation period)
                          sigma = 1/1,
                          ## Rate of progression from exposed to susceptible 
                          eta =  0, 
                          ## Recovery rate (inverse of infectious period)
                          gamma = 1/6,
                          ## Disease-related mortality rate
                          ## Death probability = mu / (mu + gamma)
                          ## Adult
                          mu_adult = 1/6 * (0.5/(1-0.5)), # 50% of mortality
                          ## Nestling
                          mu_nestling = 1/6 * (0.8/(1-0.8)), # 80% of mortality
                          # Mobility  parameters
                          ## Transition from colony to the sea (breeders)
                          zeta_to_sea = 1/2,
                          ## Transition from sea to the colony (breeders)
                          zeta_to_colony = 1/2,
                          ## Transition from colony to the sea (non-breeders)
                          rho_to_sea = 1/2,
                          ## Transition from sea to the colony (non-breeders)
                          rho_to_colony = 1/40 , 
                          # Transition from breeder to non-breeder (reproductive failure)
                          # psi = 0,  
                          # Demographic parameters
                          ## Hatching date of the chicks
                          hatching_date = 10
                          
                          ) {
  
  # Parameters --------------------------------------------------------------
  
  # Initial state
  
  ## Nestlings
  ## In colony A
  N = 0                
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_a_N = c(S = initial_susceptible,
                        E = initial_exposed,
                        I = initial_infected,
                        R = initial_recovered,
                        D = initial_dead)
  ## In colony B
  N = 0                
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_b_N = c(S = initial_susceptible,
                        E = initial_exposed,
                        I = initial_infected,
                        R = initial_recovered,
                        D = initial_dead)
  
  ## In colony C
  N = 0                
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_c_N = c(S = initial_susceptible,
                        E = initial_exposed,
                        I = initial_infected,
                        R = initial_recovered,
                        D = initial_dead)
  
  
  ## Breeders
  ## In colony A
  N = initial_number_breeders_A/2                
  initial_infected = initial_number_infected_breeders_A
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_a_B = c(S = initial_susceptible,
                      E = initial_exposed,
                      I = initial_infected,
                      R = initial_recovered,
                      D = initial_dead)
  
  ## At sea A
  N = initial_number_breeders_A/2                   
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_sea_a_B = c(S = initial_susceptible,
                          E = initial_exposed,
                          I = initial_infected,
                          R = initial_recovered,
                          D = initial_dead)
  
  
  ## In colony B
  N = initial_number_breeders_B/2                  
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_b_B <- c(S = initial_susceptible,
                       E = initial_exposed,
                       I = initial_infected,
                       R = initial_recovered,
                       D = initial_dead)
  
  ## At sea B
  N = initial_number_breeders_B/2                
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_sea_b_B = c(S = initial_susceptible,
                          E = initial_exposed,
                          I = initial_infected,
                          R = initial_recovered,
                          D = initial_dead)
  
  ## In colony C
  N = initial_number_breeders_C/2               
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_c_B <- c(S = initial_susceptible,
                       E = initial_exposed,
                       I = initial_infected,
                       R = initial_recovered,
                       D = initial_dead)
  
  ## At sea C
  N = initial_number_breeders_C/2                
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_sea_c_B = c(S = initial_susceptible,
                          E = initial_exposed,
                          I = initial_infected,
                          R = initial_recovered,
                          D = initial_dead)
  
  ## Non-Breeders
  ## In colony A
  N = 0                  
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_a_NB = c(S = initial_susceptible,
                         E = initial_exposed,
                         I = initial_infected,
                         R = initial_recovered,
                         D = initial_dead)
  
  ## In colony B
  N = 0                  
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_b_NB <- c(S = initial_susceptible,
                          E = initial_exposed,
                          I = initial_infected,
                          R = initial_recovered,
                          D = initial_dead)
  
  ## In colony C
  N = 0                  
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_c_NB <- c(S = initial_susceptible,
                          E = initial_exposed,
                          I = initial_infected,
                          R = initial_recovered,
                          D = initial_dead)
  
  ## At sea 
  N = 0                  
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_sea_NB = c(S = initial_susceptible,
                             E = initial_exposed,
                             I = initial_infected,
                             R = initial_recovered,
                             D = initial_dead)
  
  
  initial_state = matrix(data = c(initial_state_a_N,
                                  initial_state_a_B,
                                  initial_state_sea_a_B,
                                  initial_state_a_NB,
                                  
                                  initial_state_b_N,
                                  initial_state_b_B,
                                  initial_state_sea_b_B,
                                  initial_state_b_NB,
                                  
                                  initial_state_c_N,
                                  initial_state_c_B,
                                  initial_state_sea_c_B,
                                  initial_state_c_NB,
  
                                  initial_state_sea_NB
                                  ), 
                         nrow = 13, ncol = 5, 
                         byrow = T)
  
  
  # Initialization
  times = c(0)
  states = array(dim = c(13,5,1), data = initial_state)
  already_dispersed = F
  simulated_dispersal_date = NA
  already_hatched = F
  first_death = F
  first_death_date = NA
  
  
  # Next event
  while (times[length(times)] < total_time) {
    
    S_a_N = states[1, 1, dim(states)[3]]
    E_a_N = states[1, 2, dim(states)[3]]
    I_a_N = states[1, 3, dim(states)[3]]
    R_a_N = states[1, 4, dim(states)[3]]
    D_a_N = states[1, 5, dim(states)[3]]
    
    S_a_B = states[2, 1, dim(states)[3]]
    E_a_B= states[2, 2, dim(states)[3]]
    I_a_B = states[2, 3, dim(states)[3]]
    R_a_B = states[2, 4, dim(states)[3]]
    D_a_B = states[2, 5, dim(states)[3]]
    
    S_sea_a_B = states[3, 1, dim(states)[3]]
    E_sea_a_B = states[3, 2, dim(states)[3]]
    I_sea_a_B = states[3, 3, dim(states)[3]]
    R_sea_a_B = states[3, 4, dim(states)[3]]
    D_sea_a_B = states[3, 5, dim(states)[3]]
  
    S_a_NB = states[4, 1, dim(states)[3]]
    E_a_NB = states[4, 2, dim(states)[3]]
    I_a_NB = states[4, 3, dim(states)[3]]
    R_a_NB = states[4, 4, dim(states)[3]]
    D_a_NB = states[4, 5, dim(states)[3]]
    
    S_b_N = states[5, 1, dim(states)[3]]
    E_b_N = states[5, 2, dim(states)[3]]
    I_b_N = states[5, 3, dim(states)[3]]
    R_b_N = states[5, 4, dim(states)[3]]
    D_b_N = states[5, 5, dim(states)[3]]
    
    S_b_B = states[6, 1, dim(states)[3]]
    E_b_B = states[6, 2, dim(states)[3]]
    I_b_B = states[6, 3, dim(states)[3]]
    R_b_B = states[6, 4, dim(states)[3]]
    D_b_B = states[6, 5, dim(states)[3]]
    
    S_sea_b_B = states[7, 1, dim(states)[3]]
    E_sea_b_B = states[7, 2, dim(states)[3]]
    I_sea_b_B = states[7, 3, dim(states)[3]]
    R_sea_b_B = states[7, 4, dim(states)[3]]
    D_sea_b_B = states[7, 5, dim(states)[3]]
    
    S_b_NB = states[8, 1, dim(states)[3]]
    E_b_NB = states[8, 2, dim(states)[3]]
    I_b_NB = states[8, 3, dim(states)[3]]
    R_b_NB = states[8, 4, dim(states)[3]]
    D_b_NB = states[8, 5, dim(states)[3]]
    
    S_c_N = states[9, 1, dim(states)[3]]
    E_c_N = states[9, 2, dim(states)[3]]
    I_c_N = states[9, 3, dim(states)[3]]
    R_c_N = states[9, 4, dim(states)[3]]
    D_c_N = states[9, 5, dim(states)[3]]
    
    S_c_B = states[10, 1, dim(states)[3]]
    E_c_B = states[10, 2, dim(states)[3]]
    I_c_B = states[10, 3, dim(states)[3]]
    R_c_B = states[10, 4, dim(states)[3]]
    D_c_B = states[10, 5, dim(states)[3]]
    
    S_sea_c_B = states[11, 1, dim(states)[3]]
    E_sea_c_B = states[11, 2, dim(states)[3]]
    I_sea_c_B = states[11, 3, dim(states)[3]]
    R_sea_c_B = states[11, 4, dim(states)[3]]
    D_sea_c_B = states[11, 5, dim(states)[3]]
    
    S_c_NB = states[12, 1, dim(states)[3]]
    E_c_NB = states[12, 2, dim(states)[3]]
    I_c_NB = states[12, 3, dim(states)[3]]
    R_c_NB = states[12, 4, dim(states)[3]]
    D_c_NB = states[12, 5, dim(states)[3]]
    
    S_sea_NB = states[13, 1, dim(states)[3]]
    E_sea_NB = states[13, 2, dim(states)[3]]
    I_sea_NB = states[13, 3, dim(states)[3]]
    R_sea_NB = states[13, 4, dim(states)[3]]
    D_sea_NB = states[13, 5, dim(states)[3]]
    

    
    # Rates of each possible event
    rates = calculate_rates(beta_E_colony, beta_I_colony,
                            sigma,eta, gamma, mu_adult, mu_nestling,
                            zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                            # A
                            S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                            S_a_B, E_a_B, I_a_B, R_a_B,  D_a_B, 
                            S_sea_a_B, E_sea_a_B, I_sea_a_B, R_sea_a_B, D_sea_a_B,
                            S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                            # B
                            S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                            S_b_B, E_b_B, I_b_B, R_b_B, D_b_B,
                            S_sea_b_B, E_sea_b_B, I_sea_b_B, R_sea_b_B, D_sea_b_B,
                            S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                            # C
                            S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                            S_c_B, E_c_B, I_c_B, R_c_B, D_c_B,
                            S_sea_c_B, E_sea_c_B, I_sea_c_B, R_sea_c_B, D_sea_c_B,
                            S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                            # Non-breeders at sea
                            S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB)

    total_rate = sum(rates)
    
    if (total_rate == 0) {
      break
    }
    
    time_step <- min(tau, (total_time - times[length(times)]))
    next_time = times[length(times)] + time_step

    # Hatching
    if (next_time > hatching_date & !already_hatched) { 
      
      S_a_N = round((S_a_B + E_a_B + I_a_B + R_a_B + S_sea_a_B + E_sea_a_B + I_sea_a_B + R_sea_a_B)/2)
      S_b_N = round((S_b_B + E_b_B + I_b_B + R_b_B + S_sea_b_B + E_sea_b_B + I_sea_b_B + R_sea_b_B)/2)
      S_c_N = round((S_c_B + E_c_B + I_c_B + R_c_B + S_sea_c_B + E_sea_c_B + I_sea_c_B + R_sea_c_B)/2)
      
      already_hatched = T
      
      new_state = matrix(data = c(S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                  S_a_B, E_a_B, I_a_B, R_a_B,  D_a_B, 
                                  S_sea_a_B, E_sea_a_B, I_sea_a_B, R_sea_a_B, D_sea_a_B,
                                  S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                  
                                  S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                  S_b_B, E_b_B, I_b_B, R_b_B, D_b_B,
                                  S_sea_b_B, E_sea_b_B, I_sea_b_B, R_sea_b_B, D_sea_b_B,
                                  S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                  
                                  S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                  S_c_B, E_c_B, I_c_B, R_c_B, D_c_B,
                                  S_sea_c_B, E_sea_c_B, I_sea_c_B, R_sea_c_B, D_sea_c_B,
                                  S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                                  
                                  S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB
                                  ),
                         nrow = 13, ncol = 5, 
                         byrow = T)
      
      states = abind(states, new_state)
      times = c(times, hatching_date)
      
      
      # Rates of each possible event
      rates =     rates = calculate_rates(beta_E_colony, beta_I_colony,
                                          sigma,eta, gamma, mu_adult, mu_nestling,
                                          zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                                          # A
                                          S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                          S_a_B, E_a_B, I_a_B, R_a_B,  D_a_B, 
                                          S_sea_a_B, E_sea_a_B, I_sea_a_B, R_sea_a_B, D_sea_a_B,
                                          S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                          # B
                                          S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                          S_b_B, E_b_B, I_b_B, R_b_B, D_b_B,
                                          S_sea_b_B, E_sea_b_B, I_sea_b_B, R_sea_b_B, D_sea_b_B,
                                          S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                          # C
                                          S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                          S_c_B, E_c_B, I_c_B, R_c_B, D_c_B,
                                          S_sea_c_B, E_sea_c_B, I_sea_c_B, R_sea_c_B, D_sea_c_B,
                                          S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                                          # Non-breeders at sea
                                          S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB)
      
      total_rate = sum(rates)
      
      if (total_rate == 0) {
        break
      }
      
      time_step = rexp(1, total_rate)
      next_time = times[length(times)] + time_step
      
    }
    
    
    
    # Determining the date of first infection
    if (!first_death & D_a_B > 0){
      first_death_date = times[length(times)] 
      first_death = T
    }
    
    # Induction of dispersion
    ## Has the induced dispersal strategy been triggered?
    if (induced_dispersal){ 
      ## Has the dispersal occured?
      if (!already_dispersed){
        # Is it time to induce dispersion, according to the stochastic case and the deterministic case?
        if ((dispersal_stochastic & first_death & next_time > first_death_date + dispersal_reaction_time) 
            |
            (!dispersal_stochastic & next_time > dispersal_date)
        ){
          
          # Dispersal of breeders
          
          # Number of adults in  A
          N_a_B = S_a_B + E_a_B+ I_a_B + R_a_B + S_sea_a_B + E_sea_a_B + I_sea_a_B + R_sea_a_B
          # Number of adults A who are dispersed
          N_disp_a_B = round(N_a_B * prop_dispersal)
          # Distribution of dispersed adults by epidemiological status
          disp_a_B = sample(c(rep("S_a_B", S_a_B), rep("E_a_B", E_a_B),rep("I_a_B", I_a_B),rep("R_a_B", R_a_B),
                            rep("S_sea_a_B", S_sea_a_B), rep("E_sea_a_B", E_sea_a_B),rep("I_sea_a_B", I_sea_a_B),rep("R_sea_a_B", R_sea_a_B)),
                          size = N_disp_a_B, 
                          replace = F) %>% 
            factor(., levels = c("S_a_B","E_a_B","I_a_B","R_a_B",
                                 "S_sea_a_B","E_sea_a_B","I_sea_a_B","R_sea_a_B")) %>% 
            table()
          
          # Number of susceptible adults who are dispersed from A
          disp_S_a_B = disp_a_B["S_a_B"]
          disp_S_sea_a_B = disp_a_B["S_sea_a_B"]
          # Update of the number of susceptible adults
          S_a_B = S_a_B - disp_S_a_B
          S_sea_a_B = S_sea_a_B - disp_S_sea_a_B
          S_sea_NB = S_sea_NB + disp_S_a_B + disp_S_sea_a_B
          
          # Number of exposed  adults who are dispersed from A
          disp_E_a_B= disp_a_B["E_a_B"]
          disp_E_sea_a_B = disp_a_B["E_sea_a_B"]
          # Update of the number of exposed adults
          E_a_B= E_a_B- disp_E_a_B
          E_sea_a_B = E_sea_a_B - disp_E_sea_a_B
          E_sea_NB = E_sea_NB + disp_E_a_B+ disp_E_sea_a_B
          
          # Number of infectious adults who are dispersed from A
          disp_I_a_B = disp_a_B["I_a_B"]
          disp_I_sea_a_B = disp_a_B["I_sea_a_B"]
          # Update of the number of infectious adults
          I_a_B = I_a_B - disp_I_a_B
          I_sea_a_B = I_sea_a_B - disp_I_sea_a_B
          I_sea_NB = I_sea_NB + disp_I_a_B + disp_I_sea_a_B
          
          # Number of recovered adults who are dispersed from A
          disp_R_a_B = disp_a_B["R_a_B"]
          disp_R_sea_a_B = disp_a_B["R_sea_a_B"]
          # Update of the number of recovered adults
          R_a_B = R_a_B - disp_R_a_B
          R_sea_a_B = R_sea_a_B - disp_R_sea_a_B
          R_sea_NB = R_sea_NB + disp_R_a_B + disp_R_sea_a_B
          
          # Death of nestlings
          
          D_a_N = D_a_N + (S_a_N + E_a_N + I_a_N + R_a_N)
          S_a_N = 0
          E_a_N = 0
          I_a_N = 0
          R_a_N = 0
          
          
          already_dispersed = T
          
          new_state = matrix(data = c(S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                      S_a_B, E_a_B, I_a_B, R_a_B,  D_a_B, 
                                      S_sea_a_B, E_sea_a_B, I_sea_a_B, R_sea_a_B, D_sea_a_B,
                                      S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                      
                                      S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                      S_b_B, E_b_B, I_b_B, R_b_B, D_b_B,
                                      S_sea_b_B, E_sea_b_B, I_sea_b_B, R_sea_b_B, D_sea_b_B,
                                      S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                      
                                      S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                      S_c_B, E_c_B, I_c_B, R_c_B, D_c_B,
                                      S_sea_c_B, E_sea_c_B, I_sea_c_B, R_sea_c_B, D_sea_c_B,
                                      S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                                      
                                      S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB
          ),
          nrow = 13, ncol = 5, 
          byrow = T)
          
          states = abind(states, new_state)
          times = c(times, next_time)
          
          simulated_dispersal_date = next_time
          
          
          # Rates of each possible event
          rates =     rates = calculate_rates(beta_E_colony, beta_I_colony,
                                              sigma,eta, gamma, mu_adult, mu_nestling,
                                              zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                                              # A
                                              S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                              S_a_B, E_a_B, I_a_B, R_a_B,  D_a_B, 
                                              S_sea_a_B, E_sea_a_B, I_sea_a_B, R_sea_a_B, D_sea_a_B,
                                              S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                              # B
                                              S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                              S_b_B, E_b_B, I_b_B, R_b_B, D_b_B,
                                              S_sea_b_B, E_sea_b_B, I_sea_b_B, R_sea_b_B, D_sea_b_B,
                                              S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                              # C
                                              S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                              S_c_B, E_c_B, I_c_B, R_c_B, D_c_B,
                                              S_sea_c_B, E_sea_c_B, I_sea_c_B, R_sea_c_B, D_sea_c_B,
                                              S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                                              # Non-breeders at sea
                                              S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB)
          
          total_rate = sum(rates)
          
          if (total_rate == 0) {
            break
          }
          
          time_step = rexp(1, total_rate)
          next_time = times[length(times)] + time_step
          
          
        }
        
      }
    }
    
    # Continuation of the script if there is no special event (dispersal or hatching)
    
    times = c(times, next_time)
    
    nb_events <- rpois(1, total_rate * time_step) 
    whichevent = sample(1:length(rates), nb_events, prob = rates, replace = TRUE)
    
    transitions = factor(whichevent, levels = 1:length(rates)) %>% 
      table() %>% 
      as.matrix() %>% 
      t() %>% 
      as.data.frame()
    
    names_transitions = names(rates)
    names(transitions) <- names_transitions

    transitions_bank = c()
    for (i in 1:length(transitions)){
      if (length(transitions)>0){
        transitions_bank = c(transitions_bank, rep(names_transitions[i], times = transitions[i]))
      }
    }
    
    
    if (length(transitions_bank)>0){
      for (i in 1:length(transitions_bank)){
      
        transition = transitions_bank[i]
      
    # SEIR
    ## Nestling
    ### In A
    if (transition == "S_a_N_to_E_a_N" & S_a_N > 0){
      S_a_N = S_a_N - 1
      E_a_N = E_a_N + 1
    }else if (transition == "E_a_N_to_S_a_N" & E_a_N > 0){
      E_a_N = E_a_N - 1
      S_a_N = S_a_N + 1
    }else if (transition == "E_a_N_to_I_a_N" & E_a_N > 0){
      E_a_N = E_a_N - 1
      I_a_N = I_a_N + 1
    }else if (transition == "I_a_N_to_R_a_N" & I_a_N > 0){
      I_a_N = I_a_N - 1
      R_a_N = R_a_N + 1
    }else if (transition == "I_a_N_to_R_a_N" & I_a_N > 0){
      I_a_N = I_a_N - 1
      R_a_N = R_a_N + 1
    }else if (transition == "I_a_N_to_D_a_N"  & I_a_N > 0){
      # If a nestling dies, both parents become non-breeders.
      I_a_N = I_a_N - 1
      D_a_N = D_a_N + 1
      parent1 = sample(c(rep("S_a_B", S_a_B), rep("E_a_B", E_a_B),rep("I_a_B", I_a_B),rep("R_a_B", R_a_B),
                         rep("S_sea_a_B", S_sea_a_B), rep("E_sea_a_B", E_sea_a_B),rep("I_sea_a_B", I_sea_a_B),rep("R_sea_a_B", R_sea_a_B)),
                       size = 1)
      if (parent1 == "S_a_B"){
        S_a_B = S_a_B - 1
        S_a_NB = S_a_NB + 1
      } else if (parent1 == "E_a_B"){
        E_a_B= E_a_B- 1
        E_a_NB = E_a_NB + 1
      } else if (parent1 == "I_a_B"){
        I_a_B = I_a_B - 1
        I_a_NB = I_a_NB + 1
      } else if (parent1 == "R_a_B"){
        R_a_B = R_a_B - 1
        R_a_NB = R_a_NB + 1
      } else if (parent1 == "S_sea_a_B"){
        S_sea_a_B = S_sea_a_B - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent1 == "E_sea_a_B"){
        E_sea_a_B = E_sea_a_B - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent1 == "I_sea_a_B"){
        I_sea_a_B = I_sea_a_B - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent1 == "R_sea_a_B"){
        R_sea_a_B = R_sea_a_B - 1
        R_sea_NB = R_sea_NB + 1
      }
      parent2 = sample(c(rep("S_a_B", S_a_B), rep("E_a_B", E_a_B),rep("I_a_B", I_a_B),rep("R_a_B", R_a_B),
                         rep("S_sea_a_B", S_sea_a_B), rep("E_sea_a_B", E_sea_a_B),rep("I_sea_a_B", I_sea_a_B),rep("R_sea_a_B", R_sea_a_B)),
                       size = 1)
      if (parent2 == "S_a_B"){
        S_a_B = S_a_B - 1
        S_a_NB = S_a_NB + 1
      } else if (parent2 == "E_a_B"){
        E_a_B= E_a_B- 1
        E_a_NB = E_a_NB + 1
      } else if (parent2 == "I_a_B"){
        I_a_B = I_a_B - 1
        I_a_NB = I_a_NB + 1
      } else if (parent2 == "R_a_B"){
        R_a_B = R_a_B - 1
        R_a_NB = R_a_NB + 1
      } else if (parent2 == "S_sea_a_B"){
        S_sea_a_B = S_sea_a_B - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent2 == "E_sea_a_B"){
        E_sea_a_B = E_sea_a_B - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent2 == "I_sea_a_B"){
        I_sea_a_B = I_sea_a_B - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent2 == "R_sea_a_B"){
        R_sea_a_B = R_sea_a_B - 1
        R_sea_NB = R_sea_NB + 1
      }
      ### In B
    }else if (transition == "S_b_N_to_E_b_N" & S_b_N > 0){
      S_b_N = S_b_N - 1
      E_b_N = E_b_N + 1
    }else if (transition == "E_b_N_to_S_b_N" & E_b_N > 0){
      E_b_N = E_b_N - 1
      S_b_N = S_b_N + 1
    }else if (transition == "E_b_N_to_I_b_N" & E_b_N > 0){
      E_b_N = E_b_N - 1
      I_b_N = I_b_N + 1
    }else if (transition == "I_b_N_to_R_b_N" & I_b_N > 0){
      I_b_N = I_b_N - 1
      R_b_N = R_b_N + 1
    }else if (transition == "I_b_N_to_R_b_N" & I_b_N > 0){
      I_b_N = I_b_N - 1
      R_b_N = R_b_N + 1
    }else if (transition == "I_b_N_to_D_b_N" & I_b_N > 0){
      I_b_N = I_b_N - 1
      D_b_N = D_b_N + 1
      # If a nestling dies, both parents become non-breeders.
      parent1 = sample(c(rep("S_b_B", S_b_B), rep("E_b_B", E_b_B),rep("I_b_B", I_b_B),rep("R_b_B", R_b_B),
                         rep("S_sea_b_B", S_sea_b_B), rep("E_sea_b_B", E_sea_b_B),rep("I_sea_b_B", I_sea_b_B),rep("R_sea_b_B", R_sea_b_B)),
                       size = 1)
      if (parent1 == "S_b_B"){
        S_b_B = S_b_B - 1
        S_b_NB = S_b_NB + 1
      } else if (parent1 == "E_b_B"){
        E_b_B = E_b_B - 1
        E_b_NB = E_b_NB + 1
      } else if (parent1 == "I_b_B"){
        I_b_B = I_b_B - 1
        I_b_NB = I_b_NB + 1
      } else if (parent1 == "R_b_B"){
        R_b_B = R_b_B - 1
        R_b_NB = R_b_NB + 1
      } else if (parent1 == "S_sea_b_B"){
        S_sea_b_B = S_sea_b_B - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent1 == "E_sea_b_B"){
        E_sea_b_B = E_sea_b_B - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent1 == "I_sea_b_B"){
        I_sea_b_B = I_sea_b_B - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent1 == "R_sea_b_B"){
        R_sea_b_B = R_sea_b_B - 1
        R_sea_NB = R_sea_NB + 1
      }
      parent2 = sample(c(rep("S_b_B", S_b_B), rep("E_b_B", E_b_B),rep("I_b_B", I_b_B),rep("R_b_B", R_b_B),
                         rep("S_sea_b_B", S_sea_b_B), rep("E_sea_b_B", E_sea_b_B),rep("I_sea_b_B", I_sea_b_B),rep("R_sea_b_B", R_sea_b_B)),
                       size = 1)
      if (parent2 == "S_b_B"){
        S_b_B = S_b_B - 1
        S_b_NB = S_b_NB + 1
      } else if (parent2 == "E_b_B"){
        E_b_B = E_b_B - 1
        E_b_NB = E_b_NB + 1
      } else if (parent2 == "I_b_B"){
        I_b_B = I_b_B - 1
        I_b_NB = I_b_NB + 1
      } else if (parent2 == "R_b_B"){
        R_b_B = R_b_B - 1
        R_b_NB = R_b_NB + 1
      } else if (parent2 == "S_sea_b_B"){
        S_sea_b_B = S_sea_b_B - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent2 == "E_sea_b_B"){
        E_sea_b_B = E_sea_b_B - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent2 == "I_sea_b_B"){
        I_sea_b_B = I_sea_b_B - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent2 == "R_sea_b_B"){
        R_sea_b_B = R_sea_b_B - 1
        R_sea_NB = R_sea_NB + 1
        }
      ### In C
      }else if (transition == "S_c_N_to_E_c_N" & S_c_N > 0){
        S_c_N = S_c_N - 1
        E_c_N = E_c_N + 1
      }else if (transition == "E_c_N_to_S_c_N" & E_c_N > 0){
        E_c_N = E_c_N - 1
        S_c_N = S_c_N + 1
      }else if (transition == "E_c_N_to_I_c_N" & E_c_N > 0){
        E_c_N = E_c_N - 1
        I_c_N = I_c_N + 1
      }else if (transition == "I_c_N_to_R_c_N" & I_c_N > 0){
        I_c_N = I_c_N - 1
        R_c_N = R_c_N + 1
      }else if (transition == "I_c_N_to_R_c_N" & I_c_N > 0){
        I_c_N = I_c_N - 1
        R_c_N = R_c_N + 1
      }else if (transition == "I_c_N_to_D_c_N" & I_c_N > 0){
        I_c_N = I_c_N - 1
        D_c_N = D_c_N + 1
        # If a nestling dies, both parents become non-breeders.
        parent1 = sample(c(rep("S_c_B", S_c_B), rep("E_c_B", E_c_B),rep("I_c_B", I_c_B),rep("R_c_B", R_c_B),
                           rep("S_sea_c_B", S_sea_c_B), rep("E_sea_c_B", E_sea_c_B),rep("I_sea_c_B", I_sea_c_B),rep("R_sea_c_B", R_sea_c_B)),
                         size = 1)
        if (parent1 == "S_c_B"){
          S_c_B = S_c_B - 1
          S_c_NB = S_c_NB + 1
        } else if (parent1 == "E_c_B"){
          E_c_B = E_c_B - 1
          E_c_NB = E_c_NB + 1
        } else if (parent1 == "I_c_B"){
          I_c_B = I_c_B - 1
          I_c_NB = I_c_NB + 1
        } else if (parent1 == "R_c_B"){
          R_c_B = R_c_B - 1
          R_c_NB = R_c_NB + 1
        } else if (parent1 == "S_sea_c_B"){
          S_sea_c_B = S_sea_c_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (parent1 == "E_sea_c_B"){
          E_sea_c_B = E_sea_c_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (parent1 == "I_sea_c_B"){
          I_sea_c_B = I_sea_c_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (parent1 == "R_sea_c_B"){
          R_sea_c_B = R_sea_c_B - 1
          R_sea_NB = R_sea_NB + 1
        }
        parent2 = sample(c(rep("S_c_B", S_c_B), rep("E_c_B", E_c_B),rep("I_c_B", I_c_B),rep("R_c_B", R_c_B),
                           rep("S_sea_c_B", S_sea_c_B), rep("E_sea_c_B", E_sea_c_B),rep("I_sea_c_B", I_sea_c_B),rep("R_sea_c_B", R_sea_c_B)),
                         size = 1)
        if (parent2 == "S_c_B"){
          S_c_B = S_c_B - 1
          S_c_NB = S_c_NB + 1
        } else if (parent2 == "E_c_B"){
          E_c_B = E_c_B - 1
          E_c_NB = E_c_NB + 1
        } else if (parent2 == "I_c_B"){
          I_c_B = I_c_B - 1
          I_c_NB = I_c_NB + 1
        } else if (parent2 == "R_c_B"){
          R_c_B = R_c_B - 1
          R_c_NB = R_c_NB + 1
        } else if (parent2 == "S_sea_c_B"){
          S_sea_c_B = S_sea_c_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (parent2 == "E_sea_c_B"){
          E_sea_c_B = E_sea_c_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (parent2 == "I_sea_c_B"){
          I_sea_c_B = I_sea_c_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (parent2 == "R_sea_c_B"){
          R_sea_c_B = R_sea_c_B - 1
          R_sea_NB = R_sea_NB + 1
        }
      ## Breeders
      ### In A
      #### At colony
      } else if  (transition == "S_a_B_to_E_a_B" & S_a_B > 0) {
        S_a_B = S_a_B - 1
        E_a_B= E_a_B+ 1
      } else if (transition == "E_a_B_to_S_a_B" & E_a_B> 0) {
        E_a_B= E_a_B- 1
        S_a_B = S_a_B + 1
      } else if (transition == "E_a_B_to_I_a_B" & E_a_B> 0) {
        E_a_B= E_a_B- 1
        I_a_B = I_a_B + 1
      } else if (transition == "I_a_B_to_R_a_B" & I_a_B > 0) {
        I_a_B = I_a_B - 1
        R_a_B = R_a_B + 1
      } else if (transition == "I_a_B_to_D_a_B" & I_a_B > 0) {
        # If an adult dies, the partner becomes a non-breeder and the nestling dies
        I_a_B = I_a_B - 1
        D_a_B = D_a_B + 1
        # The nestling dies
        if (S_a_N + E_a_N + I_a_N + R_a_N > 0){
          nestling = sample(c(rep("S_a_N", S_a_N), rep("E_a_N", E_a_N),rep("I_a_N", I_a_N),rep("R_a_N", R_a_N)),
                            size = 1)
          if (nestling == "S_a_N"){
            S_a_N = S_a_N - 1
            D_a_N = D_a_N + 1
          } else if (nestling == "E_a_N"){
            E_a_N = E_a_N - 1
            D_a_N = D_a_N + 1
          } else if (nestling == "I_a_N"){
            I_a_N = I_a_N - 1
            D_a_N = D_a_N + 1
          } else if (nestling == "R_a_N"){
            R_a_N = R_a_N - 1
            D_a_N = D_a_N + 1
          }
        }
        # The partner becomes a non-breeder
        partner = sample(c(rep("S_a_B", S_a_B), rep("E_a_B", E_a_B),rep("I_a_B", I_a_B),rep("R_a_B", R_a_B),
                           rep("S_sea_a_B", S_sea_a_B), rep("E_sea_a_B", E_sea_a_B),rep("I_sea_a_B", I_sea_a_B),rep("R_sea_a_B", R_sea_a_B)),
                         size = 1)
        if (partner == "S_a_B"){
          S_a_B = S_a_B - 1
          S_a_NB = S_a_NB + 1
        } else if (partner == "E_a_B"){
          E_a_B= E_a_B- 1
          E_a_NB = E_a_NB + 1
        } else if (partner == "I_a_B"){
          I_a_B = I_a_B - 1
          I_a_NB = I_a_NB + 1
        } else if (partner == "R_a_B"){
          R_a_B = R_a_B - 1
          R_a_NB = R_a_NB + 1
        } else if (partner == "S_sea_a_B"){
          S_sea_a_B = S_sea_a_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_a_B"){
          E_sea_a_B = E_sea_a_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_a_B"){
          I_sea_a_B = I_sea_a_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_a_B"){
          R_sea_a_B = R_sea_a_B - 1
          R_sea_NB = R_sea_NB + 1
        }
      #### At sea
      } else if (transition == "S_sea_a_B_to_E_sea_a_B" & S_sea_a_B > 0) {
        S_sea_a_B = S_sea_a_B - 1
        E_sea_a_B = E_sea_a_B + 1
      } else if (transition == "E_sea_a_B_to_S_sea_a_B" & E_sea_a_B > 0) {
        E_sea_a_B = E_sea_a_B - 1
        S_sea_a_B = S_sea_a_B + 1
      } else if (transition == "E_sea_a_B_to_I_sea_a_B" & E_sea_a_B > 0) {
        E_sea_a_B = E_sea_a_B - 1
        I_sea_a_B = I_sea_a_B + 1
      } else if (transition == "I_sea_a_B_to_R_sea_a_B" & I_sea_a_B > 0) {
        I_sea_a_B = I_sea_a_B - 1
        R_sea_a_B = R_sea_a_B + 1
      } else if (transition == "I_sea_a_B_to_D_sea_a_B" & I_sea_a_B > 0) {
        I_sea_a_B = I_sea_a_B - 1
        D_sea_a_B = D_sea_a_B + 1
        # The nestling dies
        if (S_a_N + E_a_N + I_a_N + R_a_N > 0){
          nestling = sample(c(rep("S_a_N", S_a_N), rep("E_a_N", E_a_N),rep("I_a_N", I_a_N),rep("R_a_N", R_a_N)),
                            size = 1)
          if (nestling == "S_a_N"){
            S_a_N = S_a_N - 1
            D_a_N = D_a_N + 1
          } else if (nestling == "E_a_N"){
            E_a_N = E_a_N - 1
            D_a_N = D_a_N + 1
          } else if (nestling == "I_a_N"){
            I_a_N = I_a_N - 1
            D_a_N = D_a_N + 1
          } else if (nestling == "R_a_N"){
            R_a_N = R_a_N - 1
            D_a_N = D_a_N + 1
          }
        }
        # The partner becomes a non-breeder
        partner = sample(c(rep("S_a_B", S_a_B), rep("E_a_B", E_a_B),rep("I_a_B", I_a_B),rep("R_a_B", R_a_B),
                           rep("S_sea_a_B", S_sea_a_B), rep("E_sea_a_B", E_sea_a_B),rep("I_sea_a_B", I_sea_a_B),rep("R_sea_a_B", R_sea_a_B)),
                         size = 1)
        if (partner == "S_a_B"){
          S_a_B = S_a_B - 1
          S_a_NB = S_a_NB + 1
        } else if (partner == "E_a_B"){
          E_a_B= E_a_B- 1
          E_a_NB = E_a_NB + 1
        } else if (partner == "I_a_B"){
          I_a_B = I_a_B - 1
          I_a_NB = I_a_NB + 1
        } else if (partner == "R_a_B"){
          R_a_B = R_a_B - 1
          R_a_NB = R_a_NB + 1
        } else if (partner == "S_sea_a_B"){
          S_sea_a_B = S_sea_a_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_a_B"){
          E_sea_a_B = E_sea_a_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_a_B"){
          I_sea_a_B = I_sea_a_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_a_B"){
          R_sea_a_B = R_sea_a_B - 1
          R_sea_NB = R_sea_NB + 1
        }
      ### In B
      #### At colony
      } else if (transition == "S_b_B_to_E_b_B" & S_b_B > 0) {
        S_b_B = S_b_B - 1
        E_b_B = E_b_B + 1
      } else if (transition == "E_b_B_to_S_b_B" & E_b_B > 0) {
        E_b_B = E_b_B - 1
        S_b_B = S_b_B + 1
      } else if (transition == "E_b_B_to_I_b_B" & E_b_B > 0) {
        E_b_B = E_b_B - 1
        I_b_B = I_b_B + 1
      } else if (transition == "I_b_B_to_R_b_B" & I_b_B > 0) {
        I_b_B = I_b_B - 1
        R_b_B = R_b_B + 1
      } else if (transition == "I_b_B_to_D_b_B" & I_b_B > 0) {
        # If an adult dies, the partner becomes a non-breeder and the nestling dies
        I_b_B = I_b_B - 1
        D_b_B = D_b_B + 1
        # The nestling dies
        if (S_b_N + E_b_N + I_b_N + R_b_N > 0){
          nestling = sample(c(rep("S_b_N", S_b_N), rep("E_b_N", E_b_N),rep("I_b_N", I_b_N),rep("R_b_N", R_b_N)),
                            size = 1)
          if (nestling == "S_b_N"){
            S_b_N = S_b_N - 1
            D_b_N = D_b_N + 1
          } else if (nestling == "E_b_N"){
            E_b_N = E_b_N - 1
            D_b_N = D_b_N + 1
          } else if (nestling == "I_b_N"){
            I_b_N = I_b_N - 1
            D_b_N = D_b_N + 1
          } else if (nestling == "R_b_N"){
            R_b_N = R_b_N - 1
            D_b_N = D_b_N + 1
          }
        }
        # The partner becomes a non-breeder
        partner = sample(c(rep("S_b_B", S_b_B), rep("E_b_B", E_b_B),rep("I_b_B", I_b_B),rep("R_b_B", R_b_B),
                           rep("S_sea_b_B", S_sea_b_B), rep("E_sea_b_B", E_sea_b_B),rep("I_sea_b_B", I_sea_b_B),rep("R_sea_b_B", R_sea_b_B)),
                         size = 1)
       
        if (partner == "S_b_B"){
          S_b_B = S_b_B - 1
          S_b_NB = S_b_NB + 1
        } else if (partner == "E_b_B"){
          E_b_B = E_b_B - 1
          E_b_NB = E_b_NB + 1
        } else if (partner == "I_b_B"){
          I_b_B = I_b_B - 1
          I_b_NB = I_b_NB + 1
        } else if (partner == "R_b_B"){
          R_b_B = R_b_B - 1
          R_b_NB = R_b_NB + 1
        } else if (partner == "S_sea_b_B"){
          S_sea_b_B = S_sea_b_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_b_B"){
          E_sea_b_B = E_sea_b_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_b_B"){
          I_sea_b_B = I_sea_b_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_b_B"){
          R_sea_b_B = R_sea_b_B - 1
          R_sea_NB = R_sea_NB + 1
        }
      #### At sea
      } else if (transition == "S_sea_b_B_to_E_sea_b_B" & S_sea_b_B > 0) {
        S_sea_b_B = S_sea_b_B - 1
        E_sea_b_B = E_sea_b_B + 1
      } else if (transition == "E_sea_b_B_to_S_sea_b_B" & E_sea_b_B > 0) {
        E_sea_b_B = E_sea_b_B - 1
        S_sea_b_B = S_sea_b_B + 1
      } else if (transition == "E_sea_b_B_to_I_sea_b_B" & E_sea_b_B > 0) {
        E_sea_b_B = E_sea_b_B - 1
        I_sea_b_B = I_sea_b_B + 1
      } else if (transition == "I_sea_b_B_to_R_sea_b_B" & I_sea_b_B > 0) {
        I_sea_b_B = I_sea_b_B - 1
        R_sea_b_B = R_sea_b_B + 1
      } else if (transition == "I_sea_b_B_to_D_sea_b_B" & I_sea_b_B > 0) {
        I_sea_b_B = I_sea_b_B - 1
        D_sea_b_B = D_sea_b_B + 1
        # The nestling dies
        if (S_b_N + E_b_N + I_b_N + R_b_N > 0){
          nestling = sample(c(rep("S_b_N", S_b_N), rep("E_b_N", E_b_N),rep("I_b_N", I_b_N),rep("R_b_N", R_b_N)),
                            size = 1)
          if (nestling == "S_b_N"){
            S_b_N = S_b_N - 1
            D_b_N = D_b_N + 1
          } else if (nestling == "E_b_N"){
            E_b_N = E_b_N - 1
            D_b_N = D_b_N + 1
          } else if (nestling == "I_b_N"){
            I_b_N = I_b_N - 1
            D_b_N = D_b_N + 1
          } else if (nestling == "R_b_N"){
            R_b_N = R_b_N - 1
            D_b_N = D_b_N + 1
          }
        }
        # The partner becomes a non-breeder
        partner = sample(c(rep("S_b_B", S_b_B), rep("E_b_B", E_b_B),rep("I_b_B", I_b_B),rep("R_b_B", R_b_B),
                           rep("S_sea_b_B", S_sea_b_B), rep("E_sea_b_B", E_sea_b_B),rep("I_sea_b_B", I_sea_b_B),rep("R_sea_b_B", R_sea_b_B)),
                         size = 1)
        
        if (partner == "S_b_B"){
          S_b_B = S_b_B - 1
          S_b_NB = S_b_NB + 1
        } else if (partner == "E_b_B"){
          E_b_B = E_b_B - 1
          E_b_NB = E_b_NB + 1
        } else if (partner == "I_b_B"){
          I_b_B = I_b_B - 1
          I_b_NB = I_b_NB + 1
        } else if (partner == "R_b_B"){
          R_b_B = R_b_B - 1
          R_b_NB = R_b_NB + 1
        } else if (partner == "S_sea_b_B"){
          S_sea_b_B = S_sea_b_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_b_B"){
          E_sea_b_B = E_sea_b_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_b_B"){
          I_sea_b_B = I_sea_b_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_b_B"){
          R_sea_b_B = R_sea_b_B - 1
          R_sea_NB = R_sea_NB + 1
        }
      ### In C
      #### At colony
      } else if (transition == "S_c_B_to_E_c_B" & S_c_B > 0) {
        S_c_B = S_c_B - 1
        E_c_B = E_c_B + 1
      } else if (transition == "E_c_B_to_S_c_B" & E_c_B > 0) {
        E_c_B = E_c_B - 1
        S_c_B = S_c_B + 1
      } else if (transition == "E_c_B_to_I_c_B" & E_c_B > 0) {
        E_c_B = E_c_B - 1
        I_c_B = I_c_B + 1
      } else if (transition == "I_c_B_to_R_c_B" & I_c_B > 0) {
        I_c_B = I_c_B - 1
        R_c_B = R_c_B + 1
      } else if (transition == "I_c_B_to_D_c_B" & I_c_B > 0) {
        # If an adult dies, the partner becomes a non-breeder and the nestling dies
        I_c_B = I_c_B - 1
        D_c_B = D_c_B + 1
        # The nestling dies
        if (S_c_N + E_c_N + I_c_N + R_c_N > 0){
          nestling = sample(c(rep("S_c_N", S_c_N), rep("E_c_N", E_c_N),rep("I_c_N", I_c_N),rep("R_c_N", R_c_N)),
                            size = 1)
          if (nestling == "S_c_N"){
            S_c_N = S_c_N - 1
            D_c_N = D_c_N + 1
          } else if (nestling == "E_c_N"){
            E_c_N = E_c_N - 1
            D_c_N = D_c_N + 1
          } else if (nestling == "I_c_N"){
            I_c_N = I_c_N - 1
            D_c_N = D_c_N + 1
          } else if (nestling == "R_c_N"){
            R_c_N = R_c_N - 1
            D_c_N = D_c_N + 1
          }
        }
        # The partner becomes a non-breeder
        partner = sample(c(rep("S_c_B", S_c_B), rep("E_c_B", E_c_B),rep("I_c_B", I_c_B),rep("R_c_B", R_c_B),
                           rep("S_sea_c_B", S_sea_c_B), rep("E_sea_c_B", E_sea_c_B),rep("I_sea_c_B", I_sea_c_B),rep("R_sea_c_B", R_sea_c_B)),
                         size = 1)
        
        if (partner == "S_c_B"){
          S_c_B = S_c_B - 1
          S_c_NB = S_c_NB + 1
        } else if (partner == "E_c_B"){
          E_c_B = E_c_B - 1
          E_c_NB = E_c_NB + 1
        } else if (partner == "I_c_B"){
          I_c_B = I_c_B - 1
          I_c_NB = I_c_NB + 1
        } else if (partner == "R_c_B"){
          R_c_B = R_c_B - 1
          R_c_NB = R_c_NB + 1
        } else if (partner == "S_sea_c_B"){
          S_sea_c_B = S_sea_c_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_c_B"){
          E_sea_c_B = E_sea_c_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_c_B"){
          I_sea_c_B = I_sea_c_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_c_B"){
          R_sea_c_B = R_sea_c_B - 1
          R_sea_NB = R_sea_NB + 1
        }
        #### At sea
      } else if (transition == "S_sea_c_B_to_E_sea_c_B" & S_sea_c_B > 0) {
        S_sea_c_B = S_sea_c_B - 1
        E_sea_c_B = E_sea_c_B + 1
      } else if (transition == "E_sea_c_B_to_S_sea_c_B" & E_sea_c_B > 0) {
        E_sea_c_B = E_sea_c_B - 1
        S_sea_c_B = S_sea_c_B + 1
      } else if (transition == "E_sea_c_B_to_I_sea_c_B" & E_sea_c_B > 0) {
        E_sea_c_B = E_sea_c_B - 1
        I_sea_c_B = I_sea_c_B + 1
      } else if (transition == "I_sea_c_B_to_R_sea_c_B" & I_sea_c_B > 0) {
        I_sea_c_B = I_sea_c_B - 1
        R_sea_c_B = R_sea_c_B + 1
      } else if (transition == "I_sea_c_B_to_D_sea_c_B" & I_sea_c_B > 0) {
        I_sea_c_B = I_sea_c_B - 1
        D_sea_c_B = D_sea_c_B + 1
        # The nestling dies
        if (S_c_N + E_c_N + I_c_N + R_c_N > 0){
          nestling = sample(c(rep("S_c_N", S_c_N), rep("E_c_N", E_c_N),rep("I_c_N", I_c_N),rep("R_c_N", R_c_N)),
                            size = 1)
          if (nestling == "S_c_N"){
            S_c_N = S_c_N - 1
            D_c_N = D_c_N + 1
          } else if (nestling == "E_c_N"){
            E_c_N = E_c_N - 1
            D_c_N = D_c_N + 1
          } else if (nestling == "I_c_N"){
            I_c_N = I_c_N - 1
            D_c_N = D_c_N + 1
          } else if (nestling == "R_c_N"){
            R_c_N = R_c_N - 1
            D_c_N = D_c_N + 1
          }
        }
        # The partner becomes a non-breeder
        partner = sample(c(rep("S_c_B", S_c_B), rep("E_c_B", E_c_B),rep("I_c_B", I_c_B),rep("R_c_B", R_c_B),
                           rep("S_sea_c_B", S_sea_c_B), rep("E_sea_c_B", E_sea_c_B),rep("I_sea_c_B", I_sea_c_B),rep("R_sea_c_B", R_sea_c_B)),
                         size = 1)
        
        if (partner == "S_c_B"){
          S_c_B = S_c_B - 1
          S_c_NB = S_c_NB + 1
        } else if (partner == "E_c_B"){
          E_c_B = E_c_B - 1
          E_c_NB = E_c_NB + 1
        } else if (partner == "I_c_B"){
          I_c_B = I_c_B - 1
          I_c_NB = I_c_NB + 1
        } else if (partner == "R_c_B"){
          R_c_B = R_c_B - 1
          R_c_NB = R_c_NB + 1
        } else if (partner == "S_sea_c_B"){
          S_sea_c_B = S_sea_c_B - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_c_B"){
          E_sea_c_B = E_sea_c_B - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_c_B"){
          I_sea_c_B = I_sea_c_B - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_c_B"){
          R_sea_c_B = R_sea_c_B - 1
          R_sea_NB = R_sea_NB + 1
        }
      
      ## Non-Breeders
      ### In colony A
      } else if (transition == "S_a_NB_to_E_a_NB" & S_a_NB > 0) {
        S_a_NB = S_a_NB - 1
        E_a_NB = E_a_NB + 1
      } else if (transition == "E_a_NB_to_S_a_NB" & E_a_NB > 0) {
        E_a_NB = E_a_NB - 1
        S_a_NB = S_a_NB + 1
      } else if (transition == "E_a_NB_to_I_a_NB" & E_a_NB > 0) {
        E_a_NB = E_a_NB - 1
        I_a_NB = I_a_NB + 1
      } else if (transition == "I_a_NB_to_R_a_NB" & I_a_NB > 0) {
        I_a_NB = I_a_NB - 1
        R_a_NB = R_a_NB + 1
      } else if (transition == "I_a_NB_to_D_a_NB" & I_a_NB > 0) {
        I_a_NB = I_a_NB - 1
        D_a_NB = D_a_NB + 1
      ### In colony B
      } else if (transition == "S_b_NB_to_E_b_NB" & S_b_NB > 0) {
        S_b_NB = S_b_NB - 1
        E_b_NB = E_b_NB + 1
      } else if (transition == "E_b_NB_to_S_b_NB" & E_b_NB > 0) {
        E_b_NB = E_b_NB - 1
        S_b_NB = S_b_NB + 1
      } else if (transition == "E_b_NB_to_I_b_NB" & E_b_NB > 0) {
        E_b_NB = E_b_NB - 1
        I_b_NB = I_b_NB + 1
      } else if (transition == "I_b_NB_to_R_b_NB" & I_b_NB > 0) {
        I_b_NB = I_b_NB - 1
        R_b_NB = R_b_NB + 1
      } else if (transition == "I_b_NB_to_D_b_NB" & I_b_NB > 0) {
        I_b_NB = I_b_NB - 1
        D_b_NB = D_b_NB + 1
      ### In colony C
      } else if (transition == "S_c_NB_to_E_c_NB" & S_c_NB > 0) {
        S_c_NB = S_c_NB - 1
        E_c_NB = E_c_NB + 1
      } else if (transition == "E_c_NB_to_S_c_NB" & E_c_NB > 0) {
        E_c_NB = E_c_NB - 1
        S_c_NB = S_c_NB + 1
      } else if (transition == "E_c_NB_to_I_c_NB" & E_c_NB > 0) {
        E_c_NB = E_c_NB - 1
        I_c_NB = I_c_NB + 1
      } else if (transition == "I_c_NB_to_R_c_NB" & I_c_NB > 0) {
        I_c_NB = I_c_NB - 1
        R_c_NB = R_c_NB + 1
      } else if (transition == "I_c_NB_to_D_c_NB" & I_c_NB > 0) {
        I_c_NB = I_c_NB - 1
        D_c_NB = D_c_NB + 1
        ### At sea
      } else if (transition == "S_sea_NB_to_E_sea_NB" & S_sea_NB > 0) {
        S_sea_NB = S_sea_NB - 1
        E_sea_NB = E_sea_NB + 1
      } else if (transition == "E_sea_NB_to_S_sea_NB" & E_sea_NB > 0) {
        E_sea_NB = E_sea_NB - 1
        S_sea_NB = S_sea_NB + 1
      } else if (transition == "E_sea_NB_to_I_sea_NB" & E_sea_NB > 0) {
        E_sea_NB = E_sea_NB - 1
        I_sea_NB = I_sea_NB + 1
      } else if (transition == "I_sea_NB_to_R_sea_NB" & I_sea_NB > 0) {
        I_sea_NB = I_sea_NB - 1
        R_sea_NB = R_sea_NB + 1
      } else if (transition == "I_sea_NB_to_D_sea_NB" & I_sea_NB > 0) {
        I_sea_NB = I_sea_NB - 1
        D_sea_NB = D_sea_NB + 1

        # Mobility
        ## Breeders
        ### From colony A to sea A
      } else if (transition == "S_a_B_to_S_sea_a_B"  & S_a_B > 0) {
        S_a_B = S_a_B - 1
        S_sea_a_B = S_sea_a_B + 1
      } else if (transition == "E_a_B_to_E_sea_a_B"  & E_a_B> 0) {
        E_a_B= E_a_B- 1
        E_sea_a_B = E_sea_a_B + 1
      } else if (transition == "I_a_B_to_I_sea_a_B"  & I_a_B > 0) {
        I_a_B = I_a_B - 1
        I_sea_a_B = I_sea_a_B + 1
      } else if (transition == "R_a_B_to_R_sea_a_B"  & R_a_B > 0) {
        R_a_B = R_a_B - 1
        R_sea_a_B = R_sea_a_B + 1
      ### From sea A to colony A
      }else if (transition == "S_sea_a_B_to_S_a_B"  & S_sea_a_B > 0) {
        S_a_B = S_a_B + 1
        S_sea_a_B = S_sea_a_B - 1
      } else if (transition == "E_sea_a_B_to_E_a_B"  & E_sea_a_B > 0) {
        E_a_B= E_a_B+ 1
        E_sea_a_B = E_sea_a_B - 1
      } else if (transition == "I_sea_a_B_to_I_a_B"  & I_sea_a_B > 0) {
        I_a_B = I_a_B + 1
        I_sea_a_B = I_sea_a_B - 1
      } else if (transition == "R_sea_a_B_to_R_a_B" & R_sea_a_B > 0) {
        R_a_B = R_a_B + 1
        R_sea_a_B = R_sea_a_B - 1
       ### From colony B to sea B
      } else if (transition == "S_b_B_to_S_sea_b_B"  & S_b_B > 0) {
        S_b_B = S_b_B - 1
        S_sea_b_B = S_sea_b_B + 1
      } else if (transition == "E_b_B_to_E_sea_b_B"  & E_b_B > 0) {
        E_b_B = E_b_B - 1
        E_sea_b_B = E_sea_b_B + 1
      } else if (transition == "I_b_B_to_I_sea_b_B"  & I_b_B > 0) {
        I_b_B = I_b_B - 1
        I_sea_b_B = I_sea_b_B + 1
      } else if (transition == "R_b_B_to_R_sea_b_B"  & R_b_B > 0) {
        R_b_B = R_b_B - 1
        R_sea_b_B = R_sea_b_B + 1
      ### From sea B to colony B
      } else if (transition == "S_sea_b_B_to_S_b_B" & S_sea_b_B > 0) {
        S_b_B = S_b_B + 1
        S_sea_b_B = S_sea_b_B - 1
      } else if (transition == "E_sea_b_B_to_E_b_B" & E_sea_b_B > 0) {
        E_b_B = E_b_B + 1
        E_sea_b_B = E_sea_b_B - 1
      } else if (transition == "I_sea_b_B_to_I_b_B" & I_sea_b_B > 0) {
        I_b_B = I_b_B + 1
        I_sea_b_B = I_sea_b_B - 1
      } else if (transition == "R_sea_b_B_to_R_b_B" & R_sea_b_B > 0) {
        R_b_B = R_b_B + 1
        R_sea_b_B = R_sea_b_B - 1 
      ### From colony C to sea C
      } else if (transition == "S_c_B_to_S_sea_c_B"  & S_c_B > 0) {
        S_c_B = S_c_B - 1
        S_sea_c_B = S_sea_c_B + 1
      } else if (transition == "E_c_B_to_E_sea_c_B"  & E_c_B > 0) {
        E_c_B = E_c_B - 1
        E_sea_c_B = E_sea_c_B + 1
      } else if (transition == "I_c_B_to_I_sea_c_B"  & I_c_B > 0) {
        I_c_B = I_c_B - 1
        I_sea_c_B = I_sea_c_B + 1
      } else if (transition == "R_c_B_to_R_sea_c_B"  & R_c_B > 0) {
        R_c_B = R_c_B - 1
        R_sea_c_B = R_sea_c_B + 1
        ### From sea C to colony C
      } else if (transition == "S_sea_c_B_to_S_c_B" & S_sea_c_B > 0) {
        S_c_B = S_c_B + 1
        S_sea_c_B = S_sea_c_B - 1
      } else if (transition == "E_sea_c_B_to_E_c_B" & E_sea_c_B > 0) {
        E_c_B = E_c_B + 1
        E_sea_c_B = E_sea_c_B - 1
      } else if (transition == "I_sea_c_B_to_I_c_B" & I_sea_c_B > 0) {
        I_c_B = I_c_B + 1
        I_sea_c_B = I_sea_c_B - 1
      } else if (transition == "R_sea_c_B_to_R_c_B" & R_sea_c_B > 0) {
        R_c_B = R_c_B + 1
        R_sea_c_B = R_sea_c_B - 1 
        
      ## Non-Breeders
      ### From colony A to sea 
      } else if (transition == "S_a_NB_to_S_sea_NB" & S_a_NB > 0) {
        S_a_NB = S_a_NB - 1
        S_sea_NB = S_sea_NB + 1
      } else if (transition == "E_a_NB_to_E_sea_NB" & E_a_NB > 0) {
        E_a_NB = E_a_NB - 1
        E_sea_NB = E_sea_NB + 1
      } else if (transition == "I_a_NB_to_I_sea_NB" & I_a_NB > 0) {
        I_a_NB = I_a_NB - 1
        I_sea_NB = I_sea_NB + 1
      } else if (transition == "R_a_NB_to_R_sea_NB" & R_a_NB > 0) {
        R_a_NB = R_a_NB - 1
        R_sea_NB = R_sea_NB + 1
      ### From sea to colony A
      }else if (transition == "S_sea_NB_to_S_a_NB" & S_sea_NB > 0) {
        S_a_NB = S_a_NB + 1
        S_sea_NB = S_sea_NB - 1
      } else if (transition == "E_sea_NB_to_E_a_NB" & E_sea_NB > 0) {
        E_a_NB = E_a_NB + 1
        E_sea_NB = E_sea_NB - 1
      } else if (transition == "I_sea_NB_to_I_a_NB" & I_sea_NB > 0) {
        I_a_NB = I_a_NB + 1
        I_sea_NB = I_sea_NB - 1
      } else if (transition == "R_sea_NB_to_R_a_NB" & R_sea_NB > 0) {
        R_a_NB = R_a_NB + 1
        R_sea_NB = R_sea_NB - 1
      ### From colony B to sea 
      } else if (transition == "S_b_NB_to_S_sea_NB" & S_b_NB > 0) {
        S_b_NB = S_b_NB - 1
        S_sea_NB = S_sea_NB + 1
      } else if (transition == "E_b_NB_to_E_sea_NB" & E_b_NB > 0) {
        E_b_NB = E_b_NB - 1
        E_sea_NB = E_sea_NB + 1
      } else if (transition == "I_b_NB_to_I_sea_NB" & I_b_NB > 0) {
        I_b_NB = I_b_NB - 1
        I_sea_NB = I_sea_NB + 1
      } else if (transition == "R_b_NB_to_R_sea_NB" & R_b_NB > 0) {
        R_b_NB = R_b_NB - 1
        R_sea_NB = R_sea_NB + 1
      ### From sea to colony B
      } else if (transition == "S_sea_NB_to_S_b_NB" & S_sea_NB > 0) {
        S_b_NB = S_b_NB + 1
        S_sea_NB = S_sea_NB - 1
      } else if (transition == "E_sea_NB_to_E_b_NB" & E_sea_NB > 0) {
        E_b_NB = E_b_NB + 1
        E_sea_NB = E_sea_NB - 1
      } else if (transition == "I_sea_NB_to_I_b_NB" & I_sea_NB > 0) {
        I_b_NB = I_b_NB + 1
        I_sea_NB = I_sea_NB - 1
      } else if (transition == "R_sea_NB_to_R_b_NB" & R_sea_NB > 0) {
        R_b_NB = R_b_NB + 1
        R_sea_NB = R_sea_NB - 1
        ### From colony C to sea 
      } else if (transition == "S_c_NB_to_S_sea_NB" & S_c_NB > 0) {
        S_c_NB = S_c_NB - 1
        S_sea_NB = S_sea_NB + 1
      } else if (transition == "E_c_NB_to_E_sea_NB" & E_c_NB > 0) {
        E_c_NB = E_c_NB - 1
        E_sea_NB = E_sea_NB + 1
      } else if (transition == "I_c_NB_to_I_sea_NB" & I_c_NB > 0) {
        I_c_NB = I_c_NB - 1
        I_sea_NB = I_sea_NB + 1
      } else if (transition == "R_c_NB_to_R_sea_NB" & R_c_NB > 0) {
        R_c_NB = R_c_NB - 1
        R_sea_NB = R_sea_NB + 1
        ### From sea to colony C
      } else if (transition == "S_sea_NB_to_S_c_NB" & S_sea_NB > 0) {
        S_c_NB = S_c_NB + 1
        S_sea_NB = S_sea_NB - 1
      } else if (transition == "E_sea_NB_to_E_c_NB" & E_sea_NB > 0) {
        E_c_NB = E_c_NB + 1
        E_sea_NB = E_sea_NB - 1
      } else if (transition == "I_sea_NB_to_I_c_NB" & I_sea_NB > 0) {
        I_c_NB = I_c_NB + 1
        I_sea_NB = I_sea_NB - 1
      } else if (transition == "R_sea_NB_to_R_c_NB" & R_sea_NB > 0) {
        R_c_NB = R_c_NB + 1
        R_sea_NB = R_sea_NB - 1
      # } else if (transition == "S_a_B_to_S_a_NB"  & S_a_B > 0) {
      #   S_a_B = S_a_B - 2
      #   S_a_NB = S_a_NB + 2
      # } else if (transition == "E_a_B_to_E_a_NB" & E_a_B> 0) {
      #   E_a_B= E_a_B- 2
      #   E_a_NB = E_a_NB + 2
      # } else if (transition == "I_a_B_to_I_a_NB" & I_a_B > 0) {
      #   I_a_B = I_a_B - 2
      #   I_a_NB = I_a_NB + 2
      # } else if (transition == "E_a_B_to_E_a_NB" & R_a_B > 0) {
      #   R_a_B = R_a_B - 2
      #   R_a_NB = R_a_NB + 2
        
        } # transition
      } # for : transitions_bank
    } # if : transitions_bank

    
    new_state = matrix(data = c(S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                S_a_B, E_a_B, I_a_B, R_a_B,  D_a_B, 
                                S_sea_a_B, E_sea_a_B, I_sea_a_B, R_sea_a_B, D_sea_a_B,
                                S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                
                                S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                S_b_B, E_b_B, I_b_B, R_b_B, D_b_B,
                                S_sea_b_B, E_sea_b_B, I_sea_b_B, R_sea_b_B, D_sea_b_B,
                                S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                
                                S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                S_c_B, E_c_B, I_c_B, R_c_B, D_c_B,
                                S_sea_c_B, E_sea_c_B, I_sea_c_B, R_sea_c_B, D_sea_c_B,
                                S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                                
                                S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB
    ),
    nrow = 13, ncol = 5, 
    byrow = T)
    
    
    
    states = abind(states, new_state)

  } # while
  
  all_states = data.frame(
    time = times,
    
    S_a_N = states[1, 1, ],
    E_a_N = states[1, 2, ],
    I_a_N = states[1, 3, ],
    R_a_N = states[1, 4, ],
    D_a_N = states[1, 5, ],
    
    S_a_B = states[2, 1, ],
    E_a_B= states[2, 2, ],
    I_a_B = states[2, 3, ],
    R_a_B = states[2, 4, ],
    D_a_B = states[2, 5, ],
    
    S_sea_a_B = states[3, 1, ],
    E_sea_a_B = states[3, 2, ],
    I_sea_a_B = states[3, 3, ],
    R_sea_a_B = states[3, 4, ],
    D_sea_a_B = states[3, 5, ],
    
    S_a_NB = states[4, 1, ],
    E_a_NB = states[4, 2, ],
    I_a_NB = states[4, 3, ],
    R_a_NB = states[4, 4, ],
    D_a_NB = states[4, 5, ],
    
    S_b_N = states[5, 1, ],
    E_b_N = states[5, 2, ],
    I_b_N = states[5, 3, ],
    R_b_N = states[5, 4, ],
    D_b_N = states[5, 5, ],
    
    S_b_B = states[6, 1, ],
    E_b_B = states[6, 2, ],
    I_b_B = states[6, 3, ],
    R_b_B = states[6, 4, ],
    D_b_B = states[6, 5, ],
    
    S_sea_b_B = states[7, 1, ],
    E_sea_b_B = states[7, 2, ],
    I_sea_b_B = states[7, 3, ],
    R_sea_b_B = states[7, 4, ],
    D_sea_b_B = states[7, 5, ],
  
    S_b_NB = states[8, 1, ],
    E_b_NB = states[8, 2, ],
    I_b_NB = states[8, 3, ],
    R_b_NB = states[8, 4, ],
    D_b_NB = states[8, 5, ],
    
    S_c_N = states[9, 1, ],
    E_c_N = states[9, 2, ],
    I_c_N = states[9, 3, ],
    R_c_N = states[9, 4, ],
    D_c_N = states[9, 5, ],
    
    S_c_B = states[10, 1, ],
    E_c_B = states[10, 2, ],
    I_c_B = states[10, 3, ],
    R_c_B = states[10, 4, ],
    D_c_B = states[10, 5, ],
    
    S_sea_c_B = states[11, 1, ],
    E_sea_c_B = states[11, 2, ],
    I_sea_c_B = states[11, 3, ],
    R_sea_c_B = states[11, 4, ],
    D_sea_c_B = states[11, 5, ],
    
    S_c_NB = states[12, 1, ],
    E_c_NB = states[12, 2, ],
    I_c_NB = states[12, 3, ],
    R_c_NB = states[12, 4, ],
    D_c_NB = states[12, 5, ],

    S_sea_NB = states[13, 1, ],
    E_sea_NB = states[13, 2, ],
    I_sea_NB = states[13, 3, ],
    R_sea_NB = states[13, 4, ],
    D_sea_NB = states[13, 5, ]

) %>% 
    mutate(
    S_a_B_total = S_a_B + S_sea_a_B,
    E_a_B_total = E_a_B+ E_sea_a_B,
    I_a_B_total = I_a_B + I_sea_a_B,
    R_a_B_total = R_a_B + R_sea_a_B,
    D_a_B_total = D_a_B + D_sea_a_B,
    
    S_b_B_total = S_b_B + S_sea_b_B,
    E_b_B_total = E_b_B + E_sea_b_B,
    I_b_B_total = I_b_B + I_sea_b_B,
    R_b_B_total = R_b_B + R_sea_b_B,
    D_b_B_total = D_b_B + D_sea_b_B,
    
    S_c_B_total = S_c_B + S_sea_c_B,
    E_c_B_total = E_c_B + E_sea_c_B,
    I_c_B_total = I_c_B + I_sea_c_B,
    R_c_B_total = R_c_B + R_sea_c_B,
    D_c_B_total = D_c_B + D_sea_c_B
  )
  
  output = list(all_states,
       simulated_dispersal_date)
  
  return(output)
} # function


# Plot results

plot_seir = function(output_){
  
  output_1 = output_[[1]]
  output_2 = output_[[2]]
  
  
  output_long = melt(output_1[1:nrow(output_1)-1, ], id = "time")
  
  # In A
  output_a = output_long %>% filter(variable %in% c("S_a_B_total", "E_a_B_total", "I_a_B_total", "R_a_B_total", "D_a_B_total"))
  plot_a = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 , linetype= if (is.na(output_2)) "blank" else "dashed", 
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_a, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Breeders in A (colony+sea)")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_a$value == 0)) 1 else NA) 
  
  output_a_N = output_long %>% filter(variable %in% c("S_a_N", "E_a_N", "I_a_N", "R_a_N", "D_a_N"))
  plot_a_N = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 , linetype= if (is.na(output_2)) "blank" else "dashed",
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_a_N, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Nestlings in A")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_a_N$value == 0)) 1 else NA)
  
  output_a_NB = output_long %>% filter(variable %in% c("S_a_NB", "E_a_NB", "I_a_NB", "R_a_NB", "D_a_NB"))
  plot_a_NB =  ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 ,linetype= if (is.na(output_2)) "blank" else "dashed",
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_a_NB, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Non-breeder in colony A")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_a_NB$value == 0)) 1 else NA)
  
  # In B
  output_b = output_long %>% filter(variable %in% c("S_b_B_total", "E_b_B_total", "I_b_B_total", "R_b_B_total", "D_b_B_total"))
  plot_b = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 ,linetype= if (is.na(output_2)) "blank" else "dashed",
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_b, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Breeders in B (colony+sea)")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_b$value == 0)) 1 else NA)
  
  output_b_N = output_long %>% filter(variable %in% c("S_b_N", "E_b_N", "I_b_N", "R_b_N", "D_b_N"))
  plot_b_N = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 ,linetype= if (is.na(output_2)) "blank" else "dashed",
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_b_N, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Nestlings in B")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_b_N$value == 0)) 1 else NA)
  
  output_b_NB = output_long %>% filter(variable %in% c("S_b_NB", "E_b_NB", "I_b_NB", "R_b_NB", "D_b_NB"))
  plot_b_NB = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 ,linetype= if (is.na(output_2)) "blank" else "dashed", 
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_b_NB, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Non-breeder in colony B")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_b_NB$value == 0)) 1 else NA)
  
  # In C
  
  output_c = output_long %>% filter(variable %in% c("S_c_B_total", "E_c_B_total", "I_c_B_total", "R_c_B_total", "D_c_B_total"))
  plot_c = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 , linetype= if (is.na(output_2)) "blank" else "dashed", 
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_c, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Breeders in C (colony+sea)")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_c$value == 0)) 1 else NA)
  
  output_c_N = output_long %>% filter(variable %in% c("S_c_N", "E_c_N", "I_c_N", "R_c_N", "D_c_N"))
  plot_c_N = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 ,linetype= if (is.na(output_2)) "blank" else "dashed",
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_c_N, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Nestlings in C")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_c_N$value == 0)) 1 else NA)
  
  output_c_NB = output_long %>% filter(variable %in% c("S_c_NB", "E_c_NB", "I_c_NB", "R_c_NB", "D_c_NB"))
  plot_c_NB = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 ,linetype= if (is.na(output_2)) "blank" else "dashed",
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_c_NB, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Non-breeder in colony C")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_c_NB$value == 0)) 1 else NA)
  
  
  # At sea
  
  output_sea_NB = output_long %>% filter(variable %in% c("S_sea_NB", "E_sea_NB", "I_sea_NB", "R_sea_NB", "D_sea_NB"))
  plot_sea_NB = ggplot() +
    geom_vline(xintercept = if (is.na(output_2)) 1 else output_2 , linetype= if (is.na(output_2)) "blank" else "dashed", 
               color = if (is.na(output_2)) "white" else "grey", linewidth=0.65) +
    geom_line(data = output_sea_NB, aes(x = time, y = value, color = variable)) +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Non-breeder at sea")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))+
    ylim(0, if (all(output_sea_NB$value == 0)) 1 else NA)
  
  
  plot_grid_seir = plot_grid(plot_a, plot_a_N, plot_sea_NB,
                             plot_b, plot_b_N, plot_b_NB,
                             plot_c, plot_c_N, plot_c_NB,
                             labels = c("A", "B", "C",
                                        "D", "E", "F",
                                        "G", "H", "I"),
                             label_size = 12)
  
  print(plot_grid_seir)
  
}

# Summary output ----------------------------------------------------------

summary_output = function(# Time series of all states
                          output, 
                          # Probability of a nestling becoming a breeder
                          reaching.repro.prob = 0.3){
  
  output = output[[1]]

  N_a = output[1, c("S_a_B", "I_a_B", "S_sea_a_B", "I_sea_a_B")] %>% sum()
  dead_a = output[nrow(output), c("D_a_B","D_sea_a_B")] %>% sum()
  a_N = output[nrow(output), c("S_a_N", "E_a_N", "I_a_N", "R_a_N")] %>% sum()
  max_infected_a =  output[, c("E_a_B", "E_sea_a_B", "I_a_B","I_sea_a_B")] %>% 
    rowSums() %>% 
    max()
  
 
  N_b = output[1, c("S_b_B", "I_b_B", "S_sea_b_B", "I_sea_b_B")] %>% sum()
  dead_b = output[nrow(output), c("D_b_B","D_sea_b_B")] %>% sum()
  b_N = output[nrow(output), c("S_b_N", "E_b_N", "I_b_N", "R_b_N")] %>% sum()
  max_infecteD_b_B =  output[, c("E_b_B", "E_sea_b_B", "I_b_B","I_sea_b_B")] %>% 
    rowSums() %>% 
    max()
  
  N_c = output[1, c("S_c_B", "I_c_B", "S_sea_c_B", "I_sea_c_B")] %>% sum()
  dead_c = output[nrow(output), c("D_c_B","D_sea_c_B")] %>% sum()
  c_N = output[nrow(output), c("S_c_N", "E_c_N", "I_c_N", "R_c_N")] %>% sum()
  max_infecteD_c_B =  output[, c("E_c_B", "E_sea_c_B", "I_c_B","I_sea_c_B")] %>% 
    rowSums() %>% 
    max()


  nb_adults = N_a + N_b + N_c - dead_a - dead_b - dead_c
  nb_nestlings = a_N + b_N + c_N
  nb_adults_equi = nb_adults + reaching.repro.prob * nb_nestlings
  
  
  nb_infected_colonies = 
    sum(max_infected_a > 0,
        max_infecteD_b_B > 0,
        max_infecteD_c_B > 0)
  
  infected_X_time = 0
  for (t in 1:(length(output$time)-1)){
     
    infected_X_time = infected_X_time + 
      output$I_sea_NB[t] * (output$time[t+1] - output$time[t])
    
  }

  return( data.frame(

    nb_adults_equi = nb_adults_equi,
    nb_infected_colonies = nb_infected_colonies,
    infected_X_time = infected_X_time

  ))
}

# Wrapper function --------------------------------------------------------

model_wrapper = function(  
    # Parameter of the taul-leap agorithm
  tau_ = 0.05,
  # Number of simu_adultlation days
  total_time_ = 50,
  # Do we induce dispersion ?
  induced_dispersal_ = F,
  # Induced dispersion mode (deterministic or stochastic)
  dispersal_stochastic_ = T,
  # Reaction time between 1rst death and induced dispersal
  dispersal_reaction_time_ = 4,
  # Initial conditions
  initial_number_infected_breeders_A_ = 3,
  initial_number_breeders_A_ = 100,
  initial_number_breeders_B_ = 80,
  initial_number_breeders_C_ = 20,
  # Transmission rate from exposed individuals and from infectious individuals in a colony
  BETA_ = 0.02,
  # Time at sea before returning to a colony (non-breeders)
  TIME_AT_SEA_NB_ = 40,
  # Probability of a nestling becoming a breeder
  reaching.repro.prob_ = 0.3){
  
  output_ = gillespie_seir(
    tau = tau_,
    total_time = total_time_,
    induced_dispersal = induced_dispersal_,
    dispersal_stochastic = dispersal_stochastic_,
    dispersal_reaction_time = dispersal_reaction_time_,
    initial_number_infected_breeders_A = initial_number_infected_breeders_A_,
    initial_number_breeders_A = initial_number_breeders_A_,
    initial_number_breeders_B = initial_number_breeders_B_,
    initial_number_breeders_C = initial_number_breeders_C_,
    BETA = BETA_,
    TIME_AT_SEA_NB = TIME_AT_SEA_NB_)
  
  res = summary_output(output = output_, reaching.repro.prob = reaching.repro.prob_)
  
  return(res)
}




# model_wrapper(  # Parameter of the taul-leap agorithm
#   tau_ = 0.05,
#   # Number of simu_adultlation days
#   total_time_ = 50,
#   # Do we induce dispersion ?
#   induced_dispersal_ = F,
#   # Induced dispersion mode (deterministic or stochastic)
#   dispersal_stochastic_ = T,
#   # Reaction time between 1rst death and induced dispersal
#   dispersal_reaction_time_ = 4,
#   # Initial conditions
#   initial_number_infected_breeders_A_ = 3,
#   initial_number_breeders_A_ = 100,
#   initial_number_breeders_B_ = 80,
#   initial_number_breeders_C_ = 20,
#   # Transmission rate from exposed individuals and from infectious individuals in a colony
#   BETA_ = 0.02,
#   # Time at sea before returning to a colony (non-breeders)
#   TIME_AT_SEA_NB_ = 40,
#   # Probability of a nestling becoming a breeder
#   reaching.repro.prob_ = 0.3)


# Run simulation ----------------------------------------------------------

time1 <- Sys.time()
output = gillespie_seir(# Parameter of the taul-leap agorithm
  tau = 0.25,
  # Number of simu_adultlation days
  total_time = 70,
  # Initial conditions
  initial_number_infected_breeders_A = 1,
  initial_number_breeders_A = 50,
  initial_number_breeders_B = 50,
  initial_number_breeders_C = 50,
  # Induced dispersion parameters
  # Do we induce dispersion ?
  induced_dispersal = F,
  # Induced dispersion mode (deterministic or stochastic)
  dispersal_stochastic = F
  )
time2 <- Sys.time()
time2 - time1

plot_seir(output_ = output)






