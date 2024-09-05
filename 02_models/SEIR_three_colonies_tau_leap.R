

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
                             S_a, E_a, I_a, R_a, D_a,
                             S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                             S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                             # B
                             S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                             S_b, E_b, I_b, R_b, D_b,
                             S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                             S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                             # C
                             S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                             S_c, E_c, I_c, R_c, D_c,
                             S_sea_c, E_sea_c, I_sea_c, R_sea_c, D_sea_c,
                             S_c_NB, E_c_NB, I_c_NB, R_c_NB, D_c_NB,
                             # Non-breeders at sea
                             S_sea_NB, E_sea_NB, I_sea_NB, R_sea_NB, D_sea_NB
                             
                             )
                            
                             {
  
  rates = c(
    
  # SEIR transitions
  
  ## Nestlings
  ### Colony A
  "S_a_N_to_E_a_N" = beta_E_colony * S_a_N * (E_a+E_a_NB+E_a_N) + 
                     beta_I_colony * S_a_N * (I_a+I_a_NB+I_a_N),
  "E_a_N_to_S_a_N" = eta * E_a_N,
  "E_a_N_to_I_a_N" = sigma * E_a_N,
  "I_a_N_to_R_a_N" = gamma * I_a_N,
  "I_a_N_to_D_a_N" = mu_nestling * I_a_N,
  ### Colony B
  "S_b_N_to_E_b_N" = beta_E_colony * S_b_N * (E_b+E_b_NB+E_b_N) +
                     beta_I_colony * S_b_N * (I_b+I_b_NB+I_b_N),
  "E_b_N_to_S_b_N" = eta * E_b_N,
  "E_b_N_to_I_b_N" = sigma * E_b_N,
  "I_b_N_to_R_b_N" = gamma * I_b_N,
  "I_b_N_to_D_b_N" = mu_nestling * I_b_N,
  ### Colony C
  "S_c_N_to_E_c_N" = beta_E_colony * S_c_N * (E_c+E_c_NB+E_c_N) +
                     beta_I_colony * S_c_N * (I_c+I_c_NB+I_c_N),
  "E_c_N_to_S_c_N" = eta * E_c_N,
  "E_c_N_to_I_c_N" = sigma * E_c_N,
  "I_c_N_to_R_c_N" = gamma * I_c_N,
  "I_c_N_to_D_c_N" = mu_nestling * I_c_N,
  
  ## Non-Breeders
  ### Colony A
  "S_a_NB_to_E_a_NB" = beta_E_colony * S_a_NB * (E_a+E_a_NB+E_a_N) +
                       beta_I_colony * S_a_NB * (I_a+I_a_NB+I_a_N),
  "E_a_NB_to_S_a_NB" = eta * E_a_NB,
  "E_a_NB_to_I_a_NB" = sigma * E_a_NB,
  "I_a_NB_to_R_a_NB" = gamma * I_a_NB,
  "I_a_NB_to_D_a_NB" = mu_adult * I_a_NB,
  ### Colony B
  "S_b_NB_to_E_b_NB" = beta_E_colony * S_b_NB * (E_b+E_b_NB+E_b_N) + 
                       beta_I_colony * S_b_NB * (I_b+I_b_NB+I_b_N),
  "E_b_NB_to_S_b_NB" = eta * E_b_NB,
  "E_b_NB_to_I_b_NB" = sigma * E_b_NB,
  "I_b_NB_to_R_b_NB" = gamma * I_b_NB,
  "I_b_NB_to_D_b_NB" = mu_adult * I_b_NB,
  ### Colony C
  "S_c_NB_to_E_c_NB" = beta_E_colony * S_c_NB * (E_c+E_c_NB+E_c_N) + 
                       beta_I_colony * S_c_NB * (I_c+I_c_NB+I_c_N),
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
  "S_a_to_E_a" = beta_E_colony * S_a * (E_a+E_a_NB+E_a_N) + 
                 beta_I_colony * S_a * (I_a+I_a_NB+I_a_N),
  "E_a_to_S_a" = eta * E_a,
  "E_a_to_I_a" = sigma * E_a,
  "I_a_to_R_a" = gamma * I_a,
  "I_a_to_D_a" = mu_adult * I_a,
  ### Sea A
  "S_sea_a_to_E_sea_a" = 0,
  "E_sea_a_to_S_sea_a" = eta * E_sea_a,
  "E_sea_a_to_I_sea_a" = sigma * E_sea_a,
  "I_sea_a_to_R_sea_a" = gamma * I_sea_a,
  "I_sea_a_to_D_sea_a" = mu_adult * I_sea_a,
  ### Colony B
  "S_b_to_E_b" = beta_E_colony * S_b * (E_b+E_b_NB+E_b_N) + 
                 beta_I_colony * S_b * (I_b+I_b_NB+I_b_N),
  "E_b_to_S_b" = eta * E_b,
  "E_b_to_I_b" = sigma * E_b,
  "I_b_to_R_b" = gamma * I_b,
  "I_b_to_D_b" = mu_adult * I_b,
  ### Sea B
  "S_sea_b_to_E_sea_b" = 0,
  "E_sea_b_to_S_sea_b" = eta * E_sea_b,
  "E_sea_b_to_I_sea_b" = sigma * E_sea_b,
  "I_sea_b_to_R_sea_b" = gamma * I_sea_b,
  "I_sea_b_to_D_sea_b" = mu_adult * I_sea_b,
  ### Colony C
  "S_c_to_E_c" = beta_E_colony * S_c * (E_c+E_c_NB+E_c_N) + 
                 beta_I_colony * S_c * (I_c+I_c_NB+I_c_N),
  "E_c_to_S_c" = eta * E_c,
  "E_c_to_I_c" = sigma * E_c,
  "I_c_to_R_c" = gamma * I_c,
  "I_c_to_D_c" = mu_adult * I_c,
  ### Sea C
  "S_sea_c_to_E_sea_c" = 0,
  "E_sea_c_to_S_sea_c" = eta * E_sea_c,
  "E_sea_c_to_I_sea_c" = sigma * E_sea_c,
  "I_sea_c_to_R_sea_c" = gamma * I_sea_c,
  "I_sea_c_to_D_sea_c" = mu_adult * I_sea_c,

  
  
  # Mobility
  ## Non-Breeders
  ### In A
  ### From colony A to sea 
  "S_a_NB_to_S_sea_NB" = rho_to_sea * S_a_NB,
  "E_a_NB_to_E_sea_NB" = rho_to_sea * E_a_NB,
  "I_a_NB_to_I_sea_NB" = rho_to_sea * I_a_NB,
  "R_a_NB_to_R_sea_NB" = rho_to_sea * R_a_NB,
  ### From sea to colony A (conspecific attraction)
  "S_sea_NB_to_S_a_NB" = rho_to_colony * S_sea_NB * ( (S_a+E_a+I_a+R_a)
                                                     /(S_a+E_a+I_a+R_a+
                                                       S_b+E_b+I_b+R_b+
                                                       S_c+E_c+I_c+R_c+1)
                                                     ),
  
  "E_sea_NB_to_E_a_NB" = rho_to_colony * E_sea_NB * ( (S_a+E_a+I_a+R_a)
                                                     /(S_a+E_a+I_a+R_a+
                                                       S_b+E_b+I_b+R_b+
                                                       S_c+E_c+I_c+R_c+1)
                                                     ),
  "I_sea_NB_to_I_a_NB" = rho_to_colony * I_sea_NB * ( (S_a+E_a+I_a+R_a)
                                                     /(S_a+E_a+I_a+R_a+
                                                       S_b+E_b+I_b+R_b+
                                                       S_c+E_c+I_c+R_c+1)
                                                     ),
  "R_sea_NB_to_R_a_NB" = rho_to_colony * R_sea_NB * ( (S_a+E_a+I_a+R_a)
                                                      /(S_a+E_a+I_a+R_a+
                                                        S_b+E_b+I_b+R_b+
                                                        S_c+E_c+I_c+R_c+1)
                                                      ),
  ### In B
  ### From colony B to sea
  "S_b_NB_to_S_sea_NB" = rho_to_sea * S_b_NB,
  "E_b_NB_to_E_sea_NB" = rho_to_sea * E_b_NB,
  "I_b_NB_to_I_sea_NB" = rho_to_sea * I_b_NB,
  "R_b_NB_to_R_sea_NB" = rho_to_sea * R_b_NB,
  ### From sea to colony B
  "S_sea_NB_to_S_b_NB" = rho_to_colony * S_sea_NB * ( (S_b+E_b+I_b+R_b)
                                                     /(S_a+E_a+I_a+R_a+
                                                       S_b+E_b+I_b+R_b+
                                                       S_c+E_c+I_c+R_c+1)
                                                     ),
  "E_sea_NB_to_E_b_NB" = rho_to_colony * E_sea_NB * ( (S_b+E_b+I_b+R_b)
                                                     /(S_a+E_a+I_a+R_a+
                                                       S_b+E_b+I_b+R_b+
                                                       S_c+E_c+I_c+R_c+1)
                                                     ),
  "I_sea_NB_to_I_b_NB" = rho_to_colony * I_sea_NB * ( (S_b+E_b+I_b+R_b)
                                                     /(S_a+E_a+I_a+R_a+
                                                       S_b+E_b+I_b+R_b+
                                                       S_c+E_c+I_c+R_c+1)
                                                     ),
  "R_sea_NB_to_R_b_NB" = rho_to_colony * R_sea_NB * ( (S_b+E_b+I_b+R_b)
                                                     /(S_a+E_a+I_a+R_a+
                                                       S_b+E_b+I_b+R_b+
                                                       S_c+E_c+I_c+R_c+1)
                                                     ),
  ### In C
  ### From colony C to sea
  "S_c_NB_to_S_sea_NB" = rho_to_sea * S_c_NB,
  "E_c_NB_to_E_sea_NB" = rho_to_sea * E_c_NB,
  "I_c_NB_to_I_sea_NB" = rho_to_sea * I_c_NB,
  "R_c_NB_to_R_sea_NB" = rho_to_sea * R_c_NB,
  ### From sea to colony C
  "S_sea_NB_to_S_c_NB" = rho_to_colony * S_sea_NB * ( (S_c+E_c+I_c+R_c)
                                                      /(S_a+E_a+I_a+R_a+
                                                        S_b+E_b+I_b+R_b+
                                                        S_c+E_c+I_c+R_c+1)
                                                      ),
  "E_sea_NB_to_E_c_NB" = rho_to_colony * E_sea_NB * ( (S_c+E_c+I_c+R_c)
                                                      /(S_a+E_a+I_a+R_a+
                                                        S_b+E_b+I_b+R_b+
                                                        S_c+E_c+I_c+R_c+1)
                                                      ),
  "I_sea_NB_to_I_c_NB" = rho_to_colony * I_sea_NB * ( (S_c+E_c+I_c+R_c)
                                                      /(S_a+E_a+I_a+R_a+
                                                        S_b+E_b+I_b+R_b+
                                                        S_c+E_c+I_c+R_c+1)
                                                      ),
  "R_sea_NB_to_R_c_NB" = rho_to_colony * R_sea_NB * ( (S_c+E_c+I_c+R_c)
                                                      /(S_a+E_a+I_a+R_a+
                                                        S_b+E_b+I_b+R_b+
                                                        S_c+E_c+I_c+R_c+1)
                                                      ),
  
  ## Breeders
  ### In A
  ### From colony A to sea A
  "S_a_to_S_sea_a" = zeta_to_sea * S_a,
  "E_a_to_E_sea_a" = zeta_to_sea * E_a,
  "I_a_to_I_sea_a" = zeta_to_sea * I_a,
  "R_a_to_R_sea_a" = zeta_to_sea * R_a,
  ### From sea A to colony A
  "S_sea_a_to_S_a" = zeta_to_colony * S_sea_a,
  "E_sea_a_to_E_a" = zeta_to_colony * E_sea_a,
  "I_sea_a_to_I_a" = zeta_to_colony * I_sea_a,
  "R_sea_a_to_R_a" = zeta_to_colony * R_sea_a,
  ### In B
  ### From colony B to sea B
  "S_b_to_S_sea_b" = zeta_to_sea * S_b,
  "E_b_to_E_sea_b" = zeta_to_sea * E_b,
  "I_b_to_I_sea_b" = zeta_to_sea * I_b,
  "R_b_to_R_sea_b" = zeta_to_sea * R_b,
  ### From sea B to colony B
  "S_sea_b_to_S_b" = zeta_to_colony * S_sea_b,
  "E_sea_b_to_E_b" = zeta_to_colony * E_sea_b,
  "I_sea_b_to_I_b" = zeta_to_colony * I_sea_b,
  "R_sea_b_to_R_b" = zeta_to_colony * R_sea_b,
  ### In C
  ### From colony B to sea C
  "S_c_to_S_sea_c" = zeta_to_sea * S_c,
  "E_c_to_E_sea_c" = zeta_to_sea * E_c,
  "I_c_to_I_sea_c" = zeta_to_sea * I_c,
  "R_c_to_R_sea_c" = zeta_to_sea * R_c,
  ### From sea B to colony C
  "S_sea_c_to_S_c" = zeta_to_colony * S_sea_c,
  "E_sea_c_to_E_c" = zeta_to_colony * E_sea_c,
  "I_sea_c_to_I_c" = zeta_to_colony * I_sea_c,
  "R_sea_c_to_R_c" = zeta_to_colony * R_sea_c
  
  # Breeders become Non-Breeders
  
  # "S_a_to_S_a_NB" = psi * S_a,
  # "E_a_to_E_a_NB" = psi * E_a,
  # "I_a_to_I_a_NB" = psi * I_a,
  # "R_a_to_R_a_NB" = psi * R_a
  
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
                          # Do we induce dispersion ?
                          induced_dispersal = T,
                          # Induced dispersion mode (deterministic or stochastic)
                          dispersal_stochactic = T,
                          # Reaction time between 1rst death and induced dispersal
                          dispersal_reaction_time = 5,
                          # Transmission rate from exposed individuals and from infectious individuals in a colony
                          BETA = 0.5,
                          # Time at sea before returning to a colony (non-breeders)
                          TIME_AT_SEA_NB = 40
                          ) {
  
  # Parameters --------------------------------------------------------------
  
  param = list( 
    
    # Epidemiological parameters
    
    ## Transmission rate from exposed individuals and from infectious individuals in a colony
    ## First row - in colonies / Second row - at sea // First column - for exposed / Second column - for infectious
    beta = matrix(c(0.00, BETA,
                    0.00, 0.00),
                  nrow = 2, ncol = 2, byrow = T),
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
    rho_to_colony = 1/TIME_AT_SEA_NB , 
    
    
    # Transition from breeder to non-breeder (reproductive failure)
    # psi = 0,  
    
    
    # Induced dispersion parameters
    
    ## Proportion of dispersed adults
    prop_dispersal = 1,
    ## Date of induced dispersion (if deterministic)
    dispersal_date = 0,
    
    # Demographic parameters
    
    ## Hatching date of the chicks
    hatching_date = 10
    
  )
  
  # Initial state
  
  ## Nestlings
  ## In colony A
  N = 0                
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_A_N = c(S = initial_susceptible,
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
  
  initial_state_B_N = c(S = initial_susceptible,
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
  
  initial_state_C_N = c(S = initial_susceptible,
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
  
  initial_state_A = c(S = initial_susceptible,
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
  
  initial_state_sea_a = c(S = initial_susceptible,
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
  
  initial_state_B <- c(S = initial_susceptible,
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
  
  initial_state_sea_b = c(S = initial_susceptible,
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
  
  initial_state_C <- c(S = initial_susceptible,
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
  
  initial_state_sea_c = c(S = initial_susceptible,
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
  
  initial_state_A_NB = c(S = initial_susceptible,
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
  
  initial_state_B_NB <- c(S = initial_susceptible,
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
  
  initial_state_C_NB <- c(S = initial_susceptible,
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
  
  

  
  
  initial_state = matrix(data = c(initial_state_A_N,
                                  initial_state_A,
                                  initial_state_sea_a,
                                  initial_state_A_NB,
                                  
                                  initial_state_B_N,
                                  initial_state_B,
                                  initial_state_sea_b,
                                  initial_state_B_NB,
                                  
                                  initial_state_C_N,
                                  initial_state_C,
                                  initial_state_sea_c,
                                  initial_state_C_NB,
  
                                  initial_state_sea_NB
                                  ), 
                         nrow = 13, ncol = 5, 
                         byrow = T)
  
  # Parameters
  
  beta_E_colony = param$beta[1,1]
  beta_I_colony = param$beta[1,2]
  
  sigma = param$sigma
  eta = param$eta
  gamma = param$gamma
  mu_adult = param$mu_adult
  mu_nestling = param$mu_nestling
  
  zeta_to_colony = param$zeta_to_colony
  zeta_to_sea = param$zeta_to_sea
  psi = param$psi
  rho_to_colony = param$rho_to_colony
  rho_to_sea = param$rho_to_sea
  
  prop_dispersal = param$prop_dispersal
  dispersal_date = param$dispersal_date
  
  hatching_date = param$hatching_date
  
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
    
    S_a = states[2, 1, dim(states)[3]]
    E_a = states[2, 2, dim(states)[3]]
    I_a = states[2, 3, dim(states)[3]]
    R_a = states[2, 4, dim(states)[3]]
    D_a = states[2, 5, dim(states)[3]]
    
    S_sea_a = states[3, 1, dim(states)[3]]
    E_sea_a = states[3, 2, dim(states)[3]]
    I_sea_a = states[3, 3, dim(states)[3]]
    R_sea_a = states[3, 4, dim(states)[3]]
    D_sea_a = states[3, 5, dim(states)[3]]
  
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
    
    S_b = states[6, 1, dim(states)[3]]
    E_b = states[6, 2, dim(states)[3]]
    I_b = states[6, 3, dim(states)[3]]
    R_b = states[6, 4, dim(states)[3]]
    D_b = states[6, 5, dim(states)[3]]
    
    S_sea_b = states[7, 1, dim(states)[3]]
    E_sea_b = states[7, 2, dim(states)[3]]
    I_sea_b = states[7, 3, dim(states)[3]]
    R_sea_b = states[7, 4, dim(states)[3]]
    D_sea_b = states[7, 5, dim(states)[3]]
    
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
    
    S_c = states[10, 1, dim(states)[3]]
    E_c = states[10, 2, dim(states)[3]]
    I_c = states[10, 3, dim(states)[3]]
    R_c = states[10, 4, dim(states)[3]]
    D_c = states[10, 5, dim(states)[3]]
    
    S_sea_c = states[11, 1, dim(states)[3]]
    E_sea_c = states[11, 2, dim(states)[3]]
    I_sea_c = states[11, 3, dim(states)[3]]
    R_sea_c = states[11, 4, dim(states)[3]]
    D_sea_c = states[11, 5, dim(states)[3]]
    
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
                            S_a, E_a, I_a, R_a, D_a,
                            S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                            S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                            # B
                            S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                            S_b, E_b, I_b, R_b, D_b,
                            S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                            S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                            # C
                            S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                            S_c, E_c, I_c, R_c, D_c,
                            S_sea_c, E_sea_c, I_sea_c, R_sea_c, D_sea_c,
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
      
      S_a_N = round((S_a + E_a + I_a + R_a + S_sea_a + E_sea_a + I_sea_a + R_sea_a)/2)
      S_b_N = round((S_b + E_b + I_b + R_b + S_sea_b + E_sea_b + I_sea_b + R_sea_b)/2)
      S_c_N = round((S_c + E_c + I_c + R_c + S_sea_c + E_sea_c + I_sea_c + R_sea_c)/2)
      
      already_hatched = T
      
      new_state = matrix(data = c(S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                  S_a, E_a, I_a, R_a, D_a,
                                  S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                  S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                  
                                  S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                  S_b, E_b, I_b, R_b, D_b,
                                  S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                  S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                  
                                  S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                  S_c, E_c, I_c, R_c, D_c,
                                  S_sea_c, E_sea_c, I_sea_c, R_sea_c, D_sea_c,
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
                                          S_a, E_a, I_a, R_a, D_a,
                                          S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                          S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                          # B
                                          S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                          S_b, E_b, I_b, R_b, D_b,
                                          S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                          S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                          # C
                                          S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                          S_c, E_c, I_c, R_c, D_c,
                                          S_sea_c, E_sea_c, I_sea_c, R_sea_c, D_sea_c,
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
    if (!first_death & D_a > 0){
      first_death_date = times[length(times)] 
      first_death = T
    }
    
    # Induction of dispersion
    ## Has the induced dispersal strategy been triggered?
    if (induced_dispersal){ 
      ## Has the dispersal occured?
      if (!already_dispersed){
        # Is it time to induce dispersion, according to the stochastic case and the deterministic case?
        if ((dispersal_stochactic & first_death & next_time > first_death_date + dispersal_reaction_time) 
            |
            (!dispersal_stochactic & next_time > dispersal_date)
        ){
          
          # Dispersal of breeders
          
          # Number of adults in  A
          N_a = S_a + E_a + I_a + R_a + S_sea_a + E_sea_a + I_sea_a + R_sea_a
          # Number of adults A who are dispersed
          N_disp_a = round(N_a * prop_dispersal)
          # Distribution of dispersed adults by epidemiological status
          disp_a = sample(c(rep("S_a", S_a), rep("E_a", E_a),rep("I_a", I_a),rep("R_a", R_a),
                            rep("S_sea_a", S_sea_a), rep("E_sea_a", E_sea_a),rep("I_sea_a", I_sea_a),rep("R_sea_a", R_sea_a)),
                          size = N_disp_a, 
                          replace = F) %>% 
            factor(., levels = c("S_a","E_a","I_a","R_a",
                                 "S_sea_a","E_sea_a","I_sea_a","R_sea_a")) %>% 
            table()
          
          # Number of susceptible adults who are dispersed from A
          disp_S_a = disp_a["S_a"]
          disp_S_sea_a = disp_a["S_sea_a"]
          # Update of the number of susceptible adults
          S_a = S_a - disp_S_a
          S_sea_a = S_sea_a - disp_S_sea_a
          S_sea_NB = S_sea_NB + disp_S_a + disp_S_sea_a
          
          # Number of exposed  adults who are dispersed from A
          disp_E_a = disp_a["E_a"]
          disp_E_sea_a = disp_a["E_sea_a"]
          # Update of the number of exposed adults
          E_a = E_a - disp_E_a
          E_sea_a = E_sea_a - disp_E_sea_a
          E_sea_NB = E_sea_NB + disp_E_a + disp_E_sea_a
          
          # Number of infectious adults who are dispersed from A
          disp_I_a = disp_a["I_a"]
          disp_I_sea_a = disp_a["I_sea_a"]
          # Update of the number of infectious adults
          I_a = I_a - disp_I_a
          I_sea_a = I_sea_a - disp_I_sea_a
          I_sea_NB = I_sea_NB + disp_I_a + disp_I_sea_a
          
          # Number of recovered adults who are dispersed from A
          disp_R_a = disp_a["R_a"]
          disp_R_sea_a = disp_a["R_sea_a"]
          # Update of the number of recovered adults
          R_a = R_a - disp_R_a
          R_sea_a = R_sea_a - disp_R_sea_a
          R_sea_NB = R_sea_NB + disp_R_a + disp_R_sea_a
          
          # Death of nestlings
          
          D_a_N = D_a_N + (S_a_N + E_a_N + I_a_N + R_a_N)
          S_a_N = 0
          E_a_N = 0
          I_a_N = 0
          R_a_N = 0
          
          
          already_dispersed = T
          
          new_state = matrix(data = c(S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                      S_a, E_a, I_a, R_a, D_a,
                                      S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                      S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                      
                                      S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                      S_b, E_b, I_b, R_b, D_b,
                                      S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                      S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                      
                                      S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                      S_c, E_c, I_c, R_c, D_c,
                                      S_sea_c, E_sea_c, I_sea_c, R_sea_c, D_sea_c,
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
                                              S_a, E_a, I_a, R_a, D_a,
                                              S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                              S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                              # B
                                              S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                              S_b, E_b, I_b, R_b, D_b,
                                              S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                              S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                              # C
                                              S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                              S_c, E_c, I_c, R_c, D_c,
                                              S_sea_c, E_sea_c, I_sea_c, R_sea_c, D_sea_c,
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
      parent1 = sample(c(rep("S_a", S_a), rep("E_a", E_a),rep("I_a", I_a),rep("R_a", R_a),
                         rep("S_sea_a", S_sea_a), rep("E_sea_a", E_sea_a),rep("I_sea_a", I_sea_a),rep("R_sea_a", R_sea_a)),
                       size = 1)
      if (parent1 == "S_a"){
        S_a = S_a - 1
        S_a_NB = S_a_NB + 1
      } else if (parent1 == "E_a"){
        E_a = E_a - 1
        E_a_NB = E_a_NB + 1
      } else if (parent1 == "I_a"){
        I_a = I_a - 1
        I_a_NB = I_a_NB + 1
      } else if (parent1 == "R_a"){
        R_a = R_a - 1
        R_a_NB = R_a_NB + 1
      } else if (parent1 == "S_sea_a"){
        S_sea_a = S_sea_a - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent1 == "E_sea_a"){
        E_sea_a = E_sea_a - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent1 == "I_sea_a"){
        I_sea_a = I_sea_a - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent1 == "R_sea_a"){
        R_sea_a = R_sea_a - 1
        R_sea_NB = R_sea_NB + 1
      }
      parent2 = sample(c(rep("S_a", S_a), rep("E_a", E_a),rep("I_a", I_a),rep("R_a", R_a),
                         rep("S_sea_a", S_sea_a), rep("E_sea_a", E_sea_a),rep("I_sea_a", I_sea_a),rep("R_sea_a", R_sea_a)),
                       size = 1)
      if (parent2 == "S_a"){
        S_a = S_a - 1
        S_a_NB = S_a_NB + 1
      } else if (parent2 == "E_a"){
        E_a = E_a - 1
        E_a_NB = E_a_NB + 1
      } else if (parent2 == "I_a"){
        I_a = I_a - 1
        I_a_NB = I_a_NB + 1
      } else if (parent2 == "R_a"){
        R_a = R_a - 1
        R_a_NB = R_a_NB + 1
      } else if (parent2 == "S_sea_a"){
        S_sea_a = S_sea_a - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent2 == "E_sea_a"){
        E_sea_a = E_sea_a - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent2 == "I_sea_a"){
        I_sea_a = I_sea_a - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent2 == "R_sea_a"){
        R_sea_a = R_sea_a - 1
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
      parent1 = sample(c(rep("S_b", S_b), rep("E_b", E_b),rep("I_b", I_b),rep("R_b", R_b),
                         rep("S_sea_b", S_sea_b), rep("E_sea_b", E_sea_b),rep("I_sea_b", I_sea_b),rep("R_sea_b", R_sea_b)),
                       size = 1)
      if (parent1 == "S_b"){
        S_b = S_b - 1
        S_b_NB = S_b_NB + 1
      } else if (parent1 == "E_b"){
        E_b = E_b - 1
        E_b_NB = E_b_NB + 1
      } else if (parent1 == "I_b"){
        I_b = I_b - 1
        I_b_NB = I_b_NB + 1
      } else if (parent1 == "R_b"){
        R_b = R_b - 1
        R_b_NB = R_b_NB + 1
      } else if (parent1 == "S_sea_b"){
        S_sea_b = S_sea_b - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent1 == "E_sea_b"){
        E_sea_b = E_sea_b - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent1 == "I_sea_b"){
        I_sea_b = I_sea_b - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent1 == "R_sea_b"){
        R_sea_b = R_sea_b - 1
        R_sea_NB = R_sea_NB + 1
      }
      parent2 = sample(c(rep("S_b", S_b), rep("E_b", E_b),rep("I_b", I_b),rep("R_b", R_b),
                         rep("S_sea_b", S_sea_b), rep("E_sea_b", E_sea_b),rep("I_sea_b", I_sea_b),rep("R_sea_b", R_sea_b)),
                       size = 1)
      if (parent2 == "S_b"){
        S_b = S_b - 1
        S_b_NB = S_b_NB + 1
      } else if (parent2 == "E_b"){
        E_b = E_b - 1
        E_b_NB = E_b_NB + 1
      } else if (parent2 == "I_b"){
        I_b = I_b - 1
        I_b_NB = I_b_NB + 1
      } else if (parent2 == "R_b"){
        R_b = R_b - 1
        R_b_NB = R_b_NB + 1
      } else if (parent2 == "S_sea_b"){
        S_sea_b = S_sea_b - 1
        S_sea_NB = S_sea_NB + 1
      } else if (parent2 == "E_sea_b"){
        E_sea_b = E_sea_b - 1
        E_sea_NB = E_sea_NB + 1
      } else if (parent2 == "I_sea_b"){
        I_sea_b = I_sea_b - 1
        I_sea_NB = I_sea_NB + 1
      } else if (parent2 == "R_sea_b"){
        R_sea_b = R_sea_b - 1
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
        parent1 = sample(c(rep("S_c", S_c), rep("E_c", E_c),rep("I_c", I_c),rep("R_c", R_c),
                           rep("S_sea_c", S_sea_c), rep("E_sea_c", E_sea_c),rep("I_sea_c", I_sea_c),rep("R_sea_c", R_sea_c)),
                         size = 1)
        if (parent1 == "S_c"){
          S_c = S_c - 1
          S_c_NB = S_c_NB + 1
        } else if (parent1 == "E_c"){
          E_c = E_c - 1
          E_c_NB = E_c_NB + 1
        } else if (parent1 == "I_c"){
          I_c = I_c - 1
          I_c_NB = I_c_NB + 1
        } else if (parent1 == "R_c"){
          R_c = R_c - 1
          R_c_NB = R_c_NB + 1
        } else if (parent1 == "S_sea_c"){
          S_sea_c = S_sea_c - 1
          S_sea_NB = S_sea_NB + 1
        } else if (parent1 == "E_sea_c"){
          E_sea_c = E_sea_c - 1
          E_sea_NB = E_sea_NB + 1
        } else if (parent1 == "I_sea_c"){
          I_sea_c = I_sea_c - 1
          I_sea_NB = I_sea_NB + 1
        } else if (parent1 == "R_sea_c"){
          R_sea_c = R_sea_c - 1
          R_sea_NB = R_sea_NB + 1
        }
        parent2 = sample(c(rep("S_c", S_c), rep("E_c", E_c),rep("I_c", I_c),rep("R_c", R_c),
                           rep("S_sea_c", S_sea_c), rep("E_sea_c", E_sea_c),rep("I_sea_c", I_sea_c),rep("R_sea_c", R_sea_c)),
                         size = 1)
        if (parent2 == "S_c"){
          S_c = S_c - 1
          S_c_NB = S_c_NB + 1
        } else if (parent2 == "E_c"){
          E_c = E_c - 1
          E_c_NB = E_c_NB + 1
        } else if (parent2 == "I_c"){
          I_c = I_c - 1
          I_c_NB = I_c_NB + 1
        } else if (parent2 == "R_c"){
          R_c = R_c - 1
          R_c_NB = R_c_NB + 1
        } else if (parent2 == "S_sea_c"){
          S_sea_c = S_sea_c - 1
          S_sea_NB = S_sea_NB + 1
        } else if (parent2 == "E_sea_c"){
          E_sea_c = E_sea_c - 1
          E_sea_NB = E_sea_NB + 1
        } else if (parent2 == "I_sea_c"){
          I_sea_c = I_sea_c - 1
          I_sea_NB = I_sea_NB + 1
        } else if (parent2 == "R_sea_c"){
          R_sea_c = R_sea_c - 1
          R_sea_NB = R_sea_NB + 1
        }
      ## Breeders
      ### In A
      #### At colony
      } else if  (transition == "S_a_to_E_a" & S_a > 0) {
        S_a = S_a - 1
        E_a = E_a + 1
      } else if (transition == "E_a_to_S_a" & E_a > 0) {
        E_a = E_a - 1
        S_a = S_a + 1
      } else if (transition == "E_a_to_I_a" & E_a > 0) {
        E_a = E_a - 1
        I_a = I_a + 1
      } else if (transition == "I_a_to_R_a" & I_a > 0) {
        I_a = I_a - 1
        R_a = R_a + 1
      } else if (transition == "I_a_to_D_a" & I_a > 0) {
        # If an adult dies, the partner becomes a non-breeder and the nestling dies
        I_a = I_a - 1
        D_a = D_a + 1
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
        partner = sample(c(rep("S_a", S_a), rep("E_a", E_a),rep("I_a", I_a),rep("R_a", R_a),
                           rep("S_sea_a", S_sea_a), rep("E_sea_a", E_sea_a),rep("I_sea_a", I_sea_a),rep("R_sea_a", R_sea_a)),
                         size = 1)
        if (partner == "S_a"){
          S_a = S_a - 1
          S_a_NB = S_a_NB + 1
        } else if (partner == "E_a"){
          E_a = E_a - 1
          E_a_NB = E_a_NB + 1
        } else if (partner == "I_a"){
          I_a = I_a - 1
          I_a_NB = I_a_NB + 1
        } else if (partner == "R_a"){
          R_a = R_a - 1
          R_a_NB = R_a_NB + 1
        } else if (partner == "S_sea_a"){
          S_sea_a = S_sea_a - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_a"){
          E_sea_a = E_sea_a - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_a"){
          I_sea_a = I_sea_a - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_a"){
          R_sea_a = R_sea_a - 1
          R_sea_NB = R_sea_NB + 1
        }
      #### At sea
      } else if (transition == "S_sea_a_to_E_sea_a" & S_sea_a > 0) {
        S_sea_a = S_sea_a - 1
        E_sea_a = E_sea_a + 1
      } else if (transition == "E_sea_a_to_S_sea_a" & E_sea_a > 0) {
        E_sea_a = E_sea_a - 1
        S_sea_a = S_sea_a + 1
      } else if (transition == "E_sea_a_to_I_sea_a" & E_sea_a > 0) {
        E_sea_a = E_sea_a - 1
        I_sea_a = I_sea_a + 1
      } else if (transition == "I_sea_a_to_R_sea_a" & I_sea_a > 0) {
        I_sea_a = I_sea_a - 1
        R_sea_a = R_sea_a + 1
      } else if (transition == "I_sea_a_to_D_sea_a" & I_sea_a > 0) {
        I_sea_a = I_sea_a - 1
        D_sea_a = D_sea_a + 1
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
        partner = sample(c(rep("S_a", S_a), rep("E_a", E_a),rep("I_a", I_a),rep("R_a", R_a),
                           rep("S_sea_a", S_sea_a), rep("E_sea_a", E_sea_a),rep("I_sea_a", I_sea_a),rep("R_sea_a", R_sea_a)),
                         size = 1)
        if (partner == "S_a"){
          S_a = S_a - 1
          S_a_NB = S_a_NB + 1
        } else if (partner == "E_a"){
          E_a = E_a - 1
          E_a_NB = E_a_NB + 1
        } else if (partner == "I_a"){
          I_a = I_a - 1
          I_a_NB = I_a_NB + 1
        } else if (partner == "R_a"){
          R_a = R_a - 1
          R_a_NB = R_a_NB + 1
        } else if (partner == "S_sea_a"){
          S_sea_a = S_sea_a - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_a"){
          E_sea_a = E_sea_a - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_a"){
          I_sea_a = I_sea_a - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_a"){
          R_sea_a = R_sea_a - 1
          R_sea_NB = R_sea_NB + 1
        }
      ### In B
      #### At colony
      } else if (transition == "S_b_to_E_b" & S_b > 0) {
        S_b = S_b - 1
        E_b = E_b + 1
      } else if (transition == "E_b_to_S_b" & E_b > 0) {
        E_b = E_b - 1
        S_b = S_b + 1
      } else if (transition == "E_b_to_I_b" & E_b > 0) {
        E_b = E_b - 1
        I_b = I_b + 1
      } else if (transition == "I_b_to_R_b" & I_b > 0) {
        I_b = I_b - 1
        R_b = R_b + 1
      } else if (transition == "I_b_to_D_b" & I_b > 0) {
        # If an adult dies, the partner becomes a non-breeder and the nestling dies
        I_b = I_b - 1
        D_b = D_b + 1
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
        partner = sample(c(rep("S_b", S_b), rep("E_b", E_b),rep("I_b", I_b),rep("R_b", R_b),
                           rep("S_sea_b", S_sea_b), rep("E_sea_b", E_sea_b),rep("I_sea_b", I_sea_b),rep("R_sea_b", R_sea_b)),
                         size = 1)
       
        if (partner == "S_b"){
          S_b = S_b - 1
          S_b_NB = S_b_NB + 1
        } else if (partner == "E_b"){
          E_b = E_b - 1
          E_b_NB = E_b_NB + 1
        } else if (partner == "I_b"){
          I_b = I_b - 1
          I_b_NB = I_b_NB + 1
        } else if (partner == "R_b"){
          R_b = R_b - 1
          R_b_NB = R_b_NB + 1
        } else if (partner == "S_sea_b"){
          S_sea_b = S_sea_b - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_b"){
          E_sea_b = E_sea_b - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_b"){
          I_sea_b = I_sea_b - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_b"){
          R_sea_b = R_sea_b - 1
          R_sea_NB = R_sea_NB + 1
        }
      #### At sea
      } else if (transition == "S_sea_b_to_E_sea_b" & S_sea_b > 0) {
        S_sea_b = S_sea_b - 1
        E_sea_b = E_sea_b + 1
      } else if (transition == "E_sea_b_to_S_sea_b" & E_sea_b > 0) {
        E_sea_b = E_sea_b - 1
        S_sea_b = S_sea_b + 1
      } else if (transition == "E_sea_b_to_I_sea_b" & E_sea_b > 0) {
        E_sea_b = E_sea_b - 1
        I_sea_b = I_sea_b + 1
      } else if (transition == "I_sea_b_to_R_sea_b" & I_sea_b > 0) {
        I_sea_b = I_sea_b - 1
        R_sea_b = R_sea_b + 1
      } else if (transition == "I_sea_b_to_D_sea_b" & I_sea_b > 0) {
        I_sea_b = I_sea_b - 1
        D_sea_b = D_sea_b + 1
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
        partner = sample(c(rep("S_b", S_b), rep("E_b", E_b),rep("I_b", I_b),rep("R_b", R_b),
                           rep("S_sea_b", S_sea_b), rep("E_sea_b", E_sea_b),rep("I_sea_b", I_sea_b),rep("R_sea_b", R_sea_b)),
                         size = 1)
        
        if (partner == "S_b"){
          S_b = S_b - 1
          S_b_NB = S_b_NB + 1
        } else if (partner == "E_b"){
          E_b = E_b - 1
          E_b_NB = E_b_NB + 1
        } else if (partner == "I_b"){
          I_b = I_b - 1
          I_b_NB = I_b_NB + 1
        } else if (partner == "R_b"){
          R_b = R_b - 1
          R_b_NB = R_b_NB + 1
        } else if (partner == "S_sea_b"){
          S_sea_b = S_sea_b - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_b"){
          E_sea_b = E_sea_b - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_b"){
          I_sea_b = I_sea_b - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_b"){
          R_sea_b = R_sea_b - 1
          R_sea_NB = R_sea_NB + 1
        }
      ### In C
      #### At colony
      } else if (transition == "S_c_to_E_c" & S_c > 0) {
        S_c = S_c - 1
        E_c = E_c + 1
      } else if (transition == "E_c_to_S_c" & E_c > 0) {
        E_c = E_c - 1
        S_c = S_c + 1
      } else if (transition == "E_c_to_I_c" & E_c > 0) {
        E_c = E_c - 1
        I_c = I_c + 1
      } else if (transition == "I_c_to_R_c" & I_c > 0) {
        I_c = I_c - 1
        R_c = R_c + 1
      } else if (transition == "I_c_to_D_c" & I_c > 0) {
        # If an adult dies, the partner becomes a non-breeder and the nestling dies
        I_c = I_c - 1
        D_c = D_c + 1
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
        partner = sample(c(rep("S_c", S_c), rep("E_c", E_c),rep("I_c", I_c),rep("R_c", R_c),
                           rep("S_sea_c", S_sea_c), rep("E_sea_c", E_sea_c),rep("I_sea_c", I_sea_c),rep("R_sea_c", R_sea_c)),
                         size = 1)
        
        if (partner == "S_c"){
          S_c = S_c - 1
          S_c_NB = S_c_NB + 1
        } else if (partner == "E_c"){
          E_c = E_c - 1
          E_c_NB = E_c_NB + 1
        } else if (partner == "I_c"){
          I_c = I_c - 1
          I_c_NB = I_c_NB + 1
        } else if (partner == "R_c"){
          R_c = R_c - 1
          R_c_NB = R_c_NB + 1
        } else if (partner == "S_sea_c"){
          S_sea_c = S_sea_c - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_c"){
          E_sea_c = E_sea_c - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_c"){
          I_sea_c = I_sea_c - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_c"){
          R_sea_c = R_sea_c - 1
          R_sea_NB = R_sea_NB + 1
        }
        #### At sea
      } else if (transition == "S_sea_c_to_E_sea_c" & S_sea_c > 0) {
        S_sea_c = S_sea_c - 1
        E_sea_c = E_sea_c + 1
      } else if (transition == "E_sea_c_to_S_sea_c" & E_sea_c > 0) {
        E_sea_c = E_sea_c - 1
        S_sea_c = S_sea_c + 1
      } else if (transition == "E_sea_c_to_I_sea_c" & E_sea_c > 0) {
        E_sea_c = E_sea_c - 1
        I_sea_c = I_sea_c + 1
      } else if (transition == "I_sea_c_to_R_sea_c" & I_sea_c > 0) {
        I_sea_c = I_sea_c - 1
        R_sea_c = R_sea_c + 1
      } else if (transition == "I_sea_c_to_D_sea_c" & I_sea_c > 0) {
        I_sea_c = I_sea_c - 1
        D_sea_c = D_sea_c + 1
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
        partner = sample(c(rep("S_c", S_c), rep("E_c", E_c),rep("I_c", I_c),rep("R_c", R_c),
                           rep("S_sea_c", S_sea_c), rep("E_sea_c", E_sea_c),rep("I_sea_c", I_sea_c),rep("R_sea_c", R_sea_c)),
                         size = 1)
        
        if (partner == "S_c"){
          S_c = S_c - 1
          S_c_NB = S_c_NB + 1
        } else if (partner == "E_c"){
          E_c = E_c - 1
          E_c_NB = E_c_NB + 1
        } else if (partner == "I_c"){
          I_c = I_c - 1
          I_c_NB = I_c_NB + 1
        } else if (partner == "R_c"){
          R_c = R_c - 1
          R_c_NB = R_c_NB + 1
        } else if (partner == "S_sea_c"){
          S_sea_c = S_sea_c - 1
          S_sea_NB = S_sea_NB + 1
        } else if (partner == "E_sea_c"){
          E_sea_c = E_sea_c - 1
          E_sea_NB = E_sea_NB + 1
        } else if (partner == "I_sea_c"){
          I_sea_c = I_sea_c - 1
          I_sea_NB = I_sea_NB + 1
        } else if (partner == "R_sea_c"){
          R_sea_c = R_sea_c - 1
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
      } else if (transition == "S_a_to_S_sea_a"  & S_a > 0) {
        S_a = S_a - 1
        S_sea_a = S_sea_a + 1
      } else if (transition == "E_a_to_E_sea_a"  & E_a > 0) {
        E_a = E_a - 1
        E_sea_a = E_sea_a + 1
      } else if (transition == "I_a_to_I_sea_a"  & I_a > 0) {
        I_a = I_a - 1
        I_sea_a = I_sea_a + 1
      } else if (transition == "R_a_to_R_sea_a"  & R_a > 0) {
        R_a = R_a - 1
        R_sea_a = R_sea_a + 1
      ### From sea A to colony A
      }else if (transition == "S_sea_a_to_S_a"  & S_sea_a > 0) {
        S_a = S_a + 1
        S_sea_a = S_sea_a - 1
      } else if (transition == "E_sea_a_to_E_a"  & E_sea_a > 0) {
        E_a = E_a + 1
        E_sea_a = E_sea_a - 1
      } else if (transition == "I_sea_a_to_I_a"  & I_sea_a > 0) {
        I_a = I_a + 1
        I_sea_a = I_sea_a - 1
      } else if (transition == "R_sea_a_to_R_a" & R_sea_a > 0) {
        R_a = R_a + 1
        R_sea_a = R_sea_a - 1
       ### From colony B to sea B
      } else if (transition == "S_b_to_S_sea_b"  & S_b > 0) {
        S_b = S_b - 1
        S_sea_b = S_sea_b + 1
      } else if (transition == "E_b_to_E_sea_b"  & E_b > 0) {
        E_b = E_b - 1
        E_sea_b = E_sea_b + 1
      } else if (transition == "I_b_to_I_sea_b"  & I_b > 0) {
        I_b = I_b - 1
        I_sea_b = I_sea_b + 1
      } else if (transition == "R_b_to_R_sea_b"  & R_b > 0) {
        R_b = R_b - 1
        R_sea_b = R_sea_b + 1
      ### From sea B to colony B
      } else if (transition == "S_sea_b_to_S_b" & S_sea_b > 0) {
        S_b = S_b + 1
        S_sea_b = S_sea_b - 1
      } else if (transition == "E_sea_b_to_E_b" & E_sea_b > 0) {
        E_b = E_b + 1
        E_sea_b = E_sea_b - 1
      } else if (transition == "I_sea_b_to_I_b" & I_sea_b > 0) {
        I_b = I_b + 1
        I_sea_b = I_sea_b - 1
      } else if (transition == "R_sea_b_to_R_b" & R_sea_b > 0) {
        R_b = R_b + 1
        R_sea_b = R_sea_b - 1 
      ### From colony C to sea C
      } else if (transition == "S_c_to_S_sea_c"  & S_c > 0) {
        S_c = S_c - 1
        S_sea_c = S_sea_c + 1
      } else if (transition == "E_c_to_E_sea_c"  & E_c > 0) {
        E_c = E_c - 1
        E_sea_c = E_sea_c + 1
      } else if (transition == "I_c_to_I_sea_c"  & I_c > 0) {
        I_c = I_c - 1
        I_sea_c = I_sea_c + 1
      } else if (transition == "R_c_to_R_sea_c"  & R_c > 0) {
        R_c = R_c - 1
        R_sea_c = R_sea_c + 1
        ### From sea C to colony C
      } else if (transition == "S_sea_c_to_S_c" & S_sea_c > 0) {
        S_c = S_c + 1
        S_sea_c = S_sea_c - 1
      } else if (transition == "E_sea_c_to_E_c" & E_sea_c > 0) {
        E_c = E_c + 1
        E_sea_c = E_sea_c - 1
      } else if (transition == "I_sea_c_to_I_c" & I_sea_c > 0) {
        I_c = I_c + 1
        I_sea_c = I_sea_c - 1
      } else if (transition == "R_sea_c_to_R_c" & R_sea_c > 0) {
        R_c = R_c + 1
        R_sea_c = R_sea_c - 1 
        
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
      # } else if (transition == "S_a_to_S_a_NB"  & S_a > 0) {
      #   S_a = S_a - 2
      #   S_a_NB = S_a_NB + 2
      # } else if (transition == "E_a_to_E_a_NB" & E_a > 0) {
      #   E_a = E_a - 2
      #   E_a_NB = E_a_NB + 2
      # } else if (transition == "I_a_to_I_a_NB" & I_a > 0) {
      #   I_a = I_a - 2
      #   I_a_NB = I_a_NB + 2
      # } else if (transition == "E_a_to_E_a_NB" & R_a > 0) {
      #   R_a = R_a - 2
      #   R_a_NB = R_a_NB + 2
        
        } # transition
      } # for : transitions_bank
    } # if : transitions_bank

    
    new_state = matrix(data = c(S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                S_a, E_a, I_a, R_a, D_a,
                                S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                
                                S_b_N, E_b_N, I_b_N, R_b_N, D_b_N,
                                S_b, E_b, I_b, R_b, D_b,
                                S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                
                                S_c_N, E_c_N, I_c_N, R_c_N, D_c_N,
                                S_c, E_c, I_c, R_c, D_c,
                                S_sea_c, E_sea_c, I_sea_c, R_sea_c, D_sea_c,
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
    
    S_a = states[2, 1, ],
    E_a = states[2, 2, ],
    I_a = states[2, 3, ],
    R_a = states[2, 4, ],
    D_a = states[2, 5, ],
    
    S_sea_a = states[3, 1, ],
    E_sea_a = states[3, 2, ],
    I_sea_a = states[3, 3, ],
    R_sea_a = states[3, 4, ],
    D_sea_a = states[3, 5, ],
    
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
    
    S_b = states[6, 1, ],
    E_b = states[6, 2, ],
    I_b = states[6, 3, ],
    R_b = states[6, 4, ],
    D_b = states[6, 5, ],
    
    S_sea_b = states[7, 1, ],
    E_sea_b = states[7, 2, ],
    I_sea_b = states[7, 3, ],
    R_sea_b = states[7, 4, ],
    D_sea_b = states[7, 5, ],
  
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
    
    S_c = states[10, 1, ],
    E_c = states[10, 2, ],
    I_c = states[10, 3, ],
    R_c = states[10, 4, ],
    D_c = states[10, 5, ],
    
    S_sea_c = states[11, 1, ],
    E_sea_c = states[11, 2, ],
    I_sea_c = states[11, 3, ],
    R_sea_c = states[11, 4, ],
    D_sea_c = states[11, 5, ],
    
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
    S_a_total = S_a + S_sea_a,
    E_a_total = E_a + E_sea_a,
    I_a_total = I_a + I_sea_a,
    R_a_total = R_a + R_sea_a,
    D_a_total = D_a + D_sea_a,
    
    S_b_total = S_b + S_sea_b,
    E_b_total = E_b + E_sea_b,
    I_b_total = I_b + I_sea_b,
    R_b_total = R_b + R_sea_b,
    D_b_total = D_b + D_sea_b,
    
    S_c_total = S_c + S_sea_c,
    E_c_total = E_c + E_sea_c,
    I_c_total = I_c + I_sea_c,
    R_c_total = R_c + R_sea_c,
    D_c_total = D_c + D_sea_c
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
  output_a = output_long %>% filter(variable %in% c("S_a_total", "E_a_total", "I_a_total", "R_a_total", "D_a_total"))
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
  output_b = output_long %>% filter(variable %in% c("S_b_total", "E_b_total", "I_b_total", "R_b_total", "D_b_total"))
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
  
  output_c = output_long %>% filter(variable %in% c("S_c_total", "E_c_total", "I_c_total", "R_c_total", "D_c_total"))
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

# summary_output ----------------------------------------------------------

summary_output = function(# Time series of all states
                          output, 
                          # Probability of a nestling becoming a breeder
                          reaching.repro.prob = 0.3){
  
  output = output[[1]]

  N_a = output[1, c("S_a", "I_a", "S_sea_a", "I_sea_a")] %>% sum()
  dead_a = output[nrow(output), c("D_a","D_sea_a")] %>% sum()
  a_N = output[nrow(output), c("S_a_N", "E_a_N", "I_a_N", "R_a_N")] %>% sum()
  max_infected_a =  output[, c("E_a", "E_sea_a", "I_a","I_sea_a")] %>% 
    rowSums() %>% 
    max()
  
 
  N_b = output[1, c("S_b", "I_b", "S_sea_b", "I_sea_b")] %>% sum()
  dead_b = output[nrow(output), c("D_b","D_sea_b")] %>% sum()
  b_N = output[nrow(output), c("S_b_N", "E_b_N", "I_b_N", "R_b_N")] %>% sum()
  max_infected_b =  output[, c("E_b", "E_sea_b", "I_b","I_sea_b")] %>% 
    rowSums() %>% 
    max()
  
  N_c = output[1, c("S_c", "I_c", "S_sea_c", "I_sea_c")] %>% sum()
  dead_c = output[nrow(output), c("D_c","D_sea_c")] %>% sum()
  c_N = output[nrow(output), c("S_c_N", "E_c_N", "I_c_N", "R_c_N")] %>% sum()
  max_infected_c =  output[, c("E_c", "E_sea_c", "I_c","I_sea_c")] %>% 
    rowSums() %>% 
    max()


  nb_adults = N_a + N_b + N_c - dead_a - dead_b - dead_c
  nb_nestlings = a_N + b_N + c_N
  nb_adults_equi = nb_adults + reaching.repro.prob * nb_nestlings
  
  
  nb_infected_colonies = 
    sum(max_infected_a > 0,
        max_infected_b > 0,
        max_infected_c > 0)
  
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
  dispersal_stochactic_ = T,
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
    dispersal_stochactic = dispersal_stochactic_,
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

summary_output(output, reaching.repro.prob = 0.4)


model_wrapper(  # Parameter of the taul-leap agorithm
  tau_ = 0.05,
  # Number of simu_adultlation days
  total_time_ = 50,
  # Do we induce dispersion ?
  induced_dispersal_ = F,
  # Induced dispersion mode (deterministic or stochastic)
  dispersal_stochactic_ = T,
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
  reaching.repro.prob_ = 0.3)


# stat_model --------------------------------------------------------------

nb_iterations = 8
#nb_iterations = 25

stat_model = function(nb_iterations_ = nb_iterations,
                      # Do we induce dispersion ?
                      induced_dispersal_ = T,
                      # Induced dispersion mode (deterministic or stochastic)
                      dispersal_stochactic_ = T,
                      # Reaction time between 1rst death and induced dispersal 
                      dispersal_reaction_time_ = 2,
                      # Initial conditions
                      initial_number_infected_breeders_A_ = 1,
                      initial_number_breeders_A_ = 100,
                      initial_number_breeders_B_ = 100,
                      initial_number_breeders_C_ = 10,
                      # Transmission rate from exposed individuals and from infectious individuals in a colony
                      BETA_ = 0.5,
                      # Time at sea before returning to a colony (non-breeders)
                      TIME_AT_SEA_NB_ = 40,
                      # Number of simu_adultlation days
                      total_time_ = 30,
                      # Parameter of the taul-leap agorithm
                      tau_ = 0.1){

  response_list = data.frame()

  for (i in 1:nb_iterations_){

    output = gillespie_seir(induced_dispersal = induced_dispersal_,
                            dispersal_stochactic = dispersal_stochactic_,
                            dispersal_reaction_time = dispersal_reaction_time_,
                            initial_number_infected_breeders_A = initial_number_infected_breeders_A_,
                            initial_number_breeders_A = initial_number_breeders_A_,
                            initial_number_breeders_B = initial_number_breeders_B_,
                            initial_number_breeders_C = initial_number_breeders_C_,
                            TIME_AT_SEA_NB = TIME_AT_SEA_NB_,
                            BETA = BETA_,
                            total_time = total_time_,
                            tau = tau_)

    response_list = rbind(response_list, summary_output(output))

  }
  return(response_list)
}

# res = stat_model()
# res$nb_infected_colonies



# plot - 5 panels - scenario ----------------------------------------------

scenario_dt = function(beta_context,
                         time_at_sea_NB_context
                         ){
  
  no_stress = 
    stat_model(nb_iterations_ = 1,
               induced_dispersal_ = F,
               initial_number_infected_breeders_A_ = 0,
               tau_ = 0.2,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context)
  
  baseline_outbreak = 
    stat_model(induced_dispersal_ = F,
               initial_number_infected_breeders_A_ = 1,
               tau_ = 0.2,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context)
  
  proactive_strategy = 
    stat_model(nb_iterations_ = 1,
               induced_dispersal_ = T,
               initial_number_infected_breeders_A_ = 0,
               dispersal_stochactic_ = F,
               tau_ = 0.2,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context)
  
  proactive_strategy_toolate = 
    stat_model(induced_dispersal_ = T,
               initial_number_infected_breeders_A_ = 1,
               dispersal_stochactic_ = F,
               tau_ = 0.2,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context)
  
  reactive_strategy = 
    stat_model(induced_dispersal_ = T,
               initial_number_infected_breeders_A_ = 1,
               dispersal_stochactic_ = T,
               dispersal_reaction_time_ = 2,
               tau_ = 0.2,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context)
  
  dt = 
    data.frame(
      scenario = c(
        rep("Healty site", nrow(no_stress)),
        rep("Baseline outbreak",nrow(baseline_outbreak)),
        rep("Proactive strategy",nrow(proactive_strategy)),
        rep("Proactive strategy - Too late",nrow(proactive_strategy_toolate)),
        rep("Reactive strategy",nrow(reactive_strategy))
      ),
      
      equi.lost.survi.ad = c(
        0,
        no_stress$nb_adults_equi - baseline_outbreak$nb_adults_equi,
        no_stress$nb_adults_equi - proactive_strategy$nb_adults_equi,
        no_stress$nb_adults_equi - proactive_strategy_toolate$nb_adults_equi,
        no_stress$nb_adults_equi - reactive_strategy$nb_adults_equi 
      ),
      
      nb_infected_colonies = c(
        no_stress$nb_infected_colonies,
        baseline_outbreak$nb_infected_colonies,
        proactive_strategy$nb_infected_colonies,
        proactive_strategy_toolate$nb_infected_colonies,
        reactive_strategy$nb_infected_colonies
      ),
      
      infected_X_time = c(
        no_stress$infected_X_time,
        baseline_outbreak$infected_X_time,
        proactive_strategy$infected_X_time,
        proactive_strategy_toolate$infected_X_time,
        reactive_strategy$infected_X_time
        
      )
      
      
    ) %>% 
    mutate(scenario = factor(scenario, levels = c("Healty site", 
                                                  "Baseline outbreak",
                                                  "Proactive strategy",
                                                  "Proactive strategy - Too late",
                                                  "Reactive strategy")))
  return(dt)
  
}

scenario_mean = function(dt){
  
  equi.lost.survi.ad_mean = c(
  no_stress = mean(dt[dt$scenario == "Healty site", c("equi.lost.survi.ad")]),
  baseline_outbreak = mean(dt[dt$scenario == "Baseline outbreak", c("equi.lost.survi.ad")]),
  proactive_strategy = mean(dt[dt$scenario == "Proactive strategy", c("equi.lost.survi.ad")]),
  proactive_strategy_toolate = mean(dt[dt$scenario == "Proactive strategy - Too late", c("equi.lost.survi.ad")]),
  reactive_strategy = mean(dt[dt$scenario == "Reactive strategy", c("equi.lost.survi.ad")])
  )
  infected_X_time_mean = c(
    no_stress = mean(dt[dt$scenario == "Healty site", c("infected_X_time")]),
    baseline_outbreak = mean(dt[dt$scenario == "Baseline outbreak", c("infected_X_time")]),
    proactive_strategy = mean(dt[dt$scenario == "Proactive strategy", c("infected_X_time")]),
    proactive_strategy_toolate = mean(dt[dt$scenario == "Proactive strategy - Too late", c("infected_X_time")]),
    reactive_strategy = mean(dt[dt$scenario == "Reactive strategy", c("infected_X_time")])
  )
  nb_infected_colonies_mean = c(
    no_stress = mean(dt[dt$scenario == "Healty site", c("nb_infected_colonies")]),
    baseline_outbreak = mean(dt[dt$scenario == "Baseline outbreak", c("nb_infected_colonies")]),
    proactive_strategy = mean(dt[dt$scenario == "Proactive strategy", c("nb_infected_colonies")]),
    proactive_strategy_toolate = mean(dt[dt$scenario == "Proactive strategy - Too late", c("nb_infected_colonies")]),
    reactive_strategy = mean(dt[dt$scenario == "Reactive strategy", c("nb_infected_colonies")])
  )
  
  dt_mean = data.frame(equi.lost.survi.ad_mean,
                       infected_X_time_mean,
                       nb_infected_colonies_mean)
  return(dt_mean)
}



scenario_plot = function(dt){
  
  dt_panel_5 = dt
  
  sc_mean = scenario_mean(dt)
  
  p_equi.survi.ad = ggplot()+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site")), 
                 aes(x = scenario, y = equi.lost.survi.ad),
                 binaxis='y', stackdir='center',
                 dotsize = 2,
                 fill = "grey",
                 color = "white",
                 stackratio=0.05)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Baseline outbreak")), 
                aes(x = scenario, y = equi.lost.survi.ad),
                fill = "grey",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy")), 
                 aes(x = scenario, y = equi.lost.survi.ad),
                 binaxis='y', stackdir='center',
                 dotsize = 2,
                 fill = if (sc_mean["baseline_outbreak","equi.lost.survi.ad_mean"]
                            >sc_mean["proactive_strategy","equi.lost.survi.ad_mean"]) "lightgreen"
                 else "darksalmon",
                 color = "white",
                 stackratio=0.05)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy - Too late")), 
                aes(x = scenario, y = equi.lost.survi.ad),
                fill = if (sc_mean["baseline_outbreak","equi.lost.survi.ad_mean"]
                           >sc_mean["proactive_strategy_toolate","equi.lost.survi.ad_mean"]) "lightgreen"
                       else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Reactive strategy")), 
                aes(x = scenario, y = equi.lost.survi.ad),
                fill = if (sc_mean["baseline_outbreak","equi.lost.survi.ad_mean"]
                           >sc_mean["reactive_strategy","equi.lost.survi.ad_mean"]) "lightgreen"
                else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data = dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy - Too late")), 
                 aes(x = scenario, y = equi.lost.survi.ad),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 color = "darkgrey", alpha = 0.5,
                 stackratio=0.25)+
    geom_dotplot(data = dt_panel_5 %>% subset(., !(scenario %in% c("Proactive strategy - Too late"))), 
                 aes(x = scenario, y = equi.lost.survi.ad),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 color = "darkgrey", alpha = 0.5,
                 stackratio=0.45)+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site", "Proactive strategy")), 
                 aes(x = scenario, y = equi.lost.survi.ad),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 fill = "antiquewhite4",
                 color = "antiquewhite4",
                 stackratio=0.05)+
    ggthemes::theme_clean() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_blank(), # Enlever la bordure du panel
      axis.title = element_text(size = 11),  # Thicken axis titles
      axis.text = element_text(size = 10),  # Thicken axis text
      axis.line = element_line(linewidth = 2),  # Thicken axis lines
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.position =  "none"
    )+
    ylim(0, NA)+
    labs(x = "Scenario", y = "ENLA") 

  p_infected_X_time = ggplot()+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site")), 
                 aes(x = scenario, y = infected_X_time),
                 binaxis='y', stackdir='center',
                 dotsize = 2,
                 fill = "grey",
                 color = "white",
                 stackratio=0.05)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Baseline outbreak")), 
                aes(x = scenario, y = infected_X_time),
                fill = "grey",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy")), 
                 aes(x = scenario, y = infected_X_time),
                 binaxis='y', stackdir='center',
                 dotsize = 2,
                 fill = if (sc_mean["baseline_outbreak","infected_X_time_mean"]
                            >sc_mean["proactive_strategy","infected_X_time_mean"]) "lightgreen"
                 else "darksalmon",
                 color = "white",
                 stackratio=0.05)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy - Too late")), 
                aes(x = scenario, y = infected_X_time),
                fill = if (sc_mean["baseline_outbreak","infected_X_time_mean"]
                           >sc_mean["proactive_strategy_toolate","infected_X_time_mean"]) "lightgreen"
                else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Reactive strategy")), 
                aes(x = scenario, y = infected_X_time),
                fill = if (sc_mean["baseline_outbreak","infected_X_time_mean"]
                           >sc_mean["reactive_strategy","infected_X_time_mean"]) "lightgreen"
                else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data = dt_panel_5 %>% subset(., !(scenario %in% c("Proactive strategy - Too late"))),
                 aes(x = scenario, y = infected_X_time),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 color = "darkgrey", alpha = 0.5,
                 stackratio=0.5)+
    geom_dotplot(data = dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy - Too late")),
                 aes(x = scenario, y = infected_X_time),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 color = "darkgrey", alpha = 0.55,
                 stackratio=0.15)+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site", "Proactive strategy")), 
                 aes(x = scenario, y = infected_X_time),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 fill = "antiquewhite4",
                 color = "antiquewhite4",
                 stackratio=0.05)+
    ggthemes::theme_clean() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_blank(), # Enlever la bordure du panel
      axis.title = element_text(size = 11),  # Thicken axis titles
      axis.text = element_text(size = 10),  # Thicken axis text
      axis.line = element_line(linewidth = 2),  # Thicken axis lines
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.position =  "none"
    )+
    ylim(0, NA)+
    labs(x = "Scenario", y = "Infected x Time") 
  
  p_nb_infected_colonies = ggplot()+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site")), 
                 aes(x = scenario, y = nb_infected_colonies),
                 binaxis='y', stackdir='center',
                 dotsize = 2,
                 fill = "grey",
                 color = "white",
                 stackratio=0.05)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Baseline outbreak")), 
                aes(x = scenario, y = nb_infected_colonies),
                fill = "grey",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy")), 
                 aes(x = scenario, y = nb_infected_colonies),
                 binaxis='y', stackdir='center',
                 dotsize = 2,
                 fill = if (sc_mean["baseline_outbreak","nb_infected_colonies_mean"]
                            >sc_mean["proactive_strategy","nb_infected_colonies_mean"]) "lightgreen"
                 else "darksalmon",
                 color = "white",
                 stackratio=0.05)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy - Too late")), 
                aes(x = scenario, y = nb_infected_colonies),
                fill = if (sc_mean["baseline_outbreak","nb_infected_colonies_mean"]
                           >sc_mean["proactive_strategy_toolate","nb_infected_colonies_mean"]) "lightgreen"
                else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Reactive strategy")), 
                aes(x = scenario, y = nb_infected_colonies),
                fill = if (sc_mean["baseline_outbreak","nb_infected_colonies_mean"]
                           >sc_mean["reactive_strategy","nb_infected_colonies_mean"]) "lightgreen"
                else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data = dt_panel_5, aes(x = scenario, y = nb_infected_colonies),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 color = "darkgrey", alpha = 0.5,
                 stackratio=0.25)+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site", "Proactive strategy")), 
                 aes(x = scenario, y = nb_infected_colonies),
                 binaxis='y', stackdir='center',
                 dotsize = 0.5,
                 fill = "antiquewhite4",
                 color = "antiquewhite4",
                 stackratio=0.05)+
    ggthemes::theme_clean() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_blank(), # Enlever la bordure du panel
      axis.title = element_text(size = 11),  # Thicken axis titles
      axis.text = element_text(size = 10),  # Thicken axis text
      axis.line = element_line(linewidth = 2),  # Thicken axis lines
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.position =  "none"
    )+
    ylim(0, NA)+
    labs(x = "Scenario", y = "Infected colonies") 
  
  
  p = plot_grid(p_equi.survi.ad,
                p_infected_X_time,
                p_nb_infected_colonies,
            labels = c("A", "B", "C"),
            ncol = 3)

  
  print(p)
  
}

dt = scenario_dt(beta_context = 0.5,
                 time_at_sea_NB_context = 40)

scenario_plot(dt)

dt = scenario_dt(beta_context = 0.03,
                 time_at_sea_NB_context = 40)

scenario_plot(dt)

