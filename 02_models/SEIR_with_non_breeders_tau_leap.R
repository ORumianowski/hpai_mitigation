








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
# the localisation : colony A (A), colony (B), sea A (sea_a), sea B (sea_B)
# 
# Thus, an exposed individual, non-breeder and at sea B, is noted E_sea_b_NB


# Parameters --------------------------------------------------------------

param = list( 
  
  # Epidemiological parameters
  
  # Transmission rate from exposed individuals and from infectious individuals in a colony
  beta = matrix(c(0.10, 0.15,
                  0.00, 0.00),
                nrow = 2, ncol = 2, byrow = T),
  # Rate of progression from exposed to infectious (inverse of incubation period)
  sigma = 1/10,
  # Rate of progression from infectious to exposed
  eta =  0.9, 
  # Recovery rate (inverse of infectious period)
  gamma = 1/3,
  # Disease-related mortality rate
  mu = 1/3 * 0.5,
  
  
  # Mobility  parameters
  
  # Transition from sea to the colony (breeders)
  zeta_to_colony = 1,
  # Transition from colony to the sea (breeders)
  zeta_to_sea = 1,
  # Transition from breeder to non-breeder (reproductive failure)
  psi = 0, #Attention départ en couple
  # Transition from sea to the colony (non-breeders)
  rho_to_colony = 1/7,
  # Transition from colony to the sea (non-breeders)
  rho_to_sea = 1/7 * 7,

  
  # Induced dispersion parameters
  
  # Proportion of dispersed adults
  prop_dispersal = 1,
  # Proportion of prospectors among dispersed adults
  prop_prospecting = 1/8,
  # Date of induced dispersion
  dispersal_date = 30,
  # Reaction time between 1rst death and induced dispersal
  dispersal_reaction_time = 5,
  
  # Demographic parameters
  hatching_date = 10

)




# Event rates -------------------------------------------------------------


calculate_rates = function(  beta_E_colony, beta_I_colony,
                             sigma,eta, gamma, mu,
                             zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                             prop_dispersal, prop_prospecting, dispersal_date,
                             hatching_date,
                             S_a, E_a, I_a, R_a, D_a,
                             S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                             S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                             S_b, E_b, I_b, R_b, D_b,
                             S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                             S_sea_a_NB, E_sea_a_NB, I_sea_a_NB, R_sea_a_NB, D_sea_a_NB,
                             S_sea_b_NB,  E_sea_b_NB, I_sea_b_NB, R_sea_b_NB, D_sea_b_NB,
                             S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                             S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                             S_b_N, E_b_N, I_b_N, R_b_N, D_b_N){
  
  rates = c(
    
  # SEIR
  
  ## Nestlings
  ### Colony A
  "S_a_N_to_E_a_N" = beta_E_colony * S_a_N * (E_a+E_a_NB+E_a_N) + 
                     beta_I_colony * S_a_N * (I_a+I_a_NB+I_a_N),
  "E_a_N_to_S_a_N" = eta * E_a_N,
  "E_a_N_to_I_a_N" = sigma * E_a_N,
  "I_a_N_to_R_a_N" = gamma * I_a_N,
  "I_a_N_to_D_a_N" = mu * I_a_N,
  ### Colony B
  "S_b_N_to_E_b_N" = beta_E_colony * S_b_N * (E_b+E_b_NB+E_b_N) +
                     beta_I_colony * S_b_N * (I_b+I_b_NB+I_b_N),
  "E_b_N_to_S_b_N" = eta * E_b_N,
  "E_b_N_to_I_b_N" = sigma * E_b_N,
  "I_b_N_to_R_b_N" = gamma * I_b_N,
  "I_b_N_to_D_b_N" = mu * I_b_N,
  
  ## Non-Breeders
  ### Colony A
  "S_a_NB_to_E_a_NB" = beta_E_colony * S_a_NB * (E_a+E_a_NB+E_a_N) +
                       beta_I_colony * S_a_NB * (I_a+I_a_NB+I_a_N),
  "E_a_NB_to_S_a_NB" = eta * E_a_NB,
  "E_a_NB_to_I_a_NB" = sigma * E_a_NB,
  "I_a_NB_to_R_a_NB" = gamma * I_a_NB,
  "I_a_NB_to_D_a_NB" = mu * I_a_NB,
  ### Sea A
  "S_sea_a_NB_to_E_sea_a_NB" = 0,
  "E_sea_a_NB_to_S_sea_a_NB" = eta * E_sea_a_NB,
  "E_sea_a_NB_to_I_sea_a_NB" = sigma * E_sea_a_NB,
  "I_sea_a_NB_to_R_sea_a_NB" = gamma * I_sea_a_NB,
  "I_sea_a_NB_to_D_sea_a_NB" = mu * I_sea_a_NB,
  ### Sea B
  "S_sea_b_NB_to_E_sea_b_NB" = 0,
  "E_sea_b_NB_to_S_sea_b_NB" = eta * E_sea_b_NB,
  "E_sea_b_NB_to_I_sea_b_NB" = sigma * E_sea_b_NB,
  "I_sea_b_NB_to_R_sea_b_NB" = gamma * I_sea_b_NB,
  "I_sea_b_NB_to_D_sea_b_NB" = mu * I_sea_b_NB,
  ### Colony B
  "S_b_NB_to_E_b_NB" = beta_E_colony * S_b_NB * (E_b+E_b_NB+E_b_N) + 
                       beta_I_colony * S_b_NB * (I_b+I_b_NB+I_b_N),
  "E_b_NB_to_S_b_NB" = eta * E_b_NB,
  "E_b_NB_to_I_b_NB" = sigma * E_b_NB,
  "I_b_NB_to_R_b_NB" = gamma * I_b_NB,
  "I_b_NB_to_D_b_NB" = mu * I_b_NB,
  
  ## Breeders
  ### Colony A
  "S_a_to_E_a" = beta_E_colony * S_a * (E_a+E_a_NB+E_a_N) + 
                 beta_I_colony * S_a * (I_a+I_a_NB+I_a_N),
  "E_a_to_S_a" = eta * E_a,
  "E_a_to_I_a" = sigma * E_a,
  "I_a_to_R_a" = gamma * I_a,
  "I_a_to_D_a" = mu * I_a,
  ### Sea A
  "S_sea_a_to_E_sea_a" = 0,
  "E_sea_a_to_S_sea_a" = eta * E_sea_a,
  "E_sea_a_to_I_sea_a" = sigma * E_sea_a,
  "I_sea_a_to_R_sea_a" = gamma * I_sea_a,
  "I_sea_a_to_D_sea_a" = mu * I_sea_a,
  ### Sea B
  "S_sea_b_to_E_sea_b" = 0,
  "E_sea_b_to_S_sea_b" = eta * E_sea_b,
  "E_sea_b_to_I_sea_b" = sigma * E_sea_b,
  "I_sea_b_to_R_sea_b" = gamma * I_sea_b,
  "I_sea_b_to_D_sea_b" = mu * I_sea_b,
  ### Colony B
  "S_b_to_E_b" = beta_E_colony * S_b * (E_b+E_b_NB+E_b_N) + 
                 beta_I_colony * S_b * (I_b+I_b_NB+I_b_N),
  "E_b_to_S_b" = eta * E_b,
  "E_b_to_I_b" = sigma * E_b,
  "I_b_to_R_b" = gamma * I_b,
  "I_b_to_D_b" = mu * I_b,
  
  
  # Mobility
  ## Non-Breeders
  ### From colony A to sea A
  "S_a_NB_to_S_sea_a_NB" = rho_to_sea * S_a_NB,
  "E_a_NB_to_E_sea_a_NB" = rho_to_sea * E_a_NB,
  "I_a_NB_to_I_sea_a_NB" = rho_to_sea * I_a_NB,
  "R_a_NB_to_R_sea_a_NB" = rho_to_sea * R_a_NB,
  ### From sea A to colony A
  "S_sea_a_NB_to_S_a_NB" = rho_to_colony * S_sea_a_NB,
  "E_sea_a_NB_to_E_a_NB" = rho_to_colony * E_sea_a_NB,
  "I_sea_a_NB_to_I_a_NB" = rho_to_colony * I_sea_a_NB,
  "R_sea_a_NB_to_R_a_NB" = rho_to_colony * R_sea_a_NB,
  ### From colony B to sea B
  "S_b_NB_to_S_sea_b_NB" = rho_to_sea * S_b_NB,
  "E_b_NB_to_E_sea_b_NB" = rho_to_sea * E_b_NB,
  "I_b_NB_to_I_sea_b_NB" = rho_to_sea * I_b_NB,
  "R_b_NB_to_R_sea_b_NB" = rho_to_sea * R_b_NB,
  ### From sea B to colony B
  "S_sea_b_NB_to_S_b_NB" = rho_to_colony * S_sea_b_NB,
  "E_sea_b_NB_to_E_b_NB" = rho_to_colony * E_sea_b_NB,
  "I_sea_b_NB_to_I_b_NB" = rho_to_colony * I_sea_b_NB,
  "R_sea_b_NB_to_R_b_NB" = rho_to_colony * R_sea_b_NB,
  
  ## Breeders
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
  
  ## Prospecting
  ### From A to B
  "S_sea_a_NB_to_S_sea_b_NB" = prop_prospecting * S_sea_a_NB,
  "E_sea_a_NB_to_E_sea_b_NB" = prop_prospecting * E_sea_a_NB,
  "I_sea_a_NB_to_I_sea_b_NB" = prop_prospecting * I_sea_a_NB,
  "R_sea_a_NB_to_R_sea_b_NB" = prop_prospecting * R_sea_a_NB
  
  # Breeders become Non-Breeders
  
  # "S_a_to_S_a_NB" = psi * S_a,
  # "E_a_to_E_a_NB" = psi * E_a,
  # "I_a_to_I_a_NB" = psi * I_a,
  # "R_a_to_R_a_NB" = psi * R_a
  
  )
  return(rates)
}


# Gillespie SEIR model function -------------------------------------------

gillespie_seir = function(param = param, 
                          induced_dispersal = T,
                          initial_number_breeders_A = 50,
                          initial_number_infected_breeders_A = 1, 
                          initial_number_breeders_B = 50,
                          total_time = 70,
                          dispersal_stochactic = T,
                          tau = 0.25) {
  
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
  ## Nestlings
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
  
  
  ## Breeders
  ## In colony A
  N = initial_number_breeders_A                 
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
  N = 50                  
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
  
  ## At sea B
  N = 50                
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
  
  ## In colony B
  N = initial_number_breeders_B                  
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
  
  ## At sea A
  N = 0                  
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_sea_a_NB = c(S = initial_susceptible,
                             E = initial_exposed,
                             I = initial_infected,
                             R = initial_recovered,
                             D = initial_dead)
  
  ## At sea B
  N = 0                  
  initial_infected = 0
  initial_exposed = 0
  initial_recovered = 0
  initial_susceptible = N - initial_infected - initial_exposed - initial_recovered
  initial_dead = 0
  
  initial_state_sea_b_NB = c(S = initial_susceptible,
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
  
  
  initial_state = matrix(data = c(initial_state_A,
                                  initial_state_sea_a,
                                  initial_state_sea_b,
                                  initial_state_B,
                                  initial_state_A_NB,
                                  initial_state_sea_a_NB,
                                  initial_state_sea_b_NB,
                                  initial_state_B_NB,
                                  initial_state_A_N,
                                  initial_state_B_N), 
                         nrow = 10, ncol = 5, 
                         byrow = T)
  
  # Parameters
  
  beta_E_colony = param$beta[1,1]
  beta_I_colony = param$beta[1,2]
  
  sigma = param$sigma
  eta = param$eta
  gamma = param$gamma
  mu = param$mu
  
  zeta_to_colony = param$zeta_to_colony
  zeta_to_sea = param$zeta_to_sea
  psi = param$psi
  rho_to_colony = param$rho_to_colony
  rho_to_sea = param$rho_to_sea
  
  prop_dispersal = param$prop_dispersal
  prop_prospecting = param$prop_prospecting
  dispersal_date = param$dispersal_date
  dispersal_reaction_time = param$dispersal_reaction_time
  
  hatching_date = param$hatching_date
  
  # Initialization
  times = c(0)
  states = array(dim = c(10,5,1), data = initial_state)
  already_dispersed = F
  already_hatched = F
  first_death = F
  first_death_date = NA
  
  # Next event
  while (times[length(times)] < total_time) {
    
    S_a = states[1, 1, dim(states)[3]]
    E_a = states[1, 2, dim(states)[3]]
    I_a = states[1, 3, dim(states)[3]]
    R_a = states[1, 4, dim(states)[3]]
    D_a = states[1, 5, dim(states)[3]]
    
    S_sea_a = states[2, 1, dim(states)[3]]
    E_sea_a = states[2, 2, dim(states)[3]]
    I_sea_a = states[2, 3, dim(states)[3]]
    R_sea_a = states[2, 4, dim(states)[3]]
    D_sea_a = states[2, 5, dim(states)[3]]
    
    S_sea_b = states[3, 1, dim(states)[3]]
    E_sea_b = states[3, 2, dim(states)[3]]
    I_sea_b = states[3, 3, dim(states)[3]]
    R_sea_b = states[3, 4, dim(states)[3]]
    D_sea_b = states[3, 5, dim(states)[3]]
    
    S_b = states[4, 1, dim(states)[3]]
    E_b = states[4, 2, dim(states)[3]]
    I_b = states[4, 3, dim(states)[3]]
    R_b = states[4, 4, dim(states)[3]]
    D_b = states[4, 5, dim(states)[3]]
    
    S_a_NB = states[5, 1, dim(states)[3]]
    E_a_NB = states[5, 2, dim(states)[3]]
    I_a_NB = states[5, 3, dim(states)[3]]
    R_a_NB = states[5, 4, dim(states)[3]]
    D_a_NB = states[5, 5, dim(states)[3]]
    
    S_sea_a_NB = states[6, 1, dim(states)[3]]
    E_sea_a_NB = states[6, 2, dim(states)[3]]
    I_sea_a_NB = states[6, 3, dim(states)[3]]
    R_sea_a_NB = states[6, 4, dim(states)[3]]
    D_sea_a_NB = states[6, 5, dim(states)[3]]

    S_sea_b_NB = states[7, 1, dim(states)[3]]
    E_sea_b_NB = states[7, 2, dim(states)[3]]
    I_sea_b_NB = states[7, 3, dim(states)[3]]
    R_sea_b_NB = states[7, 4, dim(states)[3]]
    D_sea_b_NB = states[7, 5, dim(states)[3]]
    
    S_b_NB = states[8, 1, dim(states)[3]]
    E_b_NB = states[8, 2, dim(states)[3]]
    I_b_NB = states[8, 3, dim(states)[3]]
    R_b_NB = states[8, 4, dim(states)[3]]
    D_b_NB = states[8, 5, dim(states)[3]]
    
    S_a_N = states[9, 1, dim(states)[3]]
    E_a_N = states[9, 2, dim(states)[3]]
    I_a_N = states[9, 3, dim(states)[3]]
    R_a_N = states[9, 4, dim(states)[3]]
    D_a_N = states[9, 5, dim(states)[3]]
    
    S_b_N = states[10, 1, dim(states)[3]]
    E_b_N = states[10, 2, dim(states)[3]]
    I_b_N = states[10, 3, dim(states)[3]]
    R_b_N = states[10, 4, dim(states)[3]]
    D_b_N = states[10, 5, dim(states)[3]]
    
    # Rates of each possible event
    rates = calculate_rates(beta_E_colony, beta_I_colony,
                            sigma,eta, gamma, mu,
                            zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                            prop_dispersal, prop_prospecting, dispersal_date,
                            hatching_date,
                            S_a, E_a, I_a, R_a, D_a,
                            S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                            S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                            S_b, E_b, I_b, R_b, D_b,
                            S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                            S_sea_a_NB, E_sea_a_NB, I_sea_a_NB, R_sea_a_NB, D_sea_a_NB,
                            S_sea_b_NB,  E_sea_b_NB, I_sea_b_NB, R_sea_b_NB, D_sea_b_NB,
                            S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                            S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                            S_b_N, E_b_N, I_b_N, R_b_N, D_b_N)

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
      
      already_hatched = T
      
      new_state = matrix(data = c(S_a, E_a, I_a, R_a, D_a,
                                  S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                  S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                  S_b, E_b, I_b, R_b, D_b,
                                  
                                  S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                  S_sea_a_NB, E_sea_a_NB, I_sea_a_NB, R_sea_a_NB, D_sea_a_NB,
                                  S_sea_b_NB, E_sea_b_NB, I_sea_b_NB, R_sea_b_NB, D_sea_b_NB,
                                  S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                  
                                  S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                  S_b_N, E_b_N, I_b_N, R_b_N, D_b_N),
                         nrow = 10, ncol = 5, 
                         byrow = T)
      
      states = abind(states, new_state)
      times = c(times, hatching_date)
      
      
      # Rates of each possible event
      rates =     rates = calculate_rates(beta_E_colony, beta_I_colony,
                                          sigma,eta, gamma, mu,
                                          zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                                          prop_dispersal, prop_prospecting, dispersal_date,
                                          hatching_date,
                                          S_a, E_a, I_a, R_a, D_a,
                                          S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                          S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                          S_b, E_b, I_b, R_b, D_b,
                                          S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                          S_sea_a_NB, E_sea_a_NB, I_sea_a_NB, R_sea_a_NB, D_sea_a_NB,
                                          S_sea_b_NB,  E_sea_b_NB, I_sea_b_NB, R_sea_b_NB, D_sea_b_NB,
                                          S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                          S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                          S_b_N, E_b_N, I_b_N, R_b_N, D_b_N)
      
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
      ## Did the dispersal occur?
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
          S_sea_a_NB = S_sea_a_NB + disp_S_a + disp_S_sea_a
          
          # Number of exposed  adults who are dispersed from A
          disp_E_a = disp_a["E_a"]
          disp_E_sea_a = disp_a["E_sea_a"]
          # Update of the number of exposed adults
          E_a = E_a - disp_E_a
          E_sea_a = E_sea_a - disp_E_sea_a
          E_sea_a_NB = E_sea_a_NB + disp_E_a + disp_E_sea_a
          
          # Number of infectious adults who are dispersed from A
          disp_I_a = disp_a["I_a"]
          disp_I_sea_a = disp_a["I_sea_a"]
          # Update of the number of infectious adults
          I_a = I_a - disp_I_a
          I_sea_a = I_sea_a - disp_I_sea_a
          I_sea_a_NB = I_sea_a_NB + disp_I_a + disp_I_sea_a
          
          # Number of recovered adults who are dispersed from A
          disp_R_a = disp_a["R_a"]
          disp_R_sea_a = disp_a["R_sea_a"]
          # Update of the number of recovered adults
          R_a = R_a - disp_R_a
          R_sea_a = R_sea_a - disp_R_sea_a
          R_sea_a_NB = R_sea_a_NB + disp_R_a + disp_R_sea_a
          
          # Death of nestlings
          
          D_a_N = D_a_N + S_a_N + E_a_N + I_a_N + R_a_N
          S_a_N = 0
          E_a_N = 0
          I_a_N = 0
          R_a_N = 0
          
          
          already_dispersed = T
          
          new_state = matrix(data = c(S_a, E_a, I_a, R_a, D_a,
                                      S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                      S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                      S_b, E_b, I_b, R_b, D_b,
                                      
                                      S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                      S_sea_a_NB, E_sea_a_NB, I_sea_a_NB, R_sea_a_NB, D_sea_a_NB,
                                      S_sea_b_NB, E_sea_b_NB, I_sea_b_NB, R_sea_b_NB, D_sea_b_NB,
                                      S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                      
                                      S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                      S_b_N, E_b_N, I_b_N, R_b_N, D_b_N),
                             nrow = 10, ncol = 5, 
                             byrow = T)
          
          states = abind(states, new_state)
          times = c(times, next_time)
          
          
          # Rates of each possible event
          rates =     rates = calculate_rates(beta_E_colony, beta_I_colony,
                                              sigma,eta, gamma, mu,
                                              zeta_to_colony, zeta_to_sea, psi, rho_to_colony, rho_to_sea,
                                              prop_dispersal, prop_prospecting, dispersal_date,
                                              hatching_date,
                                              S_a, E_a, I_a, R_a, D_a,
                                              S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                              S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                              S_b, E_b, I_b, R_b, D_b,
                                              S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                              S_sea_a_NB, E_sea_a_NB, I_sea_a_NB, R_sea_a_NB, D_sea_a_NB,
                                              S_sea_b_NB,  E_sea_b_NB, I_sea_b_NB, R_sea_b_NB, D_sea_b_NB,
                                              S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                              S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                              S_b_N, E_b_N, I_b_N, R_b_N, D_b_N)
          
          total_rate = sum(rates)
          
          if (total_rate == 0) {
            break
          }
          
          time_step = rexp(1, total_rate)
          next_time = times[length(times)] + time_step
          
          
        }
        
      }
    }
    
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
      

      
      if (transition == "S_a_to_E_a" & S_a > 0) {
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
          } else if (partner == "E_a_N"){
            E_a_N = E_a_N - 1
            D_a_N = D_a_N + 1
          } else if (partner == "I_a_N"){
            I_a_N = I_a_N - 1
            D_a_N = D_a_N + 1
          } else if (partner == "R_a_N"){
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
          S_sea_a_NB = S_sea_a_NB + 1
        } else if (partner == "E_sea_a"){
          E_sea_a = E_sea_a - 1
          E_sea_a_NB = E_sea_a_NB + 1
        } else if (partner == "I_sea_a"){
          I_sea_a = I_sea_a - 1
          I_sea_a_NB = I_sea_a_NB + 1
        } else if (partner == "R_sea_a"){
          R_sea_a = R_sea_a - 1
          R_sea_a_NB = R_sea_a_NB + 1
        }
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
          } else if (partner == "E_b_N"){
            E_b_N = E_b_N - 1
            D_b_N = D_b_N + 1
          } else if (partner == "I_b_N"){
            I_b_N = I_b_N - 1
            D_b_N = D_b_N + 1
          } else if (partner == "R_b_N"){
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
          S_sea_b_NB = S_sea_b_NB + 1
        } else if (partner == "E_sea_b"){
          E_sea_b = E_sea_b - 1
          E_sea_b_NB = E_sea_b_NB + 1
        } else if (partner == "I_sea_b"){
          I_sea_b = I_sea_b - 1
          I_sea_b_NB = I_sea_b_NB + 1
        } else if (partner == "R_sea_b"){
          R_sea_b = R_sea_b - 1
          R_sea_b_NB = R_sea_b_NB + 1
        }
        
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
      } else if (transition == "S_sea_a_NB_to_E_sea_a_NB" & S_sea_a_NB > 0) {
        S_sea_a_NB = S_sea_a_NB - 1
        E_sea_a_NB = E_sea_a_NB + 1
      } else if (transition == "E_sea_a_NB_to_S_sea_a_NB" & E_sea_a_NB > 0) {
        E_sea_a_NB = E_sea_a_NB - 1
        S_sea_a_NB = S_sea_a_NB + 1
      } else if (transition == "E_sea_a_NB_to_I_sea_a_NB" & E_sea_a_NB > 0) {
        E_sea_a_NB = E_sea_a_NB - 1
        I_sea_a_NB = I_sea_a_NB + 1
      } else if (transition == "I_sea_a_NB_to_R_sea_a_NB" & I_sea_a_NB > 0) {
        I_sea_a_NB = I_sea_a_NB - 1
        R_sea_a_NB = R_sea_a_NB + 1
      } else if (transition == "I_sea_a_NB_to_D_sea_a_NB" & I_sea_a_NB > 0) {
        I_sea_a_NB = I_sea_a_NB - 1
        D_sea_a_NB = D_sea_a_NB + 1
      } else if (transition == "S_sea_b_NB_to_E_sea_b_NB" & S_sea_b_NB > 0) {
        S_sea_b_NB = S_sea_b_NB - 1
        E_sea_b_NB = E_sea_b_NB + 1
      } else if (transition == "E_sea_b_NB_to_S_sea_b_NB" & E_sea_b_NB > 0) {
        E_sea_b_NB = E_sea_b_NB - 1
        S_sea_b_NB = S_sea_b_NB + 1
      } else if (transition == "E_sea_b_NB_to_I_sea_b_NB" & E_sea_b_NB > 0) {
        E_sea_b_NB = E_sea_b_NB - 1
        I_sea_b_NB = I_sea_b_NB + 1
      } else if (transition == "I_sea_b_NB_to_R_sea_b_NB" & I_sea_b_NB > 0) {
        I_sea_b_NB = I_sea_b_NB - 1
        R_sea_b_NB = R_sea_b_NB + 1
      } else if (transition == "I_sea_b_NB_to_D_sea_b_NB" & I_sea_b_NB > 0) {
        I_sea_b_NB = I_sea_b_NB - 1
        D_sea_b_NB = D_sea_b_NB + 1
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
      } else if (transition == "S_a_NB_to_S_sea_a_NB" & S_a_NB > 0) {
        S_a_NB = S_a_NB - 1
        S_sea_a_NB = S_sea_a_NB + 1
      } else if (transition == "E_a_NB_to_E_sea_a_NB" & E_a_NB > 0) {
        E_a_NB = E_a_NB - 1
        E_sea_a_NB = E_sea_a_NB + 1
      } else if (transition == "I_a_NB_to_I_sea_a_NB" & I_a_NB > 0) {
        I_a_NB = I_a_NB - 1
        I_sea_a_NB = I_sea_a_NB + 1
      } else if (transition == "R_a_NB_to_R_sea_a_NB" & R_a_NB > 0) {
        R_a_NB = R_a_NB - 1
        R_sea_a_NB = R_sea_a_NB + 1
      } else if (transition == "S_b_NB_to_S_sea_b_NB" & S_b_NB > 0) {
        S_b_NB = S_b_NB - 1
        S_sea_b_NB = S_sea_b_NB + 1
      } else if (transition == "E_b_NB_to_E_sea_b_NB" & E_b_NB > 0) {
        E_b_NB = E_b_NB - 1
        E_sea_b_NB = E_sea_b_NB + 1
      } else if (transition == "I_b_NB_to_I_sea_b_NB" & I_b_NB > 0) {
        I_b_NB = I_b_NB - 1
        I_sea_b_NB = I_sea_b_NB + 1
      } else if (transition == "R_b_NB_to_R_sea_b_NB" & R_b_NB > 0) {
        R_b_NB = R_b_NB - 1
        R_sea_b_NB = R_sea_b_NB + 1
      }else if (transition == "S_sea_a_NB_to_S_a_NB" & S_sea_a_NB > 0) {
        S_a_NB = S_a_NB + 1
        S_sea_a_NB = S_sea_a_NB - 1
      } else if (transition == "E_sea_a_NB_to_E_a_NB" & E_sea_a_NB > 0) {
        E_a_NB = E_a_NB + 1
        E_sea_a_NB = E_sea_a_NB - 1
      } else if (transition == "I_sea_a_NB_to_I_a_NB" & I_sea_a_NB > 0) {
        I_a_NB = I_a_NB + 1
        I_sea_a_NB = I_sea_a_NB - 1
      } else if (transition == "R_sea_a_NB_to_R_a_NB" & R_sea_a_NB > 0) {
        R_a_NB = R_a_NB + 1
        R_sea_a_NB = R_sea_a_NB - 1
      } else if (transition == "S_sea_b_NB_to_S_b_NB" & S_sea_b_NB > 0) {
        S_b_NB = S_b_NB + 1
        S_sea_b_NB = S_sea_b_NB - 1
      } else if (transition == "E_sea_b_NB_to_E_b_NB" & E_sea_b_NB > 0) {
        E_b_NB = E_b_NB + 1
        E_sea_b_NB = E_sea_b_NB - 1
      } else if (transition == "I_sea_b_NB_to_I_b_NB" & I_sea_b_NB > 0) {
        I_b_NB = I_b_NB + 1
        I_sea_b_NB = I_sea_b_NB - 1
      } else if (transition == "R_sea_b_NB_to_R_b_NB" & R_sea_b_NB > 0) {
        R_b_NB = R_b_NB + 1
        R_sea_b_NB = R_sea_b_NB - 1
      } else if (transition == "S_a_to_S_a_NB"  & S_a > 0) {
        S_a = S_a - 2
        S_a_NB = S_a_NB + 2
      } else if (transition == "E_a_to_E_a_NB" & E_a > 0) {
        E_a = E_a - 2
        E_a_NB = E_a_NB + 2
      } else if (transition == "I_a_to_I_a_NB" & I_a > 0) {
        I_a = I_a - 2
        I_a_NB = I_a_NB + 2
      } else if (transition == "E_a_to_E_a_NB" & R_a > 0) {
        R_a = R_a - 2
        R_a_NB = R_a_NB + 2
      } else if (transition == "S_sea_a_NB_to_S_sea_b_NB" & S_sea_a_NB > 0) {
        S_sea_a_NB = S_sea_a_NB - 1
        S_sea_b_NB = S_sea_b_NB + 1
      }else if (transition == "E_sea_a_NB_to_E_sea_b_NB" & E_sea_a_NB > 0) {
        E_sea_a_NB = E_sea_a_NB - 1
        E_sea_b_NB = E_sea_b_NB + 1
      }else if (transition == "I_sea_a_NB_to_I_sea_b_NB" & I_sea_a_NB > 0) {
        I_sea_a_NB = I_sea_a_NB - 1
        I_sea_b_NB = I_sea_b_NB + 1
      }else if (transition == "R_sea_a_NB_to_R_sea_b_NB" & R_sea_a_NB > 0) {
        R_sea_a_NB = R_sea_a_NB - 1
        R_sea_b_NB = R_sea_b_NB + 1
      }else if (transition == "S_a_N_to_E_a_N" & S_a_N > 0){
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
          S_sea_a_NB = S_sea_a_NB + 1
        } else if (parent1 == "E_sea_a"){
          E_sea_a = E_sea_a - 1
          E_sea_a_NB = E_sea_a_NB + 1
        } else if (parent1 == "I_sea_a"){
          I_sea_a = I_sea_a - 1
          I_sea_a_NB = I_sea_a_NB + 1
        } else if (parent1 == "R_sea_a"){
          R_sea_a = R_sea_a - 1
          R_sea_a_NB = R_sea_a_NB + 1
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
          S_sea_a_NB = S_sea_a_NB + 1
        } else if (parent2 == "E_sea_a"){
          E_sea_a = E_sea_a - 1
          E_sea_a_NB = E_sea_a_NB + 1
        } else if (parent2 == "I_sea_a"){
          I_sea_a = I_sea_a - 1
          I_sea_a_NB = I_sea_a_NB + 1
        } else if (parent2 == "R_sea_a"){
          R_sea_a = R_sea_a - 1
          R_sea_a_NB = R_sea_a_NB + 1
        }
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
          S_sea_b_NB = S_sea_b_NB + 1
        } else if (parent1 == "E_sea_b"){
          E_sea_b = E_sea_b - 1
          E_sea_b_NB = E_sea_b_NB + 1
        } else if (parent1 == "I_sea_b"){
          I_sea_b = I_sea_b - 1
          I_sea_b_NB = I_sea_b_NB + 1
        } else if (parent1 == "R_sea_b"){
          R_sea_b = R_sea_b - 1
          R_sea_b_NB = R_sea_b_NB + 1
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
          S_sea_b_NB = S_sea_b_NB + 1
        } else if (parent2 == "E_sea_b"){
          E_sea_b = E_sea_b - 1
          E_sea_b_NB = E_sea_b_NB + 1
        } else if (parent2 == "I_sea_b"){
          I_sea_b = I_sea_b - 1
          I_sea_b_NB = I_sea_b_NB + 1
        } else if (parent2 == "R_sea_b"){
          R_sea_b = R_sea_b - 1
          R_sea_b_NB = R_sea_b_NB + 1
        }
      }
    }
  }
    
    new_state = matrix(data = c(S_a, E_a, I_a, R_a, D_a,
                                S_sea_a, E_sea_a, I_sea_a, R_sea_a, D_sea_a,
                                S_sea_b, E_sea_b, I_sea_b, R_sea_b, D_sea_b,
                                S_b, E_b, I_b, R_b, D_b,
                                
                                S_a_NB, E_a_NB, I_a_NB, R_a_NB, D_a_NB,
                                S_sea_a_NB, E_sea_a_NB, I_sea_a_NB, R_sea_a_NB, D_sea_a_NB,
                                S_sea_b_NB, E_sea_b_NB, I_sea_b_NB, R_sea_b_NB, D_sea_b_NB,
                                S_b_NB, E_b_NB, I_b_NB, R_b_NB, D_b_NB,
                                
                                S_a_N, E_a_N, I_a_N, R_a_N, D_a_N,
                                S_b_N, E_b_N, I_b_N, R_b_N, D_b_N),
                           nrow = 10, ncol = 5, 
                           byrow = T)
    
    
    
    states = abind(states, new_state)

  }
  
  output = data.frame(
    time = times,
    S_a = states[1, 1, ],
    E_a = states[1, 2, ],
    I_a = states[1, 3, ],
    R_a = states[1, 4, ],
    D_a = states[1, 5, ],
    S_sea_a = states[2, 1, ],
    E_sea_a = states[2, 2, ],
    I_sea_a = states[2, 3, ],
    R_sea_a = states[2, 4, ],
    D_sea_a = states[2, 5, ],
    S_sea_b = states[3, 1, ],
    E_sea_b = states[3, 2, ],
    I_sea_b = states[3, 3, ],
    R_sea_b = states[3, 4, ],
    D_sea_b = states[3, 5, ],
    S_b = states[4, 1, ],
    E_b = states[4, 2, ],
    I_b = states[4, 3, ],
    R_b = states[4, 4, ],
    D_b = states[4, 5, ],
    
    S_a_NB = states[5, 1, ],
    E_a_NB = states[5, 2, ],
    I_a_NB = states[5, 3, ],
    R_a_NB = states[5, 4, ],
    D_a_NB = states[5, 5, ],
    S_sea_a_NB = states[6, 1, ],
    E_sea_a_NB = states[6, 2, ],
    I_sea_a_NB = states[6, 3, ],
    R_sea_a_NB = states[6, 4, ],
    D_sea_a_NB = states[6, 5, ],
    S_sea_b_NB = states[7, 1, ],
    E_sea_b_NB = states[7, 2, ],
    I_sea_b_NB = states[7, 3, ],
    R_sea_b_NB = states[7, 4, ],
    D_sea_b_NB = states[7, 5, ],
    S_b_NB = states[8, 1, ],
    E_b_NB = states[8, 2, ],
    I_b_NB = states[8, 3, ],
    R_b_NB = states[8, 4, ],
    D_b_NB = states[8, 5, ],
    
    S_a_N = states[9, 1, ],
    E_a_N = states[9, 2, ],
    I_a_N = states[9, 3, ],
    R_a_N = states[9, 4, ],
    D_a_N = states[9, 5, ],
    S_b_N = states[10, 1, ],
    E_b_N = states[10, 2, ],
    I_b_N = states[10, 3, ],
    R_b_N = states[10, 4, ],
    D_b_N = states[10, 5, ]) %>% 
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
    
    S_a_NB_total = S_a_NB + S_sea_a_NB,
    E_a_NB_total = E_a_NB + E_sea_a_NB,
    I_a_NB_total = I_a_NB + I_sea_a_NB,
    R_a_NB_total = R_a_NB + R_sea_a_NB,
    D_a_NB_total = D_a_NB + D_sea_a_NB,
    
    S_b_NB_total = S_b_NB + S_sea_b_NB,
    E_b_NB_total = E_b_NB + E_sea_b_NB,
    I_b_NB_total = I_b_NB + I_sea_b_NB,
    R_b_NB_total = R_b_NB + R_sea_b_NB,
    D_b_NB_total = D_b_NB + D_sea_b_NB
  )
  
  return(output)
}


# Plot results

plot_seir = function(output_ = output){
  
  output_long = melt(output_[1:nrow(output_)-1, ], id = "time")
  
  output_a = output_long %>% filter(variable %in% c("S_a_total", "E_a_total", "I_a_total", "R_a_total", "D_a_total"))
  plot_a = ggplot(output_a, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Breeders in A (colony+Sea)")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))
  
  
  output_a_N = output_long %>% filter(variable %in% c("S_a_N", "E_a_N", "I_a_N", "R_a_N", "D_a_N"))
  plot_a_N = ggplot(output_a_N, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Nestlings in A")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))
  
  
  output_a_NB = output_long %>% filter(variable %in% c("S_a_NB_total", "E_a_NB_total", "I_a_NB_total", "R_a_NB_total", "D_a_NB_total"))
  plot_a_NB = ggplot(output_a_NB, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Non-Breeders in A (colony+Sea)")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))
  
  
  output_b = output_long %>% filter(variable %in% c("S_b_total", "E_b_total", "I_b_total", "R_b_total", "D_b_total"))
  plot_b = ggplot(output_b, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Breeders in B (colony+Sea)")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))
  
  
  output_b_N = output_long %>% filter(variable %in% c("S_b_N", "E_b_N", "I_b_N", "R_b_N", "D_b_N"))
  plot_b_N = ggplot(output_b_N, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Nestlings in B")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))
  
  
  output_b_NB = output_long %>% filter(variable %in% c("S_b_NB_total", "E_b_NB_total", "I_b_NB_total", "R_b_NB_total", "D_b_NB_total"))
  plot_b_NB = ggplot(output_b_NB, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Number of individuals", color = "Status") +
    theme_minimal() +
    ggtitle("Non-Breeders in B (colony+Sea)")+
    scale_color_brewer(palette="Set2", labels = c("S", "E", "I", "R", "D"))
  
  plot_grid_seir = plot_grid(plot_a, plot_a_N, plot_a_NB,
                             plot_b, plot_b_N, plot_b_NB,
                             labels = c("A", "B", "C",
                                        "D", "E", "F"),
                             label_size = 12)
  
  print(plot_grid_seir)
  
}

# Run simulation

time1 <- Sys.time()
output = gillespie_seir(param = param, 
                        induced_dispersal = T,
                        initial_number_breeders_A = 50,
                        initial_number_infected_breeders_A = 1, 
                        initial_number_breeders_B = 50,
                        total_time = 70,
                        dispersal_stochactic = T,
                        tau = 0.2)

time2 <- Sys.time()
time2 - time1
plot_seir()




summary_output = function(output){

  N_a = output[1, c("S_a", "I_a")] %>% sum()

  max_infected_a = max(output[, "I_a"])
  prop_max_infected_a = max_infected_a / N_a
  dead_a = output[nrow(output), "D_a"]
  prop_dead_a = dead_a / N_a

  max_infected_sea = max(output[, "I_sea_a"])
  dead_sea = output[nrow(output), "D_sea_a"]

  prop_non_exposed_from_a = (output[nrow(output), "S_a"] + output[nrow(output), "S_sea"] - output[1, "S_sea"] ) / N_a

  return( data.frame(
    N_a = N_a,
    max_infected_a = max_infected_a,
    prop_max_infected_a = prop_max_infected_a,
    dead_a = dead_a,
    prop_dead_a = prop_dead_a,
    max_infected_sea = max_infected_sea,
    dead_sea = dead_sea,
    prop_non_exposed_from_a = prop_non_exposed_from_a

  ))
}


output_long_list = data.frame()
response_list = data.frame()

nb_iterations = 8

for (i in 1:nb_iterations){

  output = gillespie_seir(param = param,
                          induced_dispersal = T,
                          initial_number_breeders_A = 50,
                          initial_number_infected_breeders_A = 1,
                          initial_number_breeders_B = 50,
                          total_time = 70,
                          dispersal_stochactic = T,
                          tau = 0.2)

  output_long = melt(output, id = "time")

  output_long_i = cbind(output_long,
                         data.frame(simulation = rep(i, times = nrow(output_long))))

  output_long_list = rbind(output_long_list, output_long_i)
  #response_list = rbind(response_list, summary_output(output))

}


output_a = output_long_list %>% filter(variable %in% c("S_a", "E_a", "I_a", "R_a", "D_a"))

p = ggplot()
for(i in 1:nb_iterations){
  p = p + geom_line(data = output_a %>% subset(., simulation == i )
                    , aes(x = time, y = value, color = variable))
}
p = p +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("Stochastic SEIR Model Simulation (Gillespie Algorithm)")


p


data_long = pivot_longer(response_list, cols = -N_a, names_to = "variable", values_to = "value")

# Créer les diagrammes en violon pour chaque variable
ggplot(data_long %>% subset(., variable %in% c("dead_a", "max_infected_a")),
       aes(x = variable, y = value)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)


