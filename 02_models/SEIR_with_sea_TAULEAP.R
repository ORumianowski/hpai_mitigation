# -------------------------------------------------------------------------
# transition A to sea
# -------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)

# Parameters
total_time <- 160          # Simulation time

param <- list(
  sigma = 1/5.2,            # Rate of progression from exposed to infectious (inverse of incubation period)
  gamma = 1/2.9,            # Recovery rate (inverse of infectious period)
  eta =  (1/5.2)*(0.1)*50,  # Rate of progression from infectious to exposed
  mu = (1/2.9) * 0.5,       # Rate of mortality 
  beta = c(0.01, 0.001),    # Transmission rate in a colony, at sea
  zeta = 0.001              # Transition from colony A to the sea
)

# Initial state

## In colony A
N_A <- 1000                # Total population in colony A
initial_infected_A <- 1
initial_state_A <- c(S = N_A - initial_infected_A,
                     E = 0,
                     I = initial_infected_A,
                     R = 0,
                     D = 0)

## At sea
N_sea <- 50                # Total population at sea
initial_infected_sea <- 0
initial_state_sea <- c(S = N_sea - initial_infected_sea,
                       E = 0,
                       I = initial_infected_sea,
                       R = 0,
                       D = 0)

initial_state <- matrix(c(initial_state_A,
                          initial_state_sea),
                        nrow = 2, byrow = TRUE)

# τ-leap SEIR model function
tau_leap_seir <- function(param, initial_state, total_time, tau) {
  
  sigma <- param$sigma
  gamma <- param$gamma
  eta <- param$eta
  mu <- param$mu
  beta_colony <- param$beta[1]
  beta_sea <- param$beta[2]
  zeta <- param$zeta
  
  times <- c(0)
  states <- array(dim = c(2, 5, 1), data = initial_state)
  
  while (times[length(times)] < total_time) {
    
    S_a <- states[1, 1, dim(states)[3]]
    E_a <- states[1, 2, dim(states)[3]]
    I_a <- states[1, 3, dim(states)[3]]
    R_a <- states[1, 4, dim(states)[3]]
    D_a <- states[1, 5, dim(states)[3]]
    
    S_sea <- states[2, 1, dim(states)[3]]
    E_sea <- states[2, 2, dim(states)[3]]
    I_sea <- states[2, 3, dim(states)[3]]
    R_sea <- states[2, 4, dim(states)[3]]
    D_sea <- states[2, 5, dim(states)[3]]
    
    rates <- c(
      "S_a_to_E_a" = beta_colony * S_a * I_a,
      "E_a_to_S_a" = eta * E_a,
      "E_a_to_I_a" = sigma * E_a,
      "I_a_to_R_a" = gamma * I_a,
      "I_a_to_D_a" = mu * I_a,
      
      "S_sea_to_E_sea" = beta_sea * S_sea * I_sea,
      "E_sea_to_S_sea" = eta * E_sea,
      "E_sea_to_I_sea" = sigma * E_sea,
      "I_sea_to_R_sea" = gamma * I_sea,
      "I_sea_to_D_sea" = mu * I_sea,
      
      "S_a_to_S_sea" = zeta * S_a,
      "E_a_to_E_sea" = zeta * E_a,
      "I_a_to_I_sea" = zeta * I_a,
      "R_a_to_R_sea" = zeta * R_a
    )
    
    total_rate <- sum(rates)
    
    if (total_rate == 0) {
      break
    }
    
    time_step <- min(tau, (total_time - times[length(times)]))
    times <- c(times, times[length(times)] + time_step)
    
    transitions <- rpois(length(rates), rates * time_step) # juste ??
    names(transitions) <- names(rates)
    

    S_a <- S_a - transitions["S_a_to_E_a"] + transitions["E_a_to_S_a"] - transitions["S_a_to_S_sea"]
    E_a <- E_a + transitions["S_a_to_E_a"] - transitions["E_a_to_S_a"] - transitions["E_a_to_I_a"] - transitions["E_a_to_E_sea"]
    I_a <- I_a + transitions["E_a_to_I_a"] - transitions["I_a_to_R_a"] - transitions["I_a_to_D_a"] - transitions["I_a_to_I_sea"]
    R_a <- R_a + transitions["I_a_to_R_a"] - transitions["R_a_to_R_sea"]
    D_a <- D_a + transitions["I_a_to_D_a"]
    
    S_sea <- S_sea - transitions["S_sea_to_E_sea"] + transitions["E_sea_to_S_sea"] + transitions["S_a_to_S_sea"]
    E_sea <- E_sea + transitions["S_sea_to_E_sea"] - transitions["E_sea_to_S_sea"] - transitions["E_sea_to_I_sea"] + transitions["E_a_to_E_sea"]
    I_sea <- I_sea + transitions["E_sea_to_I_sea"] - transitions["I_sea_to_R_sea"] - transitions["I_sea_to_D_sea"] + transitions["I_a_to_I_sea"]
    R_sea <- R_sea + transitions["I_sea_to_R_sea"] + transitions["R_a_to_R_sea"]
    D_sea <- D_sea + transitions["I_sea_to_D_sea"]
    
    new_state <- matrix(c(S_a, E_a, I_a, R_a, D_a,
                          S_sea, E_sea, I_sea, R_sea, D_sea),
                        nrow = 2, byrow = TRUE)
    
    states <- abind(states, new_state)
  }
  
  return(list(times = times, states = states))
}

# Run simulation
tau <- 0.1 # Set the tau value for τ-leap
result <- tau_leap_seir(param, initial_state, total_time, tau)

# Convert results to data frame for plotting
output <- data.frame(
  time = result$times,
  S_a = result$states[1, 1, ],
  E_a = result$states[1, 2, ],
  I_a = result$states[1, 3, ],
  R_a = result$states[1, 4, ],
  D_a = result$states[1, 5, ],
  S_sea = result$states[2, 1, ],
  E_sea = result$states[2, 2, ],
  I_sea = result$states[2, 3, ],
  R_sea = result$states[2, 4, ],
  D_sea = result$states[2, 5, ]
)

output_long <- melt(output, id = "time")

# Plot results

output_a <- output_long %>% filter(variable %in% c("S_a", "E_a", "I_a", "R_a", "D_a"))
plot_a <- ggplot(output_a, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("Stochastic SEIR Model Simulation (τ-leap Algorithm)")

output_sea <- output_long %>% filter(variable %in% c("S_sea", "E_sea", "I_sea", "R_sea", "D_sea"))
plot_sea <- ggplot(output_sea, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("Stochastic SEIR Model Simulation (τ-leap Algorithm)")

# Display plots
plot_a 
plot_sea
