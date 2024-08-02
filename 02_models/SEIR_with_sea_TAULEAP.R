library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)

# Parameters
total_time <- 160          # Simulation time

sigma <- 1/5.2             # Rate of progression from exposed to infectious (inverse of incubation period)
gamma <- 1/2.9             # Recovery rate (inverse of infectious period)
eta <-  (1/5.2)*(0.1)*50   # Rate of progression from infectious to exposed
mu <- (1/2.9) * 0.5        # Rate of mortality 

beta_colony <- 0.01        # Transmission rate in a colony
beta_sea <- 0.001          # Transmission rate at sea
beta = c(beta_colony, beta_sea)


tau = 0.000000             # Transition from colony A to the sea

# Initial state

## In colony A
N_A <- 1000                  # Total population in colony A
initial_infected_A <- 1

initial_exposed_A <- 0
initial_recovered_A <- 0
initial_susceptible_A <- N_A - initial_infected_A - initial_exposed_A - initial_recovered_A
initial_dead_A <- 0

initial_state_A <- c(S = initial_susceptible_A,
                     E = initial_exposed_A,
                     I = initial_infected_A,
                     R = initial_recovered_A,
                     D = initial_dead_A)

## At sea
N_sea <- 50                  # Total population at sea
initial_infected_sea <- 0

initial_exposed_sea <- 0
initial_recovered_sea <- 0
initial_susceptible_sea <- N_sea - initial_infected_sea - initial_exposed_sea - initial_recovered_sea
initial_dead_sea <- 0

initial_state_sea <- c(S = initial_susceptible_sea,
                       E = initial_exposed_sea,
                       I = initial_infected_sea,
                       R = initial_recovered_sea,
                       D = initial_dead_sea)

initial_state = matrix(data = c(initial_state_A,
                                initial_state_sea), 
                       nrow = 2, ncol = 5, 
                       byrow = T)

# τ-leap SEIR model function
tau_leap_seir <- function(N, beta, sigma, gamma, initial_state, total_time, tau) {
  times <- seq(0, total_time, by = tau)
  num_steps <- length(times)
  
  states <- array(dim = c(2, 5, num_steps), data = NA)
  states[,,1] <- initial_state
  
  for (t in 2:num_steps) {
    S_a <- states[1, 1, t-1]
    E_a <- states[1, 2, t-1]
    I_a <- states[1, 3, t-1]
    R_a <- states[1, 4, t-1]
    D_a <- states[1, 5, t-1]
    
    S_sea <- states[2, 1, t-1]
    E_sea <- states[2, 2, t-1]
    I_sea <- states[2, 3, t-1]
    R_sea <- states[2, 4, t-1]
    D_sea <- states[2, 5, t-1]
    
    rates <- c(
      "S_a_to_E_a" = beta[1] * S_a * I_a,
      "E_a_to_S_a" = eta * E_a,
      "E_a_to_I_a" = sigma * E_a,
      "I_a_to_R_a" = gamma * I_a,
      "I_a_to_D_a" = mu * I_a,
      
      "S_sea_to_E_sea" = beta[2] * S_sea * I_sea,
      "E_sea_to_S_sea" = eta * E_sea,
      "E_sea_to_I_sea" = sigma * E_sea,
      "I_sea_to_R_sea" = gamma * I_sea,
      "I_sea_to_D_sea" = mu * I_sea,
      
      "S_a_to_S_sea" = tau * S_a,
      "E_a_to_E_sea" = tau * E_a,
      "I_a_to_I_sea" = tau * I_a,
      "R_a_to_R_sea" = tau * R_a
    )
    
    events <- rpois(length(rates), rates * tau)
    
    S_a <- S_a - events[1] + events[2] - events[11]
    E_a <- E_a + events[1] - events[2] - events[3] - events[12]
    I_a <- I_a + events[3] - events[4] - events[5] - events[13]
    R_a <- R_a + events[4] - events[14]
    D_a <- D_a + events[5]
    
    S_sea <- S_sea + events[11] - events[6]
    E_sea <- E_sea + events[12] + events[6] - events[7] - events[8]
    I_sea <- I_sea + events[13] + events[8] - events[9] - events[10]
    R_sea <- R_sea + events[14] + events[9]
    D_sea <- D_sea + events[10]
    
    new_state = matrix(data = c(S_a, E_a, I_a, R_a, D_a,
                                S_sea, E_sea, I_sea, R_sea, D_sea), 
                       nrow = 2, ncol = 5, 
                       byrow = T)
    
    states[,,t] <- new_state
  }
  
  return(list(times = times, states = states))
}

# Run simulation
result <- tau_leap_seir(N, beta, sigma, gamma, initial_state, total_time, tau = 0.1)

# Convert results to data frame for plotting
output <- data.frame(time = result$times,
                     
                     S_a = result$states[1, 1,],
                     E_a = result$states[1, 2,],
                     I_a = result$states[1, 3,],
                     R_a = result$states[1, 4,],
                     D_a = result$states[1, 5,],
                     
                     S_sea = result$states[2, 1,],
                     E_sea = result$states[2, 2,],
                     I_sea = result$states[2, 3,],
                     R_sea = result$states[2, 4,],
                     D_sea = result$states[2, 5,])

output_long <- melt(output, id = "time")

# Plot results

output_a = output_long %>% subset(., variable %in% c("S_a", "E_a", "I_a", "R_a", "D_a"))
plot_a = ggplot(output_a, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("τ-leap SEIR Model Simulation")

output_sea = output_long %>% subset(., variable %in% c("S_sea", "E_sea", "I_sea", "R_sea", "D_sea"))
plot_sea = ggplot(output_sea, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("τ-leap SEIR Model Simulation")

plot_a 
plot_sea
