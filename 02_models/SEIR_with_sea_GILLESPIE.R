


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------



library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)

# Parameters
total_time <- 50          # Simulation time

param = list( sigma = 1/5.2,            # Rate of progression from exposed to infectious (inverse of incubation period)
              gamma = 1/2.9,             # Recovery rate (inverse of infectious period)
              eta =  (1/5.2)*(0.1)*50,   # Rate of progression from infectious to exposed
              mu = (1/2.9) * 0.5,        # Rate of mortality 
              beta = matrix(c(0.02, 0.01, 0.03, 0,
                               0.001, 0.001, 0.001, 0.001),
                            nrow = 2, ncol = 4),     # Transmission rate in a colony, at sea
              zeta = 0.001                # Transition from colony A to the sea
)


# Initial state

## In colony A
N <- 1000                  # Total population in colony A
initial_infected <- 1
initial_exposed <- 0
initial_recovered <- 0
initial_susceptible <- N - initial_infected - initial_exposed - initial_recovered
initial_dead <- 0

initial_state_A <- c(S = initial_susceptible,
                     E = initial_exposed,
                     I = initial_infected,
                     R = initial_recovered,
                     D = initial_dead)

## At sea
N <- 50                  # Total population at sea
initial_infected <- 0
initial_exposed <- 0
initial_recovered <- 0
initial_susceptible <- N - initial_infected - initial_exposed - initial_recovered
initial_dead <- 0

initial_state_sea <- c(S = initial_susceptible,
                       E = initial_exposed,
                       I = initial_infected,
                       R = initial_recovered,
                       D = initial_dead)



initial_state = matrix(data = c(initial_state_A,
                                initial_state_sea), 
                       nrow = 2, ncol = 5, 
                       byrow = T)

# Gillespie SEIR model function
gillespie_seir <- function(param, initial_state, total_time) {
  
  sigma = param$sigma
  gamma = param$gamma
  eta = param$eta
  mu = param$mu
  beta_colony = param$beta[1,]
  beta_sea = param$beta[2,]
  zeta = param$zeta

  times <- c(0)
  states <- array(dim = c(2,5,1), data = initial_state)
  
  while (times[length(times)] < total_time) {
    
    period = findInterval(times[length(times)], c(0, 8, 16, 30)) 
    
    beta_colony_t = beta_colony[period]
    beta_sea_t = beta_sea[period]
    
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
      "S_a_to_E_a" = beta_colony_t * S_a * I_a,
      "E_a_to_S_a" = eta * E_a,
      "E_a_to_I_a" = sigma * E_a,
      "I_a_to_R_a" = gamma * I_a,
      "I_a_to_D_a" = mu * I_a,
      
      "S_sea_to_E_sea" = beta_sea_t * S_sea * I_sea,
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
    
    time_step <- rexp(1, total_rate)
    times <- c(times, times[length(times)] + time_step)
    
    transition <- sample(names(rates), 1, prob = rates / total_rate)
    
    if (transition == "S_a_to_E_a") {
      S_a <- S_a - 1
      E_a <- E_a + 1
    } else if (transition == "E_a_to_S_a") {
      E_a <- E_a - 1
      S_a <- S_a + 1
    } else if (transition == "E_a_to_I_a") {
      E_a <- E_a - 1
      I_a <- I_a + 1
    } else if (transition == "I_a_to_R_a") {
      I_a <- I_a - 1
      R_a <- R_a + 1
    } else if (transition == "I_a_to_D_a") {
      I_a <- I_a - 1
      D_a <- D_a + 1
    } else if (transition == "S_sea_to_E_sea") {
      S_sea <- S_sea - 1
      E_sea <- E_sea + 1
    } else if (transition == "E_sea_to_S_sea") {
      E_sea <- E_sea - 1
      S_sea <- S_sea + 1
    } else if (transition == "E_sea_to_I_sea") {
      E_sea <- E_sea - 1
      I_sea <- I_sea + 1
    } else if (transition == "I_sea_to_R_sea") {
      I_sea <- I_sea - 1
      R_sea <- R_sea + 1
    } else if (transition == "I_sea_to_D_sea") {
      I_sea <- I_sea - 1
      D_sea <- D_sea + 1
    } else if (transition == "S_a_to_S_sea") {
      S_a <- S_a - 1
      S_sea <- S_sea + 1
    } else if (transition == "E_a_to_E_sea") {
      E_a <- E_a - 1
      E_sea <- E_sea + 1
    } else if (transition == "I_a_to_I_sea") {
      I_a <- I_a - 1
      I_sea <- I_sea + 1
    } else if (transition == "R_a_to_R_sea") {
      R_a <- R_a - 1
      R_sea <- R_sea + 1
    }
    
    
    new_state = matrix(data = c(S_a, E_a, I_a, R_a, D_a,
                                S_sea, E_sea, I_sea, R_sea, D_sea), 
                           nrow = 2, ncol = 5, 
                           byrow = T)
    
    states <- abind(states, new_state)
  }
  
  output <- data.frame(
    time = times,
    S_a = states[1, 1, ],
    E_a = states[1, 2, ],
    I_a = states[1, 3, ],
    R_a = states[1, 4, ],
    D_a = states[1, 5, ],
    S_sea = states[2, 1, ],
    E_sea = states[2, 2, ],
    I_sea = states[2, 3, ],
    R_sea = states[2, 4, ],
    D_sea = states[2, 5, ]
  )
  
  return(output)
}

# Run simulation

output <- gillespie_seir(param, initial_state, total_time)



# Plot results
output_long <- melt(output, id = "time")
output_a <- output_long %>% filter(variable %in% c("S_a", "E_a", "I_a", "R_a", "D_a"))
plot_a <- ggplot(output_a, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("")

output_sea <- output_long %>% filter(variable %in% c("S_sea", "E_sea", "I_sea", "R_sea", "D_sea"))
plot_sea <- ggplot(output_sea, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("")

# Display plots
plot_a 
plot_sea


# 
# summary_output = function(output){
#   
#   N_a = output[1, c("S_a", "I_a")] %>% sum()
#   
#   max_infected_a = max(output[, "I_a"]) 
#   prop_max_infected_a = max_infected_a / N_a
#   dead_a = output[nrow(output), "D_a"]
#   prop_dead_a = dead_a / N_a
#   
#   max_infected_sea = max(output[, "I_sea"]) 
#   dead_sea = output[nrow(output), "D_sea"]
#   
#   prop_non_exposed_from_a = (output[nrow(output), "S_a"] + output[nrow(output), "S_sea"] - output[1, "S_sea"] ) / N_a
#   
#   return( data.frame(
#     N_a = N_a,
#     max_infected_a = max_infected_a,
#     prop_max_infected_a = prop_max_infected_a,
#     dead_a = dead_a,
#     prop_dead_a = prop_dead_a,
#     max_infected_sea = max_infected_sea,
#     dead_sea = dead_sea,
#     prop_non_exposed_from_a = prop_non_exposed_from_a
#     
#   ))
# }
# 
# 
# output_long_list = data.frame()
# response_list = data.frame()
# 
# nb_iterations = 8
# 
# for (i in 1:nb_iterations){
#   
#   output = gillespie_seir(param, initial_state, total_time)
#   output_long = melt(output, id = "time")
#   
#   output_long_i <- cbind(output_long,
#                          data.frame(simulation = rep(i, times = nrow(output_long))))
#   
#   output_long_list = rbind(output_long_list, output_long_i)
#   response_list = rbind(response_list, summary_output(output))
#   
# }
# 
# 
# output_a <- output_long_list %>% filter(variable %in% c("S_a", "E_a", "I_a", "R_a", "D_a"))
# 
# p = ggplot()
# for(i in 1:nb_iterations){
#   p = p + geom_line(data = output_a %>% subset(., simulation == i )
#                     , aes(x = time, y = value, color = variable)) 
# }
# p = p +
#   labs(x = "Time", y = "Number of individuals", color = "Compartment") +
#   theme_minimal() +
#   ggtitle("Stochastic SEIR Model Simulation (Gillespie Algorithm)")
# 
# 
# p
# 
# 
# data_long <- pivot_longer(response_list, cols = -N_a, names_to = "variable", values_to = "value")
# 
# # Créer les diagrammes en violon pour chaque variable
# ggplot(data_long %>% subset(., variable %in% c("dead_a", "max_infected_a")),
#        aes(x = variable, y = value)) +
#   geom_violin() + 
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# 
# 
