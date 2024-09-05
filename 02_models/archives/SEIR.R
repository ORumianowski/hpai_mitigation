library(ggplot2)
library(reshape2)


# Parameters
N <- 1000            # Total population
beta <- 0.01          # Transmission rate
sigma <- 1/5.2       # Rate of progression from exposed to infectious (inverse of incubation period)
gamma <- 1/2.9       # Recovery rate (inverse of infectious period)
initial_infected <- 1
initial_exposed <- 0
initial_recovered <- 0
initial_susceptible <- N - initial_infected - initial_exposed - initial_recovered
total_time <- 160    # Simulation time

# Gillespie SEIR model function
gillespie_seir <- function(N, beta, sigma, gamma, initial_state, total_time) {
  times <- c(0)
  states <- matrix(ncol = 4, nrow = 1, data = initial_state)
  
  while (times[length(times)] < total_time) {
    S <- states[nrow(states), 1]
    E <- states[nrow(states), 2]
    I <- states[nrow(states), 3]
    R <- states[nrow(states), 4]
    
    rates <- c(
      "S_to_E" = beta * S * I,
      "E_to_I" = sigma * E,
      "I_to_R" = gamma * I
    )
    
    total_rate <- sum(rates)
    
    if (total_rate == 0) {
      break
    }
    
    time_step <- rexp(1, total_rate)
    times <- c(times, times[length(times)] + time_step)
    
    transition <- sample(names(rates), 1, prob = rates / total_rate)
    
    if (transition == "S_to_E") {
      S <- S - 1
      E <- E + 1
    } else if (transition == "E_to_I") {
      E <- E - 1
      I <- I + 1
    } else if (transition == "I_to_R") {
      I <- I - 1
      R <- R + 1
    }
    
    states <- rbind(states, c(S, E, I, R))
  }
  
  return(list(times = times, states = states))
}

# Initial state
initial_state <- c(S = initial_susceptible, E = initial_exposed, I = initial_infected, R = initial_recovered)

# Run simulation
result <- gillespie_seir(N, beta, sigma, gamma, initial_state, total_time)

# Convert results to data frame for plotting
output <- data.frame(time = result$times,
                     S = result$states[, 1], 
                     E = result$states[, 2], 
                     I = result$states[, 3], 
                     R = result$states[, 4])

output_long <- melt(output, id = "time")

# Plot results
ggplot(output_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "Time", y = "Number of individuals", color = "Compartment") +
  theme_minimal() +
  ggtitle("Stochastic SEIR Model Simulation (Gillespie Algorithm)")
