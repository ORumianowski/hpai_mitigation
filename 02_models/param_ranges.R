# Parameter ranges --------------------------------------------------------

param_ranges = list(
  
  # Integer parameters
  initial_number_infected_breeders_A = list(range = c(1, 3), distribution = "integer"),
  initial_number_infected_breeders_B = list(range = c(0, 0), distribution = "integer"),
  initial_number_infected_breeders_C = list(range = c(0, 0), distribution = "integer"),
  initial_number_breeders_A = list(range = c(50, 50), distribution = "integer"),
  initial_number_breeders_B = list(range = c(80, 80), distribution = "integer"),
  initial_number_breeders_C = list(range = c(20, 20), distribution = "integer"),
  dispersal_reaction_time =  list(range = c(1, 10), distribution = "integer"),
  dispersal_date = list(range = c(0, 0), distribution = "integer"),
  hatching_date = list(range = c(10, 10), distribution = "integer"),
  
  # Continuous parameters
  tau = list(range = c(0.15, 0.15), distribution = "simple_uniform"),
  total_time = list(range = c(70, 70), distribution = "simple_uniform"),
  prop_dispersal = list(range = c(1, 1), distribution = "simple_uniform"),
  beta_E_colony = list(range = c(0, 0), distribution = "simple_uniform"),
  beta_I_colony = list(range = c(0.05, 0.80), distribution = "simple_uniform"),
  sigma = list(range = c(1/1, 1/1), distribution = "logarithmic"),
  eta = list(range = c(0, 0), distribution = "simple_uniform"),
  gamma = list(range = c(1/6, 1/6), distribution = "logarithmic"),
  mu_adult = list(range = c(1/6 * (0.5 / (1 - 0.5)), 1/6 * (0.5 / (1 - 0.5))), distribution = "logarithmic"),
  mu_nestling = list(range = c(1/6 * (0.8 / (1 - 0.8)), 1/6 * (0.8 / (1 - 0.8))), distribution = "logarithmic"),
  zeta_to_sea = list(range = c(1/2, 1/2), distribution = "logarithmic"),
  zeta_to_colony = list(range = c(1/2, 1/2), distribution = "logarithmic"),
  rho_to_sea = list(range = c(1/2, 1/2), distribution = "logarithmic"),
  rho_to_colony = list(range = c(1/40, 1/2), distribution = "logarithmic"),
  psi = list(range = c(1/500, 1/500), distribution = "logarithmic"),
  hatching_sd = list(range = c(3, 3), distribution = "simple_uniform"),
  reaching_repro_prob = list(range = c(0.3, 0.7), distribution = "simple_uniform")
)
