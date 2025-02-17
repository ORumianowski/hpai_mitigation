# Parameter ranges --------------------------------------------------------

param_ranges = list(
  
  # Integer parameters
  initial_number_infected_breeders_A = list(range = c(0, 3), distribution = "integer"),
  initial_number_infected_breeders_B = list(range = c(0, 0), distribution = "integer"),
  initial_number_infected_breeders_C = list(range = c(0, 0), distribution = "integer"),
  initial_number_breeders_A = list(range = c(50, 50), distribution = "integer"),
  initial_number_breeders_B = list(range = c(80, 80), distribution = "integer"),
  initial_number_breeders_C = list(range = c(20, 20), distribution = "integer"),
  dispersal_reaction_time =  list(range = c(1, 6), distribution = "integer"),
  dispersal_date = list(range = c(0, 0), distribution = "integer"),
  hatching_date = list(range = c(0, 20), distribution = "integer"),
  
  # Continuous parameters
  tau = list(range = c(0.15, 0.15), distribution = "simple_uniform"),
  total_time = list(range = c(70, 70), distribution = "simple_uniform"),
  prop_dispersal = list(range = c(0.9, 1), distribution = "simple_uniform"),
  beta_E_colony = list(range = c(0, 0), distribution = "simple_uniform"),
  beta_I_colony = list(range = c(0.01, 0.80), distribution = "logarithmic"),
  incubation_period = list(range = c(0.8, 1.2), distribution = "simple_uniform"),
  eta = list(range = c(0, 0), distribution = "simple_uniform"),
  infectious_period = list(range = c(5, 7), distribution = "logarithmic"),
  adult_mortality = list(range = c(0.4,0.8), distribution = "simple_uniform"),
  nestling_mortality = list(range = c(0.6, 1.0), distribution = "simple_uniform"),
  avrg_stay_B_colony = list(range = c(2, 2), distribution = "logarithmic"),
  avrg_stay_B_sea = list(range = c(2, 2), distribution = "logarithmic"),
  avrg_stay_NB_colony = list(range = c(2, 2), distribution = "logarithmic"),
  avrg_stay_NB_sea = list(range = c(2, 40), distribution = "logarithmic"),
  theta = list(range = c(1/100, 1/7), distribution = "simple_uniform"),
  psi = list(range = c(1/500, 1/500), distribution = "logarithmic"),
  hatching_sd = list(range = c(3, 3), distribution = "simple_uniform"),
  reaching_repro_prob = list(range = c(0.3, 0.7), distribution = "simple_uniform"),
  prob_detection = list(range = c(0.7, 1.0), distribution = "simple_uniform")
)
