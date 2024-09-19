# Describe the 5 scenario
scenarios = data.frame(
  induced_dispersal = c(F,F,T,T,T),
  dispersal_stochastic = c(F,F,F,F,T),
  initially_infected = c(F,T,F,T,T))

rownames(scenarios) = c("HS","BO","PS","P2","RS")