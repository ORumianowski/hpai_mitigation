
# Charger les bibliothèques nécessaires
library(ggplot2)

# Paramètres
rho <- 0.29
phi_juv <- 0.5
phi_ad <- 0.9
kappa <- 0.25

# Conditions initiales
B0 <- 100  # Population initiale des adultes 
N0 <- 50   # Population initiale des jeunes
T <- 50   # Nombre d'itérations

# Initialisation des vecteurs
B <- numeric(T)
N <- numeric(T)

# Assignation des valeurs initiales
B[1] <- B0
N[1] <- N0


# Impact ponctuelle
date = 20




# Itération des suites
for (t in 1:(T - 1)) {
  N[t + 1] <- B[t] * rho * phi_juv + N[t] * phi_ad * (1 - kappa)
  B[t + 1] <- B[t] * phi_ad + N[t] * phi_ad * kappa
  if (t == date){
    B[t + 1] = B[t + 1] * 0.9
  }
}

# Création d'un data frame pour ggplot2
data <- data.frame(
  Temps = rep(1:T, 2),
  Population = c(N, B),
  Type = rep(c("Jeunes (N)", "Adultes (B)"), each = T)
)

# Tracé avec ggplot2
ggplot(data, aes(x = Temps, y = Population, color = Type)) +
  geom_line(size = 1.2) +
  geom_point() +
  labs(title = "Évolution des populations de jeunes et adultes",
       x = "Temps",
       y = "Population") +
  theme_minimal() +
  theme(legend.title = element_blank())
