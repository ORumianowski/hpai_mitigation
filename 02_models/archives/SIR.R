# Paramètres du modèle
population <- 1000  # Taille de la population totale
infectes_initiaux <- 1  # Nombre initial d'infectés
taux_transmission <- 0.3  # Taux de transmission
taux_recuperation <- 0.1  # Taux de récupération
duree_simulation <- 200  # Durée de la simulation en jours

# Initialisation des compartiments S, I et R
S <- population - infectes_initiaux
I <- infectes_initiaux
R <- 0

# Vecteurs pour stocker les résultats au fil du temps
temps <- c(0)
susceptibles <- c(S)
infectes <- c(I)
recuperes <- c(R)

# Boucle de simulation avec l'algorithme de Gillespie
while (temps[length(temps)] < duree_simulation & I!=0) {
  # Calcul des taux d'événements
  taux_infection <- taux_transmission * S * I / population
  taux_guerison <- taux_recuperation * I
  
  # Calcul du taux total d'événements
  taux_total <- taux_infection + taux_guerison
  
  # Génération d'un temps d'événement exponentiel
  temps_evenement <- rexp(1, rate = taux_total)
  
  # Mise à jour du temps de simulation
  temps <- c(temps, temps[length(temps)] + temps_evenement)
  
  # Choix de l'événement (infection ou guérison) en fonction des probabilités
  prob_infection <- taux_infection / taux_total
  if (runif(1) < prob_infection) {
    S <- S - 1
    I <- I + 1
  } else {
    I <- I - 1
    R <- R + 1
  }
  
  # Stockage des résultats
  susceptibles <- c(susceptibles, S)
  infectes <- c(infectes, I)
  recuperes <- c(recuperes, R)
}

# Tracé des courbes
plot(temps, susceptibles, type = 'l', col = 'blue', xlab = 'Jours', ylab = 'Nombre d\'individus',
     main = 'Modèle SIR Stochastique avec l\'algorithme de Gillespie')
lines(temps, infectes, col = 'red')
lines(temps, recuperes, col = 'green')
legend('topright', legend = c('Susceptibles', 'Infectés', 'Récupérés'), col = c('blue', 'red', 'green'), lty = 1)
grid()