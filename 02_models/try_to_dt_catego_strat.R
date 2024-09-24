

# Exemple de dataframe avec variables explicatives, réponses et scénario
df <- data.frame(
  var_expl1 = c(1, 2, 1, 2),
  var_expl2 = c(3, 3, 4, 4),
  reponse1 = c(10, 12, 15, 18),
  reponse2 = c(20, 22, 25, 28),
  scenario = c("A", "A", "B", "B")
)

library(tidyr)
library(dplyr)

selected_simulation_dt = 

BO_dt = simulation_dt[, 1:30] %>% 
  subset(., scenario == "BO") %>% 
  dplyr::select(-scenario)

RS_dt = simulation_dt[, 1:30] %>% 
  subset(., scenario == "RS") %>% 
  dplyr::select(-scenario)

BO_dt = simulation_dt[, 1:30] %>% 
  subset(., scenario == "") %>% 
  dplyr::select(-scenario)

merge(data_binned_1, data_binned_2, by = c("x_mid", "y_mid"))

# Étape 1 : Pivot pour mettre les scénarios côte à côte
df_wide <- selected_simulation_dt %>%
  pivot_wider(names_from = scenario, values_from = c(nb_adults_equi, nb_infected_colonies, infected_X_time))

# Étape 2 : Automatiser les différences entre tous les scénarios
# Obtenir les noms de colonnes d'outputs pour chaque scénario
output_columns <- colnames(df_wide)[grepl("output", colnames(df_wide))]

# Identifier les différents scénarios
scenarios <- unique(simulation_dt$scenario)

# Calcul dynamique des différences pour chaque output et scénario
for (output in unique(gsub("_.*", "", output_columns))) {
  for (i in 1:(length(scenarios) - 1)) {
    for (j in (i + 1):length(scenarios)) {
      # Obtenir les colonnes des scénarios à comparer
      col1 <- paste0(output, "_", scenarios[i])
      col2 <- paste0(output, "_", scenarios[j])
      
      # Créer la colonne de différence de manière dynamique
      diff_col <- paste0("diff_", output, "_", scenarios[j], "_vs_", scenarios[i])
      
      # Ajouter la colonne de différence
      df_wide <- df_wide %>%
        mutate(!!diff_col := !!sym(col2) - !!sym(col1))
    }
  }
}

# Afficher le résultat final avec les colonnes de différences dynamiques
print(df_wide)
