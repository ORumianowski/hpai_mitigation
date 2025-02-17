

setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")


library(ppcor)
library(GGally)
library(ggplot2)
library(tidyverse)
library(forcats)

load("simulation_dt/simulation_dt_1000it_cluster_1.RData")
simulation_dt1 = simulation_dt
load("simulation_dt/simulation_dt_1000it_ordi_1.RData")
simulation_dt2 = simulation_dt
load("simulation_dt/simulation_dt_10000it_cluster_1.RData")
simulation_dt3 = simulation_dt
load("simulation_dt/simulation_dt_2000it_cluster_1.RData")
simulation_dt4 = simulation_dt
simulation_dt = rbind(simulation_dt1,
                      simulation_dt2,
                      simulation_dt3,
                      simulation_dt4)


param = c("initial_number_infected_breeders_A", "initial_number_infected_breeders_B", "initial_number_infected_breeders_C", "initial_number_breeders_A",         
          "initial_number_breeders_B",          "initial_number_breeders_C",          "dispersal_reaction_time",            "dispersal_date",                    
          "hatching_date",                      "tau",                                "total_time",                         "prop_dispersal",                    
          "beta_E_colony",                      "beta_I_colony",                      "incubation_period",                  "eta",                               
          "infectious_period",                  "adult_mortality",                    "nestling_mortality",                 "avrg_stay_B_colony",                
          "avrg_stay_B_sea",                    "avrg_stay_NB_colony",                "avrg_stay_NB_sea",                   "theta",                             
          "psi",                                "hatching_sd",                        "reaching_repro_prob",                "prob_detection"  )

param_evaluated = c("initial_number_infected_breeders_A",        
                    "hatching_date",                      "prop_dispersal",                    
                    "beta_I_colony",                      "incubation_period",                                              
                    "infectious_period",                  "adult_mortality",   
                    "nestling_mortality",                         
                    "avrg_stay_NB_sea",                   "theta",                             
                    "reaching_repro_prob",                "prob_detection",
                    
                    "dispersal_reaction_time",
                    "prob_detection")



dt <- simulation_dt %>%
  dplyr::select(1:32) %>% 
  pivot_wider(
    names_from = scenario,  
    values_from = c(nb_adults_equi, nb_infected_colonies, infected_X_time), 
    names_glue = "{scenario}_{.value}" 
  ) %>% 
  unnest(cols = everything()) %>% 
  mutate(BO_P2_nb_adults_equi = BO_nb_adults_equi - P2_nb_adults_equi) %>% 
  mutate(BO_RS_nb_adults_equi = BO_nb_adults_equi - RS_nb_adults_equi) %>% 
  mutate(BO_P2_infected_X_time = BO_infected_X_time - P2_infected_X_time) %>% 
  mutate(BO_RS_infected_X_time = BO_infected_X_time - RS_infected_X_time)


partial_cor = function(data, param_evaluated, selected_output){
    
    dt_cor = data %>% 
      dplyr::select(all_of(c(param_evaluated, selected_output)))

    pcor_results  = pcor(dt_cor)

    pcor_matrix <- pcor_results$estimate %>% 
      as.matrix() %>% 
      data.frame() %>% 
      dplyr::select(selected_output) %>% 
      arrange(desc(!!sym(selected_output))) 
    
    return(pcor_matrix)
  }
  

res_ENLA = cbind(
partial_cor(data = dt,
            param_evaluated = param_evaluated,
            selected_output = "BO_nb_adults_equi")
,
partial_cor(data = dt,
            param_evaluated = param_evaluated,
            selected_output = "BO_P2_nb_adults_equi")
,
partial_cor(data = dt,
            param_evaluated = param_evaluated,
            selected_output = "BO_RS_nb_adults_equi")
) %>% 
  as.data.frame() %>% 
  slice(-1) %>%
  rownames_to_column(var = "param") %>% 
  pivot_longer(-param) %>%
  mutate(param = fct_reorder(param, value, .desc = TRUE))


res_infected = cbind(
  partial_cor(data = dt,
              param_evaluated = param_evaluated,
              selected_output = "BO_infected_X_time")
  ,
  partial_cor(data = dt,
              param_evaluated = param_evaluated,
              selected_output = "BO_P2_infected_X_time")
  ,
  partial_cor(data = dt,
              param_evaluated = param_evaluated,
              selected_output = "BO_RS_infected_X_time")
) %>% 
  as.data.frame() %>% 
  slice(-1) %>%
  rownames_to_column(var = "param") %>% 
  pivot_longer(-param) %>%
  mutate(param = fct_reorder(param, value, .desc = TRUE))





ggplot() +
  geom_point(data = res_ENLA, aes(x = param, y = value)) +
  facet_wrap(~ name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Pour incliner les labels des abscisses si nécessaire


ggplot() +
  geom_point(data = res_infected, aes(x = param, y = value)) +
  facet_wrap(~ name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Pour incliner les labels des abscisses si nécessaire



