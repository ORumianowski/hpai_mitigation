

setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")


library(ppcor)
library(GGally)
library(ggplot2)
library(tidyverse)

load("simulation_dt/simulation_dt_1000it_cluster_1.RData")
simulation_dt1 = simulation_dt
load("simulation_dt/simulation_dt_1000it_ordi_1.RData")
simulation_dt2 = simulation_dt
simulation_dt = rbind(simulation_dt1,
                      simulation_dt2)


param = c("initial_number_infected_breeders_A", "initial_number_infected_breeders_B", "initial_number_infected_breeders_C", "initial_number_breeders_A",         
          "initial_number_breeders_B",          "initial_number_breeders_C",          "dispersal_reaction_time",            "dispersal_date",                    
          "hatching_date",                      "tau",                                "total_time",                         "prop_dispersal",                    
          "beta_E_colony",                      "beta_I_colony",                      "incubation_period",                  "eta",                               
          "infectious_period",                  "adult_mortality",                    "nestling_mortality",                 "avrg_stay_B_colony",                
          "avrg_stay_B_sea",                    "avrg_stay_NB_colony",                "avrg_stay_NB_sea",                   "theta",                             
          "psi",                                "hatching_sd",                        "reaching_repro_prob",                "prob_detection"  )

param_evaluated = c("BO_P2_nb_adults_equi",
                    "initial_number_infected_breeders_A",        
                    "hatching_date",                      "prop_dispersal",                    
                    "beta_I_colony",                      "incubation_period",                                              
                    "infectious_period",                  "adult_mortality",   
                    "nestling_mortality",                         
                    "avrg_stay_NB_sea",                   "theta",                             
                    "reaching_repro_prob",                "prob_detection"  )

dt = simulation_dt[,1:32]


dt2 <- dt %>%
  pivot_wider(
    names_from = scenario,  
    values_from = c(nb_adults_equi, nb_infected_colonies, infected_X_time), 
    names_glue = "{scenario}_{.value}" 
  ) %>% 
  unnest(cols = everything()) %>% 
  mutate(BO_P2_nb_adults_equi = BO_nb_adults_equi - P2_nb_adults_equi) %>% 
  mutate(BO_RS_nb_adults_equi = BO_nb_adults_equi - RS_nb_adults_equi)


dt3 = dt2 %>% 
  dplyr::select(c(param, "BO_P2_nb_adults_equi"))%>% 
  dplyr::select(param_evaluated)
  
pcor_results  = pcor(dt3)

pcor_matrix <- pcor_results$estimate %>% 
  data.frame() %>% 
  dplyr::select(BO_P2_nb_adults_equi) %>% 
  arrange(desc(BO_P2_nb_adults_equi))
