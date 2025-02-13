
setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")

library(tidyverse)
library(rpart)				       
library(rattle)					
library(rpart.plot)			
library(RColorBrewer)				
library(party)					
library(partykit)			
library(caret)	

#load("simulation_dt/dt_cluster_3_1E4_for_v2.RData")
load("simulation_dt/dt_cluster_3_1E4_for_one_col_v1.RData")


dt = simulation_dt %>% 
  mutate(S_adult = S_a_B + S_b_B + S_c_B + S_a_NB + S_b_NB + S_c_NB) %>% 
  mutate(R_adult = R_a_B + R_b_B + R_c_B + R_a_NB + R_b_NB + R_c_NB) %>% 
  mutate(D_adult = D_a_B + D_b_B + D_c_B + D_a_NB + D_b_NB + D_c_NB) %>% 
  mutate(prop_R = R_adult / (R_adult + S_adult)) %>% 
  mutate(prop_D = D_adult / (D_adult + R_adult + S_adult)) %>% 
  subset(., scenario == "BO") %>% 
  subset(., initial_number_infected_breeders_A != 0) %>% 
  mutate(ad_mortality = cut(adult_mortality,
                            breaks = 3, 
                            labels = c("Low mortality", "Medium mortality", "High mortality")))%>% 
  mutate(transmission_rate = cut(log(beta_I_colony),
                            breaks = 3, 
                            labels = c("Low Transmission Rate", "Medium Transmission Rate", "High Transmission Rate")))




ggplot()+
  geom_point(data = dt, aes(x = prop_D, y = prop_R, color = beta_I_colony )) +
  scale_color_gradient2(trans = "log", low = "blue", mid = "yellow", high = "red", midpoint = log(0.1)) +
  facet_wrap(~ad_mortality)


ggplot()+
  geom_point(data = dt, aes(x = prop_D, y = prop_R, color = adult_mortality )) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0.5) +
  facet_wrap(~transmission_rate)



dt2 = dt %>%
  subset(., prop_R < 0.25)%>%
  subset(., prop_R > 0.00)%>%
  subset(., prop_D < 0.10)%>%
  subset(., prop_D > 0.00)


ggplot()+
  geom_point(data = dt2, aes(x = prop_D, y = prop_R, color = beta_I_colony )) +
  scale_color_gradient2(trans = "log", low = "blue", mid = "yellow", high = "red", midpoint = log(0.1)) +
  facet_wrap(~ad_mortality)


ggplot()+
  geom_point(data = dt2, aes(x = prop_D, y = prop_R, color = adult_mortality )) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0.5) +
  facet_wrap(~transmission_rate)

