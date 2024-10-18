
setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")

library(tidyverse)
library(rpart)				       
library(rattle)					
library(rpart.plot)			
library(RColorBrewer)				
library(party)					
library(partykit)			
library(caret)	

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
  mutate(BO_RS_infected_X_time = BO_infected_X_time - RS_infected_X_time) %>% 
  mutate(P2_good = BO_P2_nb_adults_equi > 0)%>% 
  mutate(P2_good = as.factor(P2_good)) %>% 
  mutate(RS_good = BO_RS_nb_adults_equi > 0)%>% 
  mutate(RS_good = as.factor(RS_good))

####################   Caret with CV

train <- createDataPartition(1:nrow(dt),p=0.80,list=FALSE) 
dt.trn <- dt[train,]
dt.tst <- dt[-train,] 
ctrl  <- trainControl(method  = "cv", number  = 50)

fit.cv <- train(RS_good ~  
                  initial_number_infected_breeders_A +
                  hatching_date +                      
                  prop_dispersal +                    
                  beta_I_colony +
                  incubation_period +                                           
                  infectious_period +                 
                  adult_mortality+   
                  nestling_mortality +                         
                  avrg_stay_NB_sea +                 
                  theta +                  
                  reaching_repro_prob +          
                  prob_detection,
                
                data = dt.trn,
                method = "rpart", 
                trControl = ctrl,  
                tuneLength = 6)

# Plot the final tree model
par(xpd = NA) # Avoid clipping the text in some device
plot(fit.cv$finalModel)
text(fit.cv$finalModel,  digits = 5)


# pred = predict(fit.cv,dt.tst) %>% 
#   as.numeric()
# 
# real = dt.tst[,"nb_adults_equi"]%>% 
#   as.numeric()
# 
# confusionMatrix(table(real,pred))

print(fit.cv) # Plot the results. See at the end the chosen value of "k"
plot(fit.cv) # Plot the Cross-validation output
# plot(fit.cv$finalModel)
# text(fit.cv$finalModel)

rpart.plot(fit.cv$finalModel, fallen.leaves = T)



# Basic method  -----------------------------------------------------------


# 
# model <- rpart(P2_good ~  
#                  initial_number_infected_breeders_A +
#                  hatching_date +                      
#                  prop_dispersal +                    
#                  beta_I_colony +
#                  incubation_period +                                           
#                  infectious_period +                 
#                  adult_mortality+   
#                  nestling_mortality +                         
#                  avrg_stay_NB_sea +                 
#                  theta +                  
#                  reaching_repro_prob +          
#                  prob_detection,
#                
#                data = dt)
# prp(model)
# 
# 
# 
# # printcp(model)
# # plotcp(model)
# # model$cptable[which.min(model$cptable[,"xerror"]),"CP"]
# # model.prune=prune(model,cp=0.04180515) #I manually selected the one that had the lowest CP and cross-validation error
# 
# 
# 
# 
# 
# # Best strategy -----------------------------------------------------------
# 
# 
# 
# 
# X = simulation_dt %>% 
#   subset(., scenario  == "BO") %>% 
#   dplyr::select(1:26)
# 
# Y_BO = simulation_dt %>% 
#   subset(., scenario  == "BO") %>% 
#   dplyr::select(28:30) %>% 
#   rename(nb_adults_equi = nb_adults_equi_BO)
# 
# Y_RS = simulation_dt %>% 
#   subset(., scenario  == "RS") %>% 
#   dplyr::select(28:30)
# 
# cbind(X, Y_BO, Y_RS)
# 
# 
# 
# # 
# # dt = best_strat_dt
# # 
# # 
# # ####################Caret with CV
# # train <- createDataPartition(1:nrow(dt),p=0.80,list=FALSE) # nb_adults_equi is the "target" class
# # dt.trn <- dt[train,] # Check with str that nrows of dt.trn is 70% of the original "dt"
# # dt.tst <- dt[-train,] # And this, just 30%
# # ctrl  <- trainControl(method  = "cv",number  = 10)
# # 
# # #nb_adults_equi <- as.factor(dt[,"nb_adults_equi"])
# # 
# # # dt.trn$nb_adults_equi<- as.factor(dt.trn$nb_adults_equi)
# # 
# # fit.cv <- train(best_strat ~ beta_I_colony 
# #                 + reaching_repro_prob 
# #                 + rho_to_colony
# #                 + initial_number_infected_breeders_A
# #                 + dispersal_reaction_time,
# #                 data = dt.trn, method = "rpart", # k nearest neighbors
# #                 trControl = ctrl,  
# #                 tuneLength = 100) # ??
# # 
# # 
# # ### or categoriacal output
# # # pred = predict(fit.cv,dt.tst) %>% 
# # #   as.numeric()
# # # 
# # # real = dt.tst[,"nb_adults_equi"]%>% 
# # #   as.numeric()
# # # 
# # # confusionMatrix(table(real,pred))
# # 
# # print(fit.cv) 
# # plot(fit.cv) 
# # library(rpart.plot)
# # rpart.plot(fit.cv$finalModel, fallen.leaves = F)
# # 
