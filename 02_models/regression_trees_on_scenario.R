
setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")

library(rpart)				        # Popular decision tree algorithm
library(rattle)					# Fancy tree plot
library(rpart.plot)				# Enhanced tree plots
library(RColorBrewer)				# Color selection for fancy tree plot
library(party)					# Alternative decision tree algorithm
library(partykit)				# Convert rpart object to BinaryTree
library(caret)	


load("simulation_dt/simulation_dt_50_2.RData")
simulation_dt1 = simulation_dt
load("simulation_dt/simulation_dt_150_2.RData")
simulation_dt2 = simulation_dt
load("simulation_dt/simulation_dt_300_2.RData")
simulation_dt3 = simulation_dt
load("simulation_dt/simulation_dt_500_2.RData")
simulation_dt4 = simulation_dt

simulation_dt = rbind(simulation_dt1,
                      simulation_dt2,
                      simulation_dt3,
                      simulation_dt4)


selected_simulation_dt = simulation_dt %>% subset(., scenario == "RS")

model <- rpart(nb_adults_equi ~ beta_I_colony 
               + reaching_repro_prob 
               + rho_to_colony
               + initial_number_infected_breeders_A
               + dispersal_reaction_time,
               data = selected_simulation_dt)
prp(model)



# printcp(model)
# plotcp(model)
# model$cptable[which.min(model$cptable[,"xerror"]),"CP"]
# model.prune=prune(model,cp=0.04180515) #I manually selected the one that had the lowest CP and cross-validation error



####################Caret with CV
train <- createDataPartition(1:nrow(selected_simulation_dt),p=0.80,list=FALSE) # nb_adults_equi is the "target" class
selected_simulation_dt.trn <- selected_simulation_dt[train,] # Check with str that nrows of selected_simulation_dt.trn is 70% of the original "selected_simulation_dt"
selected_simulation_dt.tst <- selected_simulation_dt[-train,] # And this, just 30%
ctrl  <- trainControl(method  = "cv",number  = 10)

#nb_adults_equi <- as.factor(selected_simulation_dt[,"nb_adults_equi"])

# selected_simulation_dt.trn$nb_adults_equi<- as.factor(selected_simulation_dt.trn$nb_adults_equi)

fit.cv <- train(nb_adults_equi ~ beta_I_colony 
                + reaching_repro_prob 
                + rho_to_colony
                + initial_number_infected_breeders_A
                + dispersal_reaction_time,
                data = selected_simulation_dt.trn, method = "rpart", # k nearest neighbors
                trControl = ctrl,  
                tuneLength = 100) # ??


### or categoriacal output
# pred = predict(fit.cv,selected_simulation_dt.tst) %>% 
#   as.numeric()
# 
# real = selected_simulation_dt.tst[,"nb_adults_equi"]%>% 
#   as.numeric()
# 
# confusionMatrix(table(real,pred))

print(fit.cv) # Plot the results. See at the end the chosen value of "k"
plot(fit.cv) # Plot the Cross-validation output
# plot(fit.cv$finalModel)
# text(fit.cv$finalModel)
library(rpart.plot)
rpart.plot(fit.cv$finalModel, fallen.leaves = F)



# Best strategy -----------------------------------------------------------




X = simulation_dt %>% 
  subset(., scenario  == "BO") %>% 
  dplyr::select(1:26)

Y_BO = simulation_dt %>% 
  subset(., scenario  == "BO") %>% 
  dplyr::select(28:30) %>% 
  rename(nb_adults_equi = nb_adults_equi_BO)

Y_RS = simulation_dt %>% 
  subset(., scenario  == "RS") %>% 
  dplyr::select(28:30)

cbind(X, Y_BO, Y_RS)



# 
# selected_simulation_dt = best_strat_dt
# 
# 
# ####################Caret with CV
# train <- createDataPartition(1:nrow(selected_simulation_dt),p=0.80,list=FALSE) # nb_adults_equi is the "target" class
# selected_simulation_dt.trn <- selected_simulation_dt[train,] # Check with str that nrows of selected_simulation_dt.trn is 70% of the original "selected_simulation_dt"
# selected_simulation_dt.tst <- selected_simulation_dt[-train,] # And this, just 30%
# ctrl  <- trainControl(method  = "cv",number  = 10)
# 
# #nb_adults_equi <- as.factor(selected_simulation_dt[,"nb_adults_equi"])
# 
# # selected_simulation_dt.trn$nb_adults_equi<- as.factor(selected_simulation_dt.trn$nb_adults_equi)
# 
# fit.cv <- train(best_strat ~ beta_I_colony 
#                 + reaching_repro_prob 
#                 + rho_to_colony
#                 + initial_number_infected_breeders_A
#                 + dispersal_reaction_time,
#                 data = selected_simulation_dt.trn, method = "rpart", # k nearest neighbors
#                 trControl = ctrl,  
#                 tuneLength = 100) # ??
# 
# 
# ### or categoriacal output
# # pred = predict(fit.cv,selected_simulation_dt.tst) %>% 
# #   as.numeric()
# # 
# # real = selected_simulation_dt.tst[,"nb_adults_equi"]%>% 
# #   as.numeric()
# # 
# # confusionMatrix(table(real,pred))
# 
# print(fit.cv) 
# plot(fit.cv) 
# library(rpart.plot)
# rpart.plot(fit.cv$finalModel, fallen.leaves = F)
# 
