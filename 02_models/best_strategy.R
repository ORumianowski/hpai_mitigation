


library(rpart)				        # Popular decision tree algorithm
library(rattle)					# Fancy tree plot
library(rpart.plot)				# Enhanced tree plots
library(RColorBrewer)				# Color selection for fancy tree plot
library(party)					# Alternative decision tree algorithm
library(partykit)				# Convert rpart object to BinaryTree
library(caret)	



my_tuneLength = 100

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

x = simulation_dt[, 1:26]

y = data.frame(
  nb_adults_equi_BO = simulation_dt %>% subset(., scenario == "BO") %>% dplyr::select(nb_adults_equi) %>% unlist() %>% as.vector(),
  nb_adults_equi_RS = simulation_dt %>% subset(., scenario == "RS") %>% dplyr::select(nb_adults_equi) %>% unlist() %>% as.vector(),
  nb_adults_equi_P2 = simulation_dt %>% subset(., scenario == "P2") %>% dplyr::select(nb_adults_equi) %>% unlist() %>% as.vector())%>% 
  mutate(best_strat = apply(., 1, function(row) {colnames(.)[which.min(row)]}
  ))
  



dt = cbind(x, y)

train <- createDataPartition(1:nrow(dt),p=0.80,list=FALSE) # nb_adults_equi is the "target" class
dt.trn <- dt[train,] # Check with str that nrows of dt.trn is 70% of the original "dt"
dt.tst <- dt[-train,] # And this, just 30%
ctrl  <- trainControl(method  = "cv",number  = 10)


fit.cv <- train(best_strat ~ beta_I_colony 
                + reaching_repro_prob 
                + rho_to_colony
                + initial_number_infected_breeders_A
                + dispersal_reaction_time,
                data = dt.trn, method = "rpart", # k nearest neighbors
                trControl = ctrl,
                tuneLength = my_tuneLength)

print(fit.cv) # Plot the results. See at the end the chosen value of "k"
plot(fit.cv) # Plot the Cross-validation output
plot(fit.cv$finalModel) #3000 pixels
text(fit.cv$finalModel)
library(rpart.plot)
rpart.plot(fit.cv$finalModel, fallen.leaves = F)
