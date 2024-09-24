

library(rpart)				        # Popular decision tree algorithm
library(rattle)					# Fancy tree plot
library(rpart.plot)				# Enhanced tree plots
library(RColorBrewer)				# Color selection for fancy tree plot
library(party)					# Alternative decision tree algorithm
library(partykit)				# Convert rpart object to BinaryTree
library(caret)	



model <- rpart(nb_adults_equi ~ beta_I_colony + reaching_repro_prob + initial_number_infected_breeders_A,
               data = simulation_dt)
prp(model)
printcp(model)
plotcp(model)
model$cptable[which.min(model$cptable[,"xerror"]),"CP"]
model.prune=prune(model,cp=0.04180515) #I manually selected the one that had the lowest CP and cross-validation error


#Caret with CV
train <- createDataPartition(simulation_dt[,"nb_adults_equi"],p=0.70,list=FALSE) # nb_adults_equi is the "target" class
simulation_dt.trn <- simulation_dt[train,] # Check with str that nrows of simulation_dt.trn is 70% of the original "simulation_dt"
simulation_dt.tst <- simulation_dt[-train,] # And this, just 30%
ctrl  <- trainControl(method  = "cv",number  = 10)
nb_adults_equi <- as.factor(simulation_dt[,"nb_adults_equi"])
simulation_dt.trn$nb_adults_equi<- as.factor(simulation_dt.trn$nb_adults_equi)
fit.cv <- train(nb_adults_equi ~ beta_I_colony + reaching_repro_prob + initial_number_infected_breeders_A,
                data = simulation_dt.trn, method = "rpart", # k nearest neighbors
                trControl = ctrl,  # Add the control
                # preProcess = c("center","scale"),  # preprocess the simulation_dt (center=> -mean(); scale= /standard.deviation)
                #tuneGrid =data.frame(cp=seq(0,1,by=0.0025))) # Try only these values in the CV step
                tuneLength = 100)
# tuneLength = 25) # Use 25 sequential numbers instead
pred <- predict(fit.cv,simulation_dt.tst) # predict the output classes
# pred.prb <- predict(fit.cv,simulation_dt.tst, type="prob") # predict the output classes


confusionMatrix(table(simulation_dt.tst[,"nb_adults_equi"],pred))
print(fit.cv) # Plot the results. See at the end the chosen value of "k"
plot(fit.cv) # Plot the Cross-validation output
plot(fit.cv$finalModel)
text(fit.cv$finalModel)
library(rpart.plot)
rpart.plot(fit.cv$finalModel, fallen.leaves = F)

