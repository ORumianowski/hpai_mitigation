

setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")


# Downloading model -------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)

source("model.R")

# LHS ---------------------------------------------------------------------

library(lhs)
library(ggplot2)
library(dplyr)
library(stringr)

source("param_ranges.R")
source("scenarios.R")

# Convert LHS samples  ----------------------------------------------------

convert_samples = function(lhs_samples, param_ranges, nb_params) {
  for (i in 1:nb_params) {
    if (param_ranges[[i]]$distribution == "integer"){
      lhs_samples[, i] = qinteger(lhs_samples[, i], param_ranges[[i]]$range[1], param_ranges[[i]]$range[2])
    }
    else if (param_ranges[[i]]$distribution == "logarithmic"){
      log_min <- log(param_ranges[[i]]$range[1])
      log_max <- log(param_ranges[[i]]$range[2])
      lhs_samples[, i] <- exp(qunif(lhs_samples[, i], log_min, log_max))
    }
    else if (param_ranges[[i]]$distribution == "simple_uniform"){
      lhs_samples[, i] = qunif(lhs_samples[, i], param_ranges[[i]]$range[1], param_ranges[[i]]$range[2])
    }
  }
  return(lhs_samples)
}




# Number of samples
nb_samples = 100

# Total number of parameters 
nb_params = length(param_ranges)

# Generate samples
lhs_samples = randomLHS(nb_samples, nb_params)
samples = convert_samples(lhs_samples, param_ranges, nb_params)
colnames(samples) = names(param_ranges)


induced_dispersal_ = F
dispersal_stochastic_ = F
initially_infected_ = T

library(rsample)

samples_split <- initial_split(samples %>% data.frame(), prop = 0.5)
train_samples <- training(samples_split)
test_samples <- testing(samples_split)

# sensitivity analysis

x <- soboljansen(model = gillespie_seir, X1 = train_samples, X2 = test_samples, nboot = nb_samples/2)
print(x)
plot(x)


