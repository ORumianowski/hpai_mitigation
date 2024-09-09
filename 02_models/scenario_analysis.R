

library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)
library(gridExtra)

source("model.R")



# stat_model --------------------------------------------------------------

nb_iterations = 50
#nb_iterations = 25

stat_model = function(nb_iterations_ = nb_iterations,
                      # Number of simulation days
                      total_time_ = 50,
                      # Reaction time between 1rst death and induced dispersal 
                      dispersal_reaction_time_ = 2,
                      # Initial conditions
                      initial_number_infected_breeders_A_ = 1,
                      initial_number_breeders_A_ = 100,
                      initial_number_breeders_B_ = 100,
                      initial_number_breeders_C_ = 10,
                      
                      # Do we induce dispersion ?
                      induced_dispersal_,
                      # Induced dispersion mode (deterministic or stochastic)
                      dispersal_stochactic_,
                      # Transmission rate from exposed individuals and from infectious individuals in a colony
                      BETA_,
                      # Time at sea before returning to a colony (non-breeders)
                      TIME_AT_SEA_NB_,
                      # Parameter of the tau-leap algorithm
                      tau_,
                      # Probability of a nestling becoming a breeder
                      reaching.repro.prob_){
  
  response_list = data.frame()
  
  for (i in 1:nb_iterations_){
    
    output = gillespie_seir(induced_dispersal = induced_dispersal_,
                            dispersal_stochactic = dispersal_stochactic_,
                            dispersal_reaction_time = dispersal_reaction_time_,
                            initial_number_infected_breeders_A = initial_number_infected_breeders_A_,
                            initial_number_breeders_A = initial_number_breeders_A_,
                            initial_number_breeders_B = initial_number_breeders_B_,
                            initial_number_breeders_C = initial_number_breeders_C_,
                            TIME_AT_SEA_NB = TIME_AT_SEA_NB_,
                            BETA = BETA_,
                            total_time = total_time_,
                            tau = tau_)
    
    response_list = rbind(response_list, summary_output(output, reaching.repro.prob_))
    
  }
  return(response_list)
}

# res = stat_model()
# res$nb_infected_colonies



# plot - scenario ----------------------------------------------

scenario_dt = function(beta_context,
                       time_at_sea_NB_context,
                       reaching.repro.prob,
                       tau = 0.1
){
  
  no_stress = 
    stat_model(nb_iterations_ = 1,
               induced_dispersal_ = F,
               initial_number_infected_breeders_A_ = 0,
               tau_ = tau,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context,
               reaching.repro.prob_ = reaching.repro.prob)
  
  baseline_outbreak = 
    stat_model(induced_dispersal_ = F,
               initial_number_infected_breeders_A_ = 1,
               tau_ = tau,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context,
               reaching.repro.prob_ = reaching.repro.prob)
  
  proactive_strategy = 
    stat_model(nb_iterations_ = 1,
               induced_dispersal_ = T,
               initial_number_infected_breeders_A_ = 0,
               dispersal_stochactic_ = F,
               tau_ = tau,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context,
               reaching.repro.prob_ = reaching.repro.prob)
  
  proactive_strategy_toolate = 
    stat_model(induced_dispersal_ = T,
               initial_number_infected_breeders_A_ = 1,
               dispersal_stochactic_ = F,
               tau_ = tau,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context,
               reaching.repro.prob_ = reaching.repro.prob)
  
  reactive_strategy = 
    stat_model(induced_dispersal_ = T,
               initial_number_infected_breeders_A_ = 1,
               dispersal_stochactic_ = T,
               dispersal_reaction_time_ = 2,
               tau_ = tau,
               BETA_ =  beta_context,  
               TIME_AT_SEA_NB_ = time_at_sea_NB_context,
               reaching.repro.prob_ = reaching.repro.prob)
  
  dt = 
    data.frame(
      scenario = c(
        rep("Healty site", nrow(no_stress)),
        rep("Baseline outbreak",nrow(baseline_outbreak)),
        rep("Proactive strategy",nrow(proactive_strategy)),
        rep("Proactive strategy - Too late",nrow(proactive_strategy_toolate)),
        rep("Reactive strategy",nrow(reactive_strategy))
      ),
      
      equi.lost.survi.ad = c(
        0,
        no_stress$nb_adults_equi - baseline_outbreak$nb_adults_equi,
        no_stress$nb_adults_equi - proactive_strategy$nb_adults_equi,
        no_stress$nb_adults_equi - proactive_strategy_toolate$nb_adults_equi,
        no_stress$nb_adults_equi - reactive_strategy$nb_adults_equi 
      ),
      
      nb_infected_colonies = c(
        no_stress$nb_infected_colonies,
        baseline_outbreak$nb_infected_colonies,
        proactive_strategy$nb_infected_colonies,
        proactive_strategy_toolate$nb_infected_colonies,
        reactive_strategy$nb_infected_colonies
      ),
      
      infected_X_time = c(
        no_stress$infected_X_time,
        baseline_outbreak$infected_X_time,
        proactive_strategy$infected_X_time,
        proactive_strategy_toolate$infected_X_time,
        reactive_strategy$infected_X_time
        
      )
      
      
    ) %>% 
    mutate(scenario = factor(scenario, levels = c("Healty site", 
                                                  "Baseline outbreak",
                                                  "Proactive strategy",
                                                  "Proactive strategy - Too late",
                                                  "Reactive strategy")))
  return(dt)
  
}

scenario_mean = function(dt){
  
  equi.lost.survi.ad_mean = c(
    no_stress = mean(dt[dt$scenario == "Healty site", c("equi.lost.survi.ad")]),
    baseline_outbreak = mean(dt[dt$scenario == "Baseline outbreak", c("equi.lost.survi.ad")]),
    proactive_strategy = mean(dt[dt$scenario == "Proactive strategy", c("equi.lost.survi.ad")]),
    proactive_strategy_toolate = mean(dt[dt$scenario == "Proactive strategy - Too late", c("equi.lost.survi.ad")]),
    reactive_strategy = mean(dt[dt$scenario == "Reactive strategy", c("equi.lost.survi.ad")])
  )
  infected_X_time_mean = c(
    no_stress = mean(dt[dt$scenario == "Healty site", c("infected_X_time")]),
    baseline_outbreak = mean(dt[dt$scenario == "Baseline outbreak", c("infected_X_time")]),
    proactive_strategy = mean(dt[dt$scenario == "Proactive strategy", c("infected_X_time")]),
    proactive_strategy_toolate = mean(dt[dt$scenario == "Proactive strategy - Too late", c("infected_X_time")]),
    reactive_strategy = mean(dt[dt$scenario == "Reactive strategy", c("infected_X_time")])
  )
  nb_infected_colonies_mean = c(
    no_stress = mean(dt[dt$scenario == "Healty site", c("nb_infected_colonies")]),
    baseline_outbreak = mean(dt[dt$scenario == "Baseline outbreak", c("nb_infected_colonies")]),
    proactive_strategy = mean(dt[dt$scenario == "Proactive strategy", c("nb_infected_colonies")]),
    proactive_strategy_toolate = mean(dt[dt$scenario == "Proactive strategy - Too late", c("nb_infected_colonies")]),
    reactive_strategy = mean(dt[dt$scenario == "Reactive strategy", c("nb_infected_colonies")])
  )
  
  dt_mean = data.frame(equi.lost.survi.ad_mean,
                       infected_X_time_mean,
                       nb_infected_colonies_mean)
  return(dt_mean)
}



scenario_plot3 = function(dt,
                          BETA,
                          MOB,
                          REPRO,
                          column_to_plot,
                          index_label,
                          its_dotsize = 0.5){  # Ajout de l'argument column_to_plot
  
  dt_panel_5 = dt %>% 
    mutate(scenario_short = recode(scenario,
                                   "Healty site" = "HS", 
                                   "Baseline outbreak" = "BO", 
                                   "Proactive strategy" = "PS",
                                   "Proactive strategy - Too late" = "P2", 
                                   "Reactive strategy" = "RS"))
  
  sc_mean = scenario_mean(dt)
  
  # Utilisation de !!sym(column_to_plot) pour rendre la colonne dynamique
  p_equi.survi.ad = ggplot()+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site")), 
                 aes(x = scenario_short, y = !!sym(column_to_plot)),
                 binaxis='y', stackdir='center',
                 dotsize = 4,
                 fill = "grey",
                 color = "white",
                 stackratio=0.025)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Baseline outbreak")), 
                aes(x = scenario_short, y = !!sym(column_to_plot)),
                fill = "grey",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy")), 
                 aes(x = scenario_short, y = !!sym(column_to_plot)),
                 binaxis='y', stackdir='center',
                 dotsize = 4,
                 fill = if (sc_mean["baseline_outbreak",paste0(column_to_plot, "_mean")]
                            >sc_mean["proactive_strategy",paste0(column_to_plot, "_mean")]) "lightgreen"
                 else "darksalmon",
                 color = "white",
                 stackratio=0.025)+
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy - Too late")), 
                aes(x = scenario_short, y = !!sym(column_to_plot)),
                fill = if (sc_mean["baseline_outbreak",paste0(column_to_plot, "_mean")]
                           >sc_mean["proactive_strategy_toolate",paste0(column_to_plot, "_mean")]) "lightgreen"
                else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_violin(data = dt_panel_5 %>% subset(., scenario %in% c("Reactive strategy")), 
                aes(x = scenario_short, y = !!sym(column_to_plot)),
                fill = if (sc_mean["baseline_outbreak",paste0(column_to_plot, "_mean")]
                           >sc_mean["reactive_strategy",paste0(column_to_plot, "_mean")]) "lightgreen"
                else "darksalmon",
                color = "white",
                scale = "width",
                trim=FALSE, position=position_dodge(1)) +
    geom_dotplot(data = dt_panel_5 %>% subset(., scenario %in% c("Proactive strategy - Too late")), 
                 aes(x = scenario_short, y = !!sym(column_to_plot)),
                 binaxis='y', stackdir='center',
                 dotsize = its_dotsize,
                 color = "darkgrey", alpha = 0.5,
                 stackratio=0.125)+
    geom_dotplot(data = dt_panel_5 %>% subset(., !(scenario %in% c("Proactive strategy - Too late"))), 
                 aes(x = scenario_short, y = !!sym(column_to_plot)),
                 binaxis='y', stackdir='center',
                 dotsize = its_dotsize,
                 color = "darkgrey", alpha = 0.5,
                 stackratio=0.345)+
    geom_dotplot(data =  dt_panel_5 %>% subset(., scenario %in% c("Healty site", "Proactive strategy")), 
                 aes(x = scenario_short, y = !!sym(column_to_plot)),
                 binaxis='y', stackdir='center',
                 dotsize = its_dotsize,
                 fill = "antiquewhite4",
                 color = "antiquewhite4",
                 stackratio=0.035)+
    ggthemes::theme_clean() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_blank(),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.line = element_line(linewidth = 2),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.position =  "none"
    )+
    ylim(0, NA)+
    labs(x = paste0("beta = ", as.character(BETA),
                    "  mob = ", as.character(MOB),
                    "  repro = ", as.character(REPRO)),
         y = index_label)  # Mettre à jour l'étiquette Y avec le nom de la colonne
  print(p_equi.survi.ad)
}








beta_bank = c(0.02, 0.05, 0.5, 2)
mob_bank = c(2, 8, 40)
repro_bank = c(0.4) # c(0.3, 0.5, 0.7)


plot_ENLSA <- list()
plot_infec.col <- list()
plot_infected_X_time <- list()

i = 0

for (beta in beta_bank){
  for (mob in mob_bank){
    for (repro in repro_bank){
      
      i = i + 1
      
      dt_ = scenario_dt(beta_context = beta,
                        time_at_sea_NB_context = mob,
                        reaching.repro.prob = repro)
      
      plot_ENLSA[[i]] <- scenario_plot3(dt = dt_, BETA = beta, MOB=mob, REPRO=repro,
                                        column_to_plot = "equi.lost.survi.ad",
                                        index_label = "ENLSA",
                                        its_dotsize = 0.6)
      
      plot_infec.col[[i]] <- scenario_plot3(dt = dt_, BETA = beta, MOB=mob, REPRO=repro,
                                            column_to_plot = "nb_infected_colonies",
                                            index_label = "Infected colonies",
                                            its_dotsize = 1.1)
      
      plot_infected_X_time[[i]] <- scenario_plot3(dt = dt_, BETA = beta, MOB=mob, REPRO=repro,
                                                  column_to_plot = "infected_X_time",
                                                  index_label = "Infected X Time colonies",
                                                  its_dotsize = 1.1)
      
      
    }
  }
}

all_plot_ENLSA = grid.arrange(grobs = plot_ENLSA, ncol = 3, nrow = 4)
all_plot_infec.col = grid.arrange(grobs = plot_infec.col, ncol = 3, nrow = 4)
all_plot_infected_X_time = grid.arrange(grobs = plot_infected_X_time, ncol = 3, nrow = 4)


