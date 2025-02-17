library(shiny)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(abind)
library(cowplot)
library(lhs)
library(dplyr)
library(stringr)
library(rlang)
library(gridExtra)
library(grid)
library(shinythemes)
library(shinydashboard)
library(shinyWidgets)



#setwd("C:/Data/1_Adminitratif/Emploi/cornell/02_projet/hpai_mitigation/02_models")


#source("param_ranges.R")
# Parameter ranges --------------------------------------------------------

param_ranges = list(
  
  # Integer parameters
  initial_number_infected_breeders_A = list(range = c(0, 3), distribution = "integer"),
  initial_number_infected_breeders_B = list(range = c(0, 0), distribution = "integer"),
  initial_number_infected_breeders_C = list(range = c(0, 0), distribution = "integer"),
  initial_number_breeders_A = list(range = c(50, 50), distribution = "integer"),
  initial_number_breeders_B = list(range = c(80, 80), distribution = "integer"),
  initial_number_breeders_C = list(range = c(20, 20), distribution = "integer"),
  dispersal_reaction_time =  list(range = c(1, 6), distribution = "integer"),
  dispersal_date = list(range = c(0, 0), distribution = "integer"),
  hatching_date = list(range = c(0, 20), distribution = "integer"),
  
  # Continuous parameters
  tau = list(range = c(0.15, 0.15), distribution = "simple_uniform"),
  total_time = list(range = c(70, 70), distribution = "simple_uniform"),
  prop_dispersal = list(range = c(0.9, 1), distribution = "simple_uniform"),
  beta_E_colony = list(range = c(0, 0), distribution = "simple_uniform"),
  beta_I_colony = list(range = c(0.01, 0.80), distribution = "logarithmic"),
  incubation_period = list(range = c(0.8, 1.2), distribution = "simple_uniform"),
  eta = list(range = c(0, 0), distribution = "simple_uniform"),
  infectious_period = list(range = c(5, 7), distribution = "logarithmic"),
  adult_mortality = list(range = c(0.4,0.8), distribution = "simple_uniform"),
  nestling_mortality = list(range = c(0.6, 1.0), distribution = "simple_uniform"),
  avrg_stay_B_colony = list(range = c(2, 2), distribution = "logarithmic"),
  avrg_stay_B_sea = list(range = c(2, 2), distribution = "logarithmic"),
  avrg_stay_NB_colony = list(range = c(2, 2), distribution = "logarithmic"),
  avrg_stay_NB_sea = list(range = c(2, 40), distribution = "logarithmic"),
  theta = list(range = c(1/100, 1/7), distribution = "simple_uniform"),
  psi = list(range = c(1/500, 1/500), distribution = "logarithmic"),
  hatching_sd = list(range = c(3, 3), distribution = "simple_uniform"),
  reaching_repro_prob = list(range = c(0.3, 0.7), distribution = "simple_uniform"),
  prob_detection = list(range = c(0.7, 1.0), distribution = "simple_uniform")
)
#source("scenarios.R")
# Describe the 5 scenario
scenarios = data.frame(
  induced_dispersal = c(F,F,T,T,T),
  dispersal_stochastic = c(F,F,F,F,T),
  initially_infected = c(F,T,F,T,T))

rownames(scenarios) = c("HS","BO","PS","P2","RS")
scenarios = scenarios[c(2,4,5),]




load("simulation_dt/dt_cluster_3_1E4_for_v2.RData")

param = c("initial_number_infected_breeders_A", "initial_number_infected_breeders_B", "initial_number_infected_breeders_C", "initial_number_breeders_A",         
          "initial_number_breeders_B",          "initial_number_breeders_C",          "dispersal_reaction_time",            "dispersal_date",                    
          "hatching_date",                      "tau",                                "total_time",                         "prop_dispersal",                    
          "beta_E_colony",                      "beta_I_colony",                      "incubation_period",                  "eta",                               
          "infectious_period",                  "adult_mortality",                    "nestling_mortality",                 "avrg_stay_B_colony",                
          "avrg_stay_B_sea",                    "avrg_stay_NB_colony",                "avrg_stay_NB_sea",                   "theta",                             
          "psi",                                "hatching_sd",                        "reaching_repro_prob",                "prob_detection"  )

param_evaluated = c("initial_number_infected_breeders_A",        
                    "hatching_date",                      "prop_dispersal",                    
                    "beta_I_colony",                                                                 
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


ui <- fluidPage(
  theme = shinytheme("cyborg"), 
  
  titlePanel("Interactive Heatmap"),
  
  sidebarLayout(
    sidebarPanel(
      chooseSliderSkin("Modern"),
      sliderInput("n_bins", "Number of Bins:", 
                  min = 4, 
                  max = 20, 
                  value = 10, 
                  step = 1),
      
      selectInput("selected_output", "Select Output:", 
                  choices = c(#"BO_nb_adults_equi", 
                              "BO_RS_nb_adults_equi", "BO_P2_nb_adults_equi"),
                  selected = "BO_P2_nb_adults_equi"),
      
      checkboxGroupInput("selected_params", "Select the parameters :", choices = param_evaluated, 
                         selected = c("beta_I_colony", "initial_number_infected_breeders_A", "adult_mortality")),
      
      sliderInput("initial_number_infected_breeders_A_range", 
                  "Number of infected breeders introduced initially", 
                  min = min(dt2$initial_number_infected_breeders_A), 
                  max = max(dt2$initial_number_infected_breeders_A), 
                  value = c(min(dt2$initial_number_infected_breeders_A), max(dt2$initial_number_infected_breeders_A)), 
                  step = 1),
      
      sliderInput("beta_I_colony_range", 
                  "Transmission rate", 
                  min = 0.02, 
                  max = 0.75, 
                  value = c(min(dt2$beta_I_colony), max(dt2$beta_I_colony)), 
                  step = diff(range(dt2$beta_I_colony))/20),
      
      sliderInput("hatching_date",
                  "Hatching date",
                  min = min(dt2$hatching_date),
                  max = max(dt2$hatching_date),
                  value = c(min(dt2$hatching_date), max(dt2$hatching_date)),
                  step = diff(range(dt2$hatching_date))/20),
      
      sliderInput("prop_dispersal",
                  "Proportion of breeders dispersed",
                  min = 0.9,
                  max = 0.999,
                  value = c(min(dt2$prop_dispersal), max(dt2$prop_dispersal)),
                  step = diff(range(dt2$prop_dispersal))/20),
      
      sliderInput("infectious_period",
                  "Infectious period",
                  min = min(dt2$infectious_period),
                  max = max(dt2$infectious_period),
                  value = c(min(dt2$infectious_period), max(dt2$infectious_period)),
                  step = diff(range(dt2$infectious_period))/20),
      
      sliderInput("adult_mortality",
                  "Adult mortality",
                  min = min(dt2$adult_mortality),
                  max = max(dt2$adult_mortality),
                  value = c(min(dt2$adult_mortality), max(dt2$adult_mortality)),
                  step = diff(range(dt2$adult_mortality))/20),
      
      sliderInput("nestling_mortality",
                  "Nestling mortality",
                  min = min(dt2$nestling_mortality),
                  max = max(dt2$nestling_mortality),
                  value = c(min(dt2$nestling_mortality), max(dt2$nestling_mortality)),
                  step = diff(range(dt2$nestling_mortality))/20),
      
      sliderInput("avrg_stay_NB_sea",
                  "avrg_stay_NB_sea period",
                  min = min(dt2$avrg_stay_NB_sea),
                  max = max(dt2$avrg_stay_NB_sea),
                  value = c(min(dt2$avrg_stay_NB_sea), max(dt2$avrg_stay_NB_sea)),
                  step = diff(range(dt2$avrg_stay_NB_sea))/20),
      
      sliderInput("theta",
                  "theta mortality",
                  min = min(dt2$theta),
                  max = max(dt2$theta),
                  value = c(min(dt2$theta), max(dt2$theta)),
                  step = diff(range(dt2$theta))/20),
      
      sliderInput("reaching_repro_prob",
                  "reaching_repro_prob mortality",
                  min = min(dt2$reaching_repro_prob),
                  max = max(dt2$reaching_repro_prob),
                  value = c(min(dt2$reaching_repro_prob), max(dt2$reaching_repro_prob)),
                  step = diff(range(dt2$reaching_repro_prob))/20),
      
      sliderInput("prob_detection",
                  "prob_detection mortality",
                  min = min(dt2$prob_detection),
                  max = max(dt2$prob_detection),
                  value = c(min(dt2$prob_detection), max(dt2$prob_detection)),
                  step = diff(range(dt2$prob_detection))/20),
      

      
      br(),
      helpText("Adjust the sliders to visualize changes in the heatmap.", 
               style = "font-size: 12px; text-align: center; color: grey;")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", 
                 plotOutput("heatmapPlot", width = "700px", height = "600px")),
        tabPanel("Regression tree",
                 verbatimTextOutput("summary"))
      )
    )
  )
)



server <- function(input, output) {
  output$heatmapPlot <- renderPlot({
    
    
    SELECTED_OUTPUT = input$selected_output
    
    N_BINS = input$n_bins
    
    SCENARIO = ""
    
    
    dt2 = dt2 %>% 
      subset(., beta_I_colony >= input$beta_I_colony_range[1])%>% 
      subset(., beta_I_colony <= input$beta_I_colony_range[2]) %>% 
      subset(., initial_number_infected_breeders_A >= input$initial_number_infected_breeders_A_range[1]) %>% 
      subset(., initial_number_infected_breeders_A <= input$initial_number_infected_breeders_A_range[2])%>%
      subset(., hatching_date >= input$hatching_date[1]) %>%
      subset(., hatching_date <= input$hatching_date[2])%>%
      subset(., prop_dispersal >= input$prop_dispersal[1]) %>%
      subset(., prop_dispersal <= input$prop_dispersal[2])%>%
      subset(., infectious_period >= input$infectious_period[1]) %>%
      subset(., infectious_period <= input$infectious_period[2])%>%
      subset(., adult_mortality >= input$adult_mortality[1]) %>%
      subset(., adult_mortality <= input$adult_mortality[2])%>%
      subset(., nestling_mortality >= input$nestling_mortality[1]) %>%
      subset(., nestling_mortality <= input$nestling_mortality[2])%>%
      subset(., avrg_stay_NB_sea >= input$avrg_stay_NB_sea[1]) %>%
      subset(., avrg_stay_NB_sea <= input$avrg_stay_NB_sea[2])%>%
      subset(., reaching_repro_prob >= input$reaching_repro_prob[1]) %>%
      subset(., reaching_repro_prob <= input$reaching_repro_prob[2])%>%
      subset(., prob_detection >= input$prob_detection[1]) %>%
      subset(., prob_detection <= input$prob_detection[2])%>%
      subset(., theta >= input$theta[1]) %>%
      subset(., theta <= input$theta[2])

    
    # c("initial_number_infected_breeders_A",        
    #   "hatching_date",                      "prop_dispersal",                    
    #   "beta_I_colony",                                                                 
    #   "infectious_period",                  "adult_mortality",   
    #   "nestling_mortality",                         
    #   "avrg_stay_NB_sea",                   "theta",                             
    #   "reaching_repro_prob",                "prob_detection"  )
    # 
    
    evaluated_parameter = input$selected_params
    
    graph_param_name = input$selected_params
    
    dt2 = dt2 %>% subset(., initial_number_infected_breeders_A!=0)
    
    
    
    # -------------------------------------------------------------------------
    
    simulation_dt = dt2
    
    
    # Function to create binned data with dynamic parameters and variable block sizes
    create_binned_data = function(data,
                                  params,
                                  param_ranges_,
                                  selected_output,
                                  n_bins) {
      
      # Apply log transformation if specified
      if (param_ranges_[[params[1]]][[2]] == "logarithmic") {
        data$var1 = log(data[[params[1]]])
      } else {
        data$var1 = data[[params[1]]]
      }
      
      if (param_ranges_[[params[2]]][[2]] == "logarithmic") {
        data$var2 = log(data[[params[2]]])
      } else {
        data$var2 = data[[params[2]]]
      }
      
      # Calculate min and max for both parameters
      min_x = min(data$var1, na.rm = TRUE)
      max_x = max(data$var1, na.rm = TRUE)
      min_y = min(data$var2, na.rm = TRUE)
      max_y = max(data$var2, na.rm = TRUE)
      
      # Define dynamic block sizes
      block_size_x = (max_x - min_x) / n_bins
      block_size_y = (max_y - min_y) / n_bins
      
      # Convert selected_output to a symbol
      selected_output_sym = sym(selected_output)
      
      res = data %>%
        mutate(
          V1 = cut(data$var1, breaks = seq(min_x, max_x, by = block_size_x), include.lowest = TRUE),
          V2 = cut(data$var2, breaks = seq(min_y, max_y, by = block_size_y), include.lowest = TRUE)
        ) %>%
        group_by(V1, V2) %>%
        summarise(output_avg = mean(!!selected_output_sym, na.rm = TRUE), .groups = 'drop') %>%
        ungroup() %>%
        mutate(
          # Use str_extract to extract the first numeric value from V1 and V2
          x_mid = as.numeric(str_extract(V1, "[-+]?[0-9]*\\.?[0-9]+")) + block_size_x / 2,
          y_mid = as.numeric(str_extract(V2, "[-+]?[0-9]*\\.?[0-9]+")) + block_size_y / 2
        )
      
      # Apply exp transformation if specified
      if (param_ranges_[[params[1]]][[2]] == "logarithmic") {
        res$x_mid = exp(res$x_mid)
      } else {
        res$x_mid = res$x_mid
      }
      
      if (param_ranges_[[params[2]]][[2]] == "logarithmic") {
        res$y_mid = exp(res$y_mid)
      } else {
        res$y_mid = res$y_mid
      }
      
      res = res %>% 
        dplyr::select(x_mid, y_mid, output_avg)
      
      return(res)
    }
    
    
    
    # Function to create a heatmap with dynamic block sizes
    plot_heatmap_binned_diff = function(data, params, param_ranges) {
      
      
      p = ggplot(data) +
        geom_tile(aes(x = x_mid, y = y_mid, fill = output_avg)
        ) +
        scale_fill_gradient2(low = "brown3", mid = "white", high = "chartreuse4",
                             midpoint = 0,
                             limits = c(min(data$output_avg), max(data$output_avg)))+
        labs(
          #x = params[1], y = params[2],
          x = NULL, y = NULL,
          fill = "ENLA")+
        guides(fill = "none") +
        theme_minimal()
      
      if (param_ranges[[params[1]]][[2]] == "logarithmic"){
        p = p + scale_x_log10()
      }
      if (param_ranges[[params[2]]][[2]] == "logarithmic"){
        p = p + scale_y_log10()
      }
      
      return(p)
    }
    
    
    # Empty list to store results
    all_diff_results <- list()
    
    # Loop through scenarios and calculate differences
    for (index_param_1 in 1:(length(evaluated_parameter))) {
      for (index_param_2 in 1:(length(evaluated_parameter))) {  
        
        
        
        diff_result = create_binned_data(data = simulation_dt,
                                         params = c(evaluated_parameter[index_param_1],
                                                    evaluated_parameter[index_param_2]),
                                         param_ranges_ = param_ranges ,
                                         selected_output = SELECTED_OUTPUT,
                                         n_bins = N_BINS)
        
        all_diff_results[[paste(evaluated_parameter[index_param_1], evaluated_parameter[index_param_2], sep = "__")]] <- list(value = diff_result,
                                                                                                                              param = c(evaluated_parameter[index_param_1],
                                                                                                                                        evaluated_parameter[index_param_2])                                                                                                                          )
      }
    }
    
    
    # Create a list of plots with lapply
    plots <- lapply(1:length(all_diff_results), function(i) {
      plot_heatmap_binned_diff(all_diff_results[[i]]$value,
                               all_diff_results[[i]]$param,
                               param_ranges)
    })
    
    # List of column labels
    col_labels_ <- graph_param_name
    
    # List of row labels
    row_labels_ <- col_labels_
    
    # Create text grobs for column labels (rotated 45 degrees)
    col_labels <- lapply(1:(length(evaluated_parameter)), function(i) {
      textGrob(label = col_labels_[i], rot = 45, gp = gpar(fontsize = 10, fontface = "bold"))
    })
    
    # Create text grobs for row labels (rotated 45 degrees for the left side)
    row_labels <- lapply(1:(length(evaluated_parameter)), function(i) {
      textGrob(label = row_labels_[i], rot = 45, gp = gpar(fontsize = 10, fontface = "bold"))
    })
    
    # Create layout matrix
    n <- length(evaluated_parameter)
    
    # Create an empty matrix for layout with n+1 rows and columns for the labels and plots
    layout_matrix <- matrix(0, nrow = n+1, ncol = n+1)
    
    # Fill the layout matrix
    layout_matrix[1, 2:(n+1)] <- 1:n   # Top row for column labels
    layout_matrix[2:(n+1), 1] <- (n+1):(2*n)  # First column for row labels
    layout_matrix[2:(n+1), 2:(n+1)] <- (2*n+1):(n*n+2*n)  # The rest for the plots
    
    # Combine nullGrob, row labels, col labels, and plots into a single list of grobs
    all_grobs <- c(list(nullGrob()), col_labels, row_labels, plots)
    
    
    plot_one_param = function(df, evaluated_param){
      
      
      # Dynamic plotting based on the evaluated parameter
      p = ggplot() +
        geom_point(data=df,
                   aes_string(x = evaluated_param, y = SELECTED_OUTPUT),  # Couleur dans aes()
                   size = 0.4, alpha = 0.14) +  # Transparence appliquée
        #geom_smooth(method = "loess") +
        theme_minimal() +
        labs(x = NULL, y = NULL)
      
      return(p)
    }
    
    
    plot_one_param_bank = list()
    
    for (i in 1:length(evaluated_parameter)){
      
      plot_one_param_bank[[i]] = plot_one_param(df = simulation_dt,
                                                evaluated_param = evaluated_parameter[i])
    }
    
    
    
    #   -----------------------------------------------------------------------
    
    nb_param = length(evaluated_parameter)
    
    for (i in 1:length(evaluated_parameter)){
      all_grobs[[(nb_param+1)+((nb_param+1)*i)]] = plot_one_param_bank[[i]]
    }
    
    
    p = grid.arrange(
      grobs = all_grobs,  
      layout_matrix = layout_matrix, 
      top = paste0(" ", SCENARIO)
    )
    
    plot(p)
    
  })
}

shinyApp(ui = ui, server = server)

