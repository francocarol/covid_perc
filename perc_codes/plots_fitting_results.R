################################################
################################################
################################################

# PLOT FITTING RESULTS

################################################
################################################
################################################

# New cases and deaths for all versions:

plot_srag_allv <- function(parameters, interventions, Y, params = NA,
                           base_size=18,
                           rowmean_line_color = "grey",rowmean_line_size = 1.2,
                           weekly = TRUE,
                           output = "cases",
                           plot_object, colour = "black"){
  for (i in names(params)){parameters[i] <- params[i]}
  
  out <- run_covid19(parameters, interventions, Y)
  res <- result_ts(out[[2]], parameters, parameters["startdate"], diff=F)
  data_c_fit = diff(res$C)
  data_m_fit = diff(res$CM)
  
  if(weekly) {
    model_c_xts <- xts::apply.weekly(data_c_fit, sum)
    model_cm_xts <- xts::apply.weekly(data_m_fit, sum)
    data_c_fit <- zoo(model_c_xts, as.Date(index(model_c_xts)))
    data_m_fit <- zoo(model_cm_xts, as.Date(index(model_cm_xts)))
  }
  
  if (output=="cases") {
    plot_object +
      geom_line(data_c_fit, mapping = aes(x=time(data_c_fit), y=data_c_fit, colour=colour))# +
    #xlim(as.Date('2020-3-1'), xmax_c) +
    #ylim(c(0, max(data_c, window(data_c_fit, end=xmax_c))))
  }
  else if (output=="mort") {
    plot_object +
      geom_line(data_m_fit, mapping = aes(x=time(data_m_fit), y=data_m_fit, colour=colour))# +
    #xlim(as.Date('2020-3-1'), xmax_m) +
    #ylim(c(0, max(data_m, window(data_m_fit, end=xmax_m))))
  }
}

plot_together <- function(file1, file2, file3, Scn1, Scn2, Scn3,
                          y_lim_cases, y_lim_mort){
  source("model2_least_squares_fitting_homesteep_v2.R")
  #color_seq <- colorRampPalette(brewer.pal(8, "Accent"))(3)
  cols <- c("Standard model"="red", "Standard + 30% NPIs"="blue", "Model with percolation"="black")
  base_size=18
  plot.formatos <- theme_bw(base_size=base_size)
  data_c = diff(now.srag.zoo)
  data_m = diff(now.obito.srag.zoo)
  max.date = as.Date('2020-08-31')
  
  pcases <- ggplot(data_c) + geom_point(data=data_c, aes(x=time(data_c), y=data_c)) +
    plot.formatos +
    xlab("Time") + ylab("New Cases") +
    xlim(as.Date('2020-3-1'), max.date) +
    ylim(c(0, y_lim_cases)) #max(data_c, window(1000, end=xmax_c))))
  
  pmort <-  ggplot(data_m) + geom_point(data=data_m, aes(x=time(data_m), y=data_m)) +
    plot.formatos +
    xlab("Time") + ylab("New Deaths") +
    xlim(as.Date('2020-3-1'), max.date) +
    ylim(c(0, y_lim_mort))#max(data_m,window(0, end=xmax_m))))
  
  # STANDART MODEL (version 2)
  base_parameters = c()
  parameters2 <- update_parameters(base_parameters, estado = estado, cidade = cidade)
  parameters2 <- replace(parameters2, c("perc_threshold","home_steep", "home_eff"), c(0,0,0))
  parameters2["startdate"] <- as.numeric(as.Date("2020-02-01")) ####
  interventions <- update_scenario_parameters(scenario=Scn2, parameters2)
  #
  fit.results <- read.csv(file2)
  dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]
  pars_to_update <- names(parameters2)[names(parameters2) %in% names(dopt)]
  parameters2[pars_to_update] <- unlist(dopt[pars_to_update])
  opt <- fit.results[which.max(fit.results$prob),]
  interventions <- update_scenario_parameters(parameters2, scenario = Scn2, interventions)
  
  pcases <- plot_srag_allv(parameters2, interventions, Y, params=unlist(opt[1:4]),
                           weekly = TRUE,
                           output="cases" , plot_object=pcases, colour = "Standard model")
  pmort <- plot_srag_allv(parameters2, interventions, Y, params=unlist(opt[1:4]),
                          weekly = TRUE,
                          output="mort" , plot_object=pmort, colour = "Standard model")
  
  # MANUAL FIT TO STANDART MODEL: 30% more coverage on NPIs (version 3)
  base_parameters = c()
  parameters3 <- update_parameters(base_parameters, estado = estado, cidade = cidade)
  parameters3 <- replace(parameters3, c("perc_threshold","home_steep", "home_eff"), c(0,0,0))
  parameters3["startdate"] <- as.numeric(as.Date("2020-01-01")) #("2020-02-01")
  interventions <- update_scenario_parameters(scenario=Scn3, parameters3)
  #
  fit.results <- read.csv(file3)
  dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]
  pars_to_update <- names(parameters3)[names(parameters3) %in% names(dopt)]
  parameters3[pars_to_update] <- unlist(dopt[pars_to_update])
  opt <- fit.results[which.max(fit.results$prob),]
  interventions <- update_scenario_parameters(parameters3, scenario = Scn3, interventions)
  
  pcases <- plot_srag_allv(parameters3, interventions, Y, params=unlist(opt[1:4]),
                           weekly = TRUE,
                           output="cases" , plot_object=pcases, colour = "Standard + 30% NPIs")
  pmort <- plot_srag_allv(parameters3, interventions, Y, params=unlist(opt[1:4]),
                          weekly = TRUE,
                          output="mort" , plot_object=pmort, colour = "Standard + 30% NPIs")
  
  # PERC MODEL (version 1)
  base_parameters = c()
  parameters1 <- update_parameters(base_parameters, estado = estado, cidade = cidade) # Aqui especifico cidade e estado
  parameters1["startdate"] <- as.numeric(as.Date("2020-01-20")) # A data está iniciando em 20/jan agora, e esse valor precisa bater com o data_start no script de least_squares_fitting
  interventions <- update_scenario_parameters(scenario=Scn1, parameters1)
  #
  fit.results <- read.csv(file1)
  dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]
  pars_to_update <- names(parameters1)[names(parameters1) %in% names(dopt)]
  parameters1[pars_to_update] <- unlist(dopt[pars_to_update])
  opt <- fit.results[which.max(fit.results$prob),]
  interventions <- update_scenario_parameters(parameters1, scenario = Scn1, interventions)
  
  pcases <- plot_srag_allv(parameters1, interventions, Y, params=unlist(opt[1:4]),
                           weekly = TRUE,
                           output="cases" , plot_object=pcases, colour = "Model with percolation")# +
  #scale_color_manual(values = color_seq, name = legend.title, labels = scn.labels)
  pmort <- plot_srag_allv(parameters1, interventions, Y, params=unlist(opt[1:4]),
                          weekly = TRUE,
                          output="mort" , plot_object=pmort, colour = "Model with percolation")# +
  #scale_color_manual(values = color_seq, name = legend.title, labels = scn.labels)
  
  pcases
  pmort
  
  combined_plots <- ( pcases + theme(legend.position = "none" ) ) | 
    (pmort + theme(legend.position = "right", legend.direction="vertical"))
  #print(combined_plots)
}

# # Plot model output and data for Cases and Deaths (New and Cumulative)
# 
# # MODEL V 1
# 
# #Sys.setlocale("LC_TIME", "en_GB.UTF-8")
# 
# fit.results <- read.csv("fits/SP/model1_srag.1.CM.log0-2021_03_22_SP.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_12_02_home_steep_var.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_11_11_homesteep_var.csv")
# dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]
# pars_to_update <- names(parameters1)[names(parameters1) %in% names(dopt)]
# parameters1[pars_to_update] <- unlist(dopt[pars_to_update])
# opt <- fit.results[which.max(fit.results$prob),]
# interventions <- update_scenario_parameters(parameters1, scenario = 50, interventions)
# 
# filename = paste0("fitting_homesteep_sragdata_aug20_model1_eng.png")
# png(filename, height = 480*3, width = 480*4, res = 200)
# test_modifying_parameters(parameters1, interventions, Y, unlist(opt[1:4]), flag='srag',
#                           rollwindow = 1, weekly = TRUE, max.date = as.Date('2020-08-31'))
# dev.off()
# 
# # MODEL V 2
# 
# #Sys.setlocale("LC_TIME", "en_GB.UTF-8")
# 
# fit.results <- read.csv("fits/SP/model2_srag.1.CM.log0-2021_03_22_SP.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_12_02_home_steep_var.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_11_11_homesteep_var.csv")
# dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]
# pars_to_update <- names(parameters2)[names(parameters2) %in% names(dopt)]
# parameters2[pars_to_update] <- unlist(dopt[pars_to_update])
# opt <- fit.results[which.max(fit.results$prob),]
# interventions <- update_scenario_parameters(parameters2, scenario = 50, interventions)
# 
# filename = paste0("fitting_homesteep_sragdata_aug20_model2_eng.png")
# png(filename, height = 480*3, width = 480*4, res = 200)
# test_modifying_parameters(parameters, interventions, Y, unlist(opt[1:4]), flag='srag',
#                           rollwindow = 1, weekly = TRUE, max.date = as.Date('2020-08-31'))
# dev.off()

# Plot model output and data for Cases and Deaths (Only New)

# MODEL V 1

#Sys.setlocale("LC_TIME", "en_GB.UTF-8")
# 
# fit.results <- read.csv("fits/SP/model1_srag.1.CM.log0-2021_03_22_SP.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_12_02_home_steep_var.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_11_11_homesteep_var.csv")
# dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]
# pars_to_update <- names(parameters1)[names(parameters1) %in% names(dopt)]
# parameters1[pars_to_update] <- unlist(dopt[pars_to_update])
# opt <- fit.results[which.max(fit.results$prob),]
# interventions <- update_scenario_parameters(parameters1, scenario = 50, interventions)
# 
# filename = paste0("fitting_homesteep_sragdata_aug20_model1_eng_newC_newM.png")
# png(filename, height = 480*3, width = 480*2, res = 200)
# test_modifying_parameters_new(parameters1, interventions, Y, unlist(opt[1:4]), flag='srag',
#                           rollwindow = 1, weekly = TRUE, max.date = as.Date('2020-08-31'))
# dev.off()
# 
# # MODEL V 2
# 
# #Sys.setlocale("LC_TIME", "en_GB.UTF-8")
# 
# fit.results <- read.csv("fits/SP/model2_srag.1.CM.log0-2021_03_22_SP.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_12_02_home_steep_var.csv")
# #fit.results <- read.csv("fits/srag.1.CM.log0-2020_11_11_homesteep_var.csv")
# dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]
# pars_to_update <- names(parameters2)[names(parameters2) %in% names(dopt)]
# parameters2[pars_to_update] <- unlist(dopt[pars_to_update])
# opt <- fit.results[which.max(fit.results$prob),]
# interventions <- update_scenario_parameters(parameters2, scenario = 50, interventions)
# 
# filename = paste0("fitting_homesteep_sragdata_aug20_model2_eng_newC_newM.png")
# png(filename, height = 480*3, width = 480*2, res = 200)
# test_modifying_parameters_new(parameters, interventions, Y, unlist(opt[1:4]), flag='srag',
#                           rollwindow = 1, weekly = TRUE, max.date = as.Date('2020-08-31'))
# dev.off()
