if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(zeallot)){install.packages("zeallot"); library(zeallot)} # provides %<-%
if(!exists("sheets_read",mode="function")) library(googlesheets4)


setwd("/home/caroline/Documents/paper_covid/codes/")

source('parameter_calibration_functions.R')
#source('calculate_R_eff.R')
source('plots.R')
source('scenarios_names_v2.R')

### Definir cidade:
#estado = "GO"; cidade = "Goiania"
estado =  "SP"; cidade = "Sao_Paulo"
#estado =  "RS"; cidade = "Porto_Alegre"

### Carregar modelo e parâmetros de acordo com a cidade
source("Brazilian_Model_v2b.R")

### Carregar dados de SRAG e COVID da base SIVEP
source("dados_municipios.R")

# Load srag and covid data from specified city
city.zoo.list <- carregar_dados_cidades(estado = estado, cidade = cidade)
for(i in 1:length(city.zoo.list)) assign(names(city.zoo.list)[i], city.zoo.list[[i]])

#lang <- 'pt'
lang <- 'en'
    
# function to run a simulation for a set of scenarios
run_scn <- function(scenarios, base_parameters, calculateR=FALSE, estado = "SP", cidade = "Sao_Paulo", model_version="model1"){
  res <- list()
  reffs <- list()
  outs <- list()
  base_parameters <- update_parameters(base_parameters)
  # use latest fit
  base_parameters <- update_params_from_fit(base_parameters,
                                            paste0('fits/', get_latest_fit(estado = estado, model_version = model_version)[1]))
  
  base_interventions <- update_scenario_parameters(base_parameters, estado = estado, cidade = cidade)
  
  for (scn in 1:length(scenarios)){
    parameters <- base_parameters
    if(scenarios[scn] != 34) { # 34 is the base scenario
      interventions <- reload_scenario_parameters(parameters, scenario = scenarios[scn], base_interventions)
    } else {interventions <- base_interventions}
    
    out <- run_covid19(parameters, interventions, Y)
    list_cols <- list(I=c(Iindex,CLindex,Xindex,Hindex,HCindex,
                          ICUindex,ICUHindex,ICUCindex,QIindex,QCindex)+1,
                      E=c(Eindex,QEindex)+1,
                      R=c(Rindex, QRindex)+1,
                      C=Cindex+1,
                      H=c(Hindex, HCindex)+1,
                      U=c(ICUindex, ICUHindex, ICUCindex)+1,
                      M=CMindex+1,
                      Q=c(QSindex,QEindex,QIindex,QRindex,QCindex,HCindex,ICUCindex,Xindex)+1)
    result <- result_ts(out[[2]], parameters,parameters["startdate"], diff=F,
                        cols=list_cols)
    outs[[scn]] <- out[[2]]
    #result$I <- result$I/sum(Y)
    result$R <- result$R/sum(Y)
    result$MD <- diff(result$M) # daily mortality
    res[[scn]] <- result
    if(calculateR){
      como <- calcR(result,method="CoMo")
      cori <- calcR(result,method="Cori")
      corinew <- calcR(result,method="Cori_new")
      reffs[[scn]] <- list(Reff_como = como, Reff_cori = cori, Reff_cori_sample = corinew)}
  }
  
  # put all results in a single data.frame
  results <- do.call(merge.zoo, res)
  colnames(results) <- c(paste0(c(names(list_cols), "MD"), 
                                rep(1:length(scenarios), each=length(list_cols)+1)))
  
  final.results <- list(results = results, reffs = reffs, outs = outs)
  
  return(final.results)
}

#scenarios <- c(11,14:17) #c(3,7,11:13) #c(2,3,5,7,8,9)
#scenarios <- c(20) #, 11, 19, 9, 29) # partial school reopening (October - December)
scenarios <- c(51,52,54)
parameters<-c()
result1 <- run_scn(scenarios, parameters, estado = "SP", cidade = "Sao_Paulo", model_version="model1")
result2 <- run_scn(scenarios, parameters, estado = "SP", cidade = "Sao_Paulo", model_version="model2")
#result <- run_scn(scenarios, base_parameters)

##se quiser salvar em .csv
# df <- data.frame(data=time(result))
# df$C1 <- result$C1  %>% round()
# df$MD1 <- result$MD1 %>% round()
# write.csv(df,"cenario_keep_curr.csv", row.names = FALSE)

if (lang == 'pt'){
  scn.labels <- names(scn.ALL.PT)[ match(scenarios,scn.ALL.PT) ]
  ylabels=c("Número de Infectados",
            "Número de Expostos (em incubação)",
            "Fração de recuperados",
            "Casos observ. acum.",
            "Hospitalizados",
            "UTI", "Óbitos acumulados", "Óbitos diários",'Quarentenados',
            "R Efetivo - Cori",
            "R Efetivo - CoMo")
}
if (lang == 'en'){
  scn.labels <- names(scn.ALL)[ match(scenarios,scn.ALL) ]
  ylabels=c("Infected",
            "Exposed",
            "Fraction recovered",
            "C",#"Cumulative rep. cases",#"Total observed cases",
            "New cases",#"Daily cases", #"Hospitalized",
            "ICU",
            "CM",#"Cumulative rep. deaths",#"Total deaths"
            "New deaths",#"Daily deaths",#"Total deaths",
            'Quarantined',
            "Effective R - Cori",
            "Effective R - CoMo")
}

# create all plots

# MODEL V1
options(scipen=10000)
plots <- plot_scenarios_simple_zoo(window(result1$result, start=as.Date('2020/01/20'), end=as.Date('2020/12/01'),format="%Y/%m/%d"),
                                   #plots <- plot_scenarios_simple_zoo(result2$result,
                                   scn.labels, data='srag',
                                   cols=c('I', 'E', 'R', 'C', 'H', 'U', 'M', 'MD','Q'),#c('I', 'r', 'C', 'H', 'U', 'MD',"Reff_cori","Reff_como"),
                                   ylabels=ylabels, lang=lang, end.data=as.Date('2020/08/31'),
                                   hlines=c(r=1-1/2.5,
                                            U=unname(parameters['icu_beds_available']),
                                            H=unname(parameters['beds_available'])),
                                   label_months=T, line.size=2,color.palette="Dark2")
# or plot only a certain window 
# plots <- plot_scenarios_simple_zoo(window(result, start=as.Date('2020-3-1'), end=as.Date('2020-8-1')), scn.labels)

# put the plots in a grid
#six_plots <- plots[c(1,4:7,9)] #select which 6 variables to plot
#six_plots <- plots[c(1,4,5,6,7,9)] #select which 6 variables to plot
#combined_plots <- combine.plots(six_plots, legend.pos='bottom')#, margin=0.2+0.1*length(scenarios))
four_plots <- plots[c(5,4,7,8)] #select which 6 variables to plot
combined_plots <- (( plots[[4]] + theme(legend.position = "bottom", legend.direction="vertical") | plots[[7]] + theme(legend.position = "none")                                               ) / ( plots[[5]] + theme(legend.position = "none") | plots[[8]]  + theme(legend.position = "none")                                               ) / guide_area()) + plot_layout(guide = "collect", heights = c(1, 1, margin=0.5))
  #                 (( plots[[1]] + theme(legend.position = "right")                               | plots[[2]] + theme(legend.position = "none") | plots[[3]] + theme(legend.position = "none")) / (plots[[4]]  + theme(legend.position = "none") | plots[[5]]  + theme(legend.position = "none") | plots[[6]] + theme(legend.position = "none"))               ) + plot_layout(guide = "collect")
# if running interactively, look at it:
print(combined_plots)

# Save file
ggsave(combined_plots, file = "v1_3scn.pdf", device = "pdf", dpi = 100, width = 40, height = 30, units = "cm")

# MODEL V2

options(scipen=10000)
plots <- plot_scenarios_simple_zoo(window(result2$result, start=as.Date('2020/01/20'), end=as.Date('2020/12/01'),format="%Y/%m/%d"),
#plots <- plot_scenarios_simple_zoo(result2$result,
                                   scn.labels, data='srag',
                                   cols=c('I', 'E', 'R', 'C', 'H', 'U', 'M', 'MD','Q'),#c('I', 'r', 'C', 'H', 'U', 'MD',"Reff_cori","Reff_como"),
                                   ylabels=ylabels, lang=lang, end.data=as.Date('2020/08/31'),
                                   hlines=c(r=1-1/2.5,
                                            U=unname(parameters['icu_beds_available']),
                                            H=unname(parameters['beds_available'])),
                                   label_months=T, line.size=2,color.palette="Dark2")
        #scale_x_date(#labels = date_format("%Y/%m/%d"), 
                          #breaks = date_breaks("months"), 
                          #limits = c(as.Date('2020/01/20'), as.Date('2020/12/01')))
# or plot only a certain window
# plots <- plot_scenarios_simple_zoo(window(result, start=as.Date('2020-3-1'), end=as.Date('2020-8-1')), scn.labels)
#result2$results
# put the plots in a grid
#six_plots <- plots[c(1,4:7,9)] #select which 6 variables to plot
#six_plots <- plots[c(1,4,5,6,7,9)] #select which 6 variables to plot
#combined_plots <- combine.plots(six_plots, legend.pos='bottom')#, margin=0.2+0.1*length(scenarios))
four_plots <- plots[c(5,4,7,8)] #select which 6 variables to plot
combined_plots <- (( plots[[4]] + theme(legend.position = "bottom", legend.direction="vertical") | plots[[7]] + theme(legend.position = "none")                                               ) / ( plots[[5]] + theme(legend.position = "none") | plots[[8]]  + theme(legend.position = "none")                                               ) / guide_area()) + plot_layout(guide = "collect", heights = c(1, 1, margin=0.5))
#                 (( plots[[1]] + theme(legend.position = "right")                               | plots[[2]] + theme(legend.position = "none") | plots[[3]] + theme(legend.position = "none")) / (plots[[4]]  + theme(legend.position = "none") | plots[[5]]  + theme(legend.position = "none") | plots[[6]] + theme(legend.position = "none"))               ) + plot_layout(guide = "collect")
# if running interactively, look at it:
print(combined_plots)

# Save file
ggsave(combined_plots, file = "v2_3scn.pdf", device = "pdf", dpi = 100, width = 40, height = 30, units = "cm")

































#############################################
#result_analysis(results,parameters)
combine.plots.R_eff(plots,legend.pos='bottom')

# plots com legenda dentro do quadro:

# Reff
plots[[1]] + #xlim(as.Date("2020-05-01"),as.Date("2020-07-31")) +
  theme(legend.key.size = unit(2.3, "cm"),
        legend.justification=c(0,1),
        legend.position=c(0.5, 0.96),
        legend.text=element_text(size=16))
#legend.background = element_blank(),
#legend.key = element_blank   ())
# ggsave(file="casos_1jun_20jul_analise.png", width=12, height=9)

# Plot Reff (Cori_new) for different scenarios
plot_scenarios_reff(result, scenarios, scn.labels = c(3,7,8))

### Plotar interven??es ao longo do tempo


base_parameters = c()
parameters <- update_parameters(base_parameters)
interventions <- update_scenario_parameters(scenario=49, parameters)
cidade = "Goiania"
# Plotar interven??es ao longo do tempo
pdf("Intervention_scenario_49.pdf", height = 7, width = 12)
plot.all.interventions(interventions, parameters, "Goiânia, GO", "Scenario 49",
                       line_sce_name = -3.2,cidade = cidade,
                       x_lim = c("2020-02-20", "2021-02-01"))
dev.off()
## COMMENTING PREVIOUSLY USED PLOTS:
#
# #plot mortes diarias com limiar em 400 mortes/dia
# plots[[6]] + xlim(as.Date("2020-05-01"),as.Date("2020-07-31")) +
#   theme(#legend.key.size = unit(1.5, "cm"),
#     legend.justification=c(0,1), 
#     legend.position=c(0.01, 0.99)) +
#   geom_hline(yintercept=800, linetype="dashed", color="black") +
#   geom_hline(yintercept=400, linetype="dashed", color="blue") +
#   annotate(geom="text", x=as.Date("2020-05-05"), y=400-75, label="400",
#            color="blue")+
#   annotate(geom="text", x=as.Date("2020-05-05"), y=800-75, label="800",
#            color="black")
# #ggsave(file="obitos_lock_relax_cut.png", width=12, height=9)
# 
# 
# #plot UTIs com limiar em 400 mortes/dia
# plots[[5]] + xlim(as.Date("2020-05-01"),as.Date("2020-07-31")) +
#   theme(#legend.key.size = unit(1.5, "cm"),
#     legend.justification=c(0,1), 
#     legend.position=c(0.01, 0.99)) +
#   geom_hline(yintercept=1335, linetype="dashed", color="gray")+
#   annotate(geom="text", x=as.Date("2020-05-05"), y= parameters["icu_beds_available"]-450, label=parameters["icu_beds_available"],
#            color="black") +
#   annotate(geom="text", x=as.Date("2020-05-05"), y=1335-450, label="1335",
#            color="gray")
# #ggsave(file="UTIs_lock_relax_cut.png", width=12, height=9)
