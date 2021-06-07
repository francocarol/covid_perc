################################################
################################################
################################################

# Code to run model fittings for models 1, 2 and 3

################################################
################################################
################################################

list_packages <- c("colorout")
# download an install packages that are not available locally
(new_packages <- list_packages[!(list_packages %in% installed.packages()[,'Package'])])  
if(length(new_packages)!=0) install.packages(new_packages, dependencies = TRUE, type = "binary")

# All necessary source codes are in this folder:
# REMEMBER TO ERASE IT WHEN ZIPPING IT AND SHARING IT 
#setwd("/home/caroline/Documents/paper_covid/codes/")
#IN CLUSTER
setwd("~/extraspace/Carol/paper_covid/codes")

################################################
# INPUTS:
################################################

### City and state:
estado = "SP"
cidade = "Sao_Paulo"
source("least_squares_fitting_homesteep_v2.R")
### Language and locale
Sys.setlocale("LC_TIME", "en_GB.UTF-8")
# name of the files, where we will save the fittign outputs
filename1 = paste0('fits/model1_scn50_srag_', data.base, "_",estado, ".csv")
filename2 = paste0('fits/model2_scn60_srag_', data.base, "_",estado, ".csv")
filename3 = paste0('fits/model3_scn59_30npi_srag_', data.base, "_",estado, ".csv")
# intervention scenario used for each model version
# 50 business as usual in SP
# 54 +30% NPIs
# 59 tentativa: cenario sem casos importados e travel ban + 30% NPIs
# 60 tentativa: cenario sem casos importados e travel ban
# 61 tentativa: cenario com 50% mais NPIs
Scn1 <- 50
Scn2 <- 60
Scn3 <- 59
# startdate at the begining of fitting range (end = startdate + 45 days)
startdate_0 <- as.Date("2020-01-10") # CF mudei tudo para 01 janeiro e depois para 10 janeiro # JANELA MAIOR #("2020-01-01") 

################################################
################################################
################################################

######
#Model version 1 (full model: p, T_perc, home_steep, startdate != 0)
######

# fitting algorithm for model 1 (with percolation)
source("least_squares_fitting_homesteep_v2.R")

parameters1 <- update_parameters(c(), estado = estado, cidade = cidade) # Aqui especifico cidade e estado
parameters1["startdate"] <- as.numeric(startdate_0) # CF mudei tudo para 01 janeiro # A data está iniciando em 20/jan agora, e esse valor precisa bater com o data_start no script de least_squares_fitting
parameters <- parameters1

interventions <- update_scenario_parameters(scenario=Scn1, parameters) # O último cenário para SP é o 50 (Marcelo)
# CF: aqui fitamos até 30 de agosto de 2020 (em 1 Set as restricoes comecam a ser levantadas desorganizadamente) max.date = 18505
run_the_fitting(filename = filename1, max.date = 18505)

######
#Model version 2 (nested model: T_perc = home_steep = 0; p, startdate != 0)
######

# fitting algortithm for models 2 and 3 (perc_threshold=home_steep=home_eff=0)
source("model2_least_squares_fitting_homesteep_v2.R")

parameters2 <- update_parameters(c(), estado = estado, cidade = cidade)
parameters2 <- replace(parameters2, c("perc_threshold","home_steep", "home_eff"), c(0,0,0))# all perc_threshold=home_steep=0
parameters2["startdate"] <- as.numeric(startdate_0)
parameters <- parameters2

#interventions <- update_scenario_parameters(scenario=Scn2, parameters)
#run_the_fitting(filename = filename2, max.date = 18505, data.start = startdate_0)

######
#Model version 3 (T_perc = home_steep = 0; p, startdate != 0 and 30% MORE NPI coverage (different intervention scenario))
######

# fitting algortithm for models 2 and 3 (perc_threshold=home_steep=home_eff=0)
source("model2_least_squares_fitting_homesteep_v2.R")

parameters3 <- update_parameters(c(), estado = estado, cidade = cidade)
parameters3 <- replace(parameters3, c("perc_threshold","home_steep", "home_eff"), c(0,0,0))# all perc_threshold=home_steep=0
parameters3["startdate"] <- as.numeric(startdate_0)
parameters <- parameters3

#interventions <- update_scenario_parameters(scenario=Scn3, parameters)
#run_the_fitting(filename = filename3, max.date = 18505, data.start = startdate_0)