if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(zeallot)){install.packages("zeallot"); library(zeallot)} # provides %<-%
if(!require(FME)){install.packages("FME"); library(FME)}
if(!require(zeallot)){install.packages("zeallot"); library(zeallot)} # provides %<-%/
if(!require(minpack.lm)){install.packages("minpack.lm"); library(minpack.lm)}
if(!require(doParallel)){install.packages("doParallel"); library(doParallel)}
if(!require(MASS)){install.packages("MASS"); library(MASS)}
if(!require(foreach)){install.packages("foreach"); library(foreach)}
if(!require(pse)){install.packages("pse"); library(pse)}

if(!exists("sheets_read",mode="function")) library(googlesheets4)
if(!exists("test_modifying_parameters", mode="function")) source("test_modifying_parameters.R")
if(!exists("result_ts", mode="function")) source('parameter_calibration_functions.R')

### Definir cidade:
#estado = "GO"; cidade = "Goiania"
estado =  "SP"; cidade = "Sao_Paulo"
#estado =  "RS"; cidade = "Porto_Alegre"

source("Brazilian_Model_v2b.R")
source("import_functions_v2.R")
source("dados_municipios.R")

# Load srag and covid data from specified city
city.zoo.list <- carregar_dados_cidades(estado = estado, cidade = cidade)
for(i in 1:length(city.zoo.list)) assign(names(city.zoo.list)[i], city.zoo.list[[i]])

################################################
##### Home_steep como parametro de fitting
################################################

fit.results <- read.csv(paste0('fits/', get_latest_fit(estado = estado))[i]); print(paste0('fits/', get_latest_fit(estado = estado))[1])

# Load parameters
base_parameters = c()
parameters <- update_parameters(base_parameters, estado = estado, cidade = cidade)

# Best fitting result
dopt <- fit.results[fit.results$residuo == min(fit.results$residuo),]

# Update parameters
pars_to_update <- names(parameters)[names(parameters) %in% names(dopt)]
parameters[pars_to_update] <- unlist(dopt[pars_to_update])

# Load interventions
interventions <- update_scenario_parameters(scenario=50, parameters)

# Write image output
#filename = paste0("fitting_SP.png")
#png(filename, height = 480*3, width = 480*5, res = 200)
test_modifying_parameters(parameters, interventions, Y, unlist(dopt[1:4]), flag='srag',
                          rollwindow = 1, weekly = TRUE)
#dev.off()
