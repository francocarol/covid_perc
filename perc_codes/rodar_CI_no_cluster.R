source("CI_estimator.R")

##############CHECAR SE ESTÁ NO FOLDER CODES##################

data<- read.csv("fits/new_model1_srag.1.CM.log0-2021_05_17_SP.csv") ###resultado percolacao
data2 <- read.csv("fits/model2_scn50_srag_2021_05_24_SP.csv")
data3 <- read.csv("fits/model3_scn54_30npi_srag_2021_05_24_SP.csv")
n.runs <- 1000 ####10 mil dá bom
ncores <- detectCores()-1
##########################MODEL 1####################################
# samples <- as.data.frame(generate_samples(data,n = n.runs))
# result <- generate_solutions(samples,n.cores = ncores)
# result <- data.frame(date = index(result),coredata(result)) ###shenanigans pra transformar em dataframe
# result$model <- "Percolation"
# write.csv(result, "fits/SP/resultados_CI/resultado_CI_model_percolation.csv",  row.names = FALSE)
# ##########################MODEL 2####################################
# samples <- as.data.frame(generate_samples(data2,params = c("startdate","p"),n = n.runs))
# samples$perc_threshold <- 0
# samples$home_steep <- 0
# samples$home_eff <- 0
# result <- generate_solutions(samples,n.cores = ncores)
# result <- data.frame(date = index(result),coredata(result)) ###shenanigans pra transformar em dataframe
# result$model <- "Standard"
# write.csv(result, "fits/SP/resultados_CI/resultado_CI_model_standard.csv",  row.names = FALSE)
###########################MODEL 3###################################
samples <- as.data.frame(generate_samples(data3,params = c("startdate","p"),n = n.runs))
samples$perc_threshold <- 0
samples$home_steep <- 0
samples$home_eff <- 0

result <- generate_solutions(samples,n.cores = ncores,scenario = 58)
result <- data.frame(date = index(result),coredata(result)) ###shenanigans pra transformar em dataframe
result$model <- "Standard_2"
write.csv(result, "fits/SP/resultados_CI/resultado_CI_model_standard_2.csv",  row.names = FALSE)
