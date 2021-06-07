library(plyr)
library(dplyr)
library(mvtnorm)
library(doParallel)
library(ggplot2)
library(googlesheets4)
library(ggpubr)
#######
# setwd("codes/")
# source files only if needed
if(!exists("run_covid19", mode="function")) source("Brazilian_Model_v2b.R")
if(!exists("result_ts", mode="function")) source('parameter_calibration_functions.R')
if(!exists("test_modifying_parameters", mode="function")) source("test_modifying_parameters.R")
source("dados_municipios.R")
city.zoo.list <- carregar_dados_cidades(estado = "SP", cidade = "Sao_Paulo")
for(i in 1:length(city.zoo.list)) assign(names(city.zoo.list)[i], city.zoo.list[[i]])
#####read data


generate_samples <- function(fit,params = c("startdate","p","perc_threshold","home_steep"),n = 100){
  all_samples <- c()
  
  if("startdate" %in% params){####faz sampling da data
    cdf_date <- cumsum(fit$prob)
    date_samples <- as.data.frame(plyr::count(sample(fit$startdate,size = n,replace=T,prob=fit$prob)))
    colnames(date_samples) <- c("startdate","freq")
    for(date in date_samples$startdate){
      freq <- date_samples[date_samples$startdate == date,"freq"]
      vals <- fit[fit$startdate == date,]
      if(length(params)-1 > 1){
        cov_name <-paste0("cov_matrix",rep(1:(length(params)-1)**2))
        cov_mat <- dplyr::select(vals,cov_name)
        cov_mat <- as.numeric(cov_mat)
        cov_mat <- matrix(cov_mat,nrow = length(params)-1,ncol=length(params)-1)
        mean <- as.numeric(dplyr::select(vals,params[params!="startdate"]))
        samples <- rmvnorm(freq,mean = mean,sigma = cov_mat, method="svd")
        samples <- as.data.frame(samples)
      }
      else{
        mean <- as.numeric(dplyr::select(vals,params[params!="startdate"]))
        cov_mat <- as.numeric(dplyr::select(vals,"cov_matrix"))
        samples <- rnorm(freq,mean = mean,sd = cov_mat)
        samples <- as.data.frame(samples)
      }
      samples$startdate <- date
      
      colnames(samples)<-c(params[params!="startdate"],"startdate")
      samples <- samples[c(length(params),1:length(params)-1)]
      all_samples <- rbind(all_samples,samples)
      }
  }
  else{
    cov_name <-paste0("cov_matrix",rep(1:(length(params))**2))
    
  }
  all_samples<-as.data.frame(all_samples)
  c(all_samples)
}


generate_solutions <- function(samples,base_parameters = c(),scenario = 53,enddate = as.Date("2020-08-31"),n.cores = 1,diff = T){
  base_parameters <- update_parameters(base_parameters)
  base_parameters["stopdate"] <- as.numeric(enddate)
  base_interventions <- update_scenario_parameters(base_parameters, estado = "SP", cidade = "Sao_Paulo")
  interventions <- reload_scenario_parameters(base_parameters, scenario = scenario, base_interventions)
  registerDoParallel(cores = n.cores)
  results <-foreach (i = 1:nrow(samples)) %dopar% {
    sample <- samples[i,]
    parameters2 <- base_parameters
    interventions2 <- interventions
    for (j in names(sample)){parameters2[j] <- sample[j]} ###assume startdate numerica
    parameters2<- unlist(parameters2)
    
    # Update startdate e enddate time difference
    dif = base_parameters["startdate"] - parameters2["startdate"]
    
    for(k in 1:length(interventions2)){
    not_zero <- interventions[[k]]["startdate",] != 0
    interventions2[[k]]["startdate",not_zero] <- interventions[[k]]["startdate",not_zero] + dif
    interventions2[[k]]["enddate",] <- interventions[[k]]["enddate",] + dif
    }
    
    out <- run_covid19(parameters2,interventions2,Y)
    list_cols <- list(C=Cindex+1,M=CMindex+1)
    result <- result_ts(out[[2]], parameters2,parameters2["startdate"], diff=F,
                        cols=list_cols)
    result <- as.list(result)
    c(result)
  }
  Cases <- zoo(results[[1]]$C)
  Deaths <- zoo(results[[1]]$M)
  for(i  in 2:nrow(samples)){
    Cases <- merge(Cases,zoo(results[[i]]$C))
    Deaths <- merge(Deaths,zoo(results[[i]]$M))
  }
  Cases <- na.fill(Cases,fill = 0)
  Deaths <- na.fill(Deaths,fill = 0)
  if(diff){
    Cases <- apply(Cases,2,diff)
    Deaths <- apply(Deaths, 2, diff)
  }
  Cases <- apply.weekly(Cases, function(x) apply(x,2,sum) )
  Deaths <- apply.weekly(Deaths,function(x) apply(x,2,sum) )
  C <- data.frame(mean = apply(Cases, 1, mean),
                  sd = apply(Cases, 1, sd),
                  q2.5 = apply(Cases, 1, quantile, 0.025),
                  q50 = apply(Cases, 1, quantile, 0.5),
                  q97.5 = apply(Cases, 1, quantile, 0.975))
  M <- data.frame(mean = apply(Deaths, 1, mean),
                  sd = apply(Deaths, 1, sd),
                  q2.5 = apply(Deaths, 1, quantile, 0.025),
                  q50 = apply(Deaths, 1, quantile, 0.5),
                  q97.5 = apply(Deaths, 1, quantile, 0.975))
  result <- merge(as.xts(C),as.xts(M))
  colnames(result) <- c("Cases.Mean","Cases.sd","Cases.q2.5",
                             "Cases.q50","Cases.q97.5","Deaths.Mean",
                             "Deaths.sd","Deaths.q2.5","Deaths.q50","Deaths.q97.5")
  c(result)
}


plot.CI <- function(sol, diff = FALSE){
  data.cases <- window(xts(now.srag.zoo),end=max(index(sol)))
  data.deaths <- window(xts(now.obito.srag.zoo),end=max(index(sol)))
  if(diff){
    data.cases2 <- diff(data.cases)
    data.deaths2 <- diff(data.deaths)
    data.cases2[1,] <- data.cases[1,]
    data.deaths2[1,] <- data.deaths[1,]
    sol.diff <- apply(sol,2,diff)
    sol.diff <- apply.weekly(sol.diff,function(x){c(sum(x$Cases.Mean),sum(x$Cases.sd),sum(x$Cases.q2.5),
                                       sum(x$Cases.q50),sum(x$Cases.q97.5),sum(x$Deaths.Mean),
                                       sum(x$Deaths.sd),sum(x$Deaths.q2.5),sum(x$Deaths.q50),sum(x$Deaths.q97.5))})
    df <- data.frame(coredata(sol.diff))
    colnames(df) <- c("Cases.Mean","Cases.sd","Cases.q2.5",
                           "Cases.q50","Cases.q97.5","Deaths.Mean",
                           "Deaths.sd","Deaths.q2.5","Deaths.q50","Deaths.q97.5")
    rownames(df)<- as.Date(rownames(df))
    df$Date <- as.Date(rownames(df))
    df <- df[1:dim(df)[1]-1,]
    # sol.diff <- fortify(sol.diff)
    p1 <- ggplot(df, aes(x = Date, y = Cases.Mean)) + geom_line() +
      geom_ribbon(aes(x = Date, ymin = Cases.q2.5,ymax = Cases.q97.5),alpha = 0.1)+
      labs(x = "Date", y = "New Cases") + ylim(0,NA)+geom_point(data=data.cases2, aes(x = Index, y = x))
    p2 <- ggplot(df, aes(x = Date, y = Deaths.Mean)) + geom_line() +
      geom_ribbon(aes(x = Date, ymin = Deaths.q2.5,ymax = Deaths.q97.5),alpha = 0.1)+
      labs(x = "Date", y = " New Deaths")+ylim(0,NA)+geom_point(data=data.deaths2, aes(x = Index, y = x))
  }
  else{
    p1 <- ggplot(sol, aes(x = Index, y = Cases.Mean)) + geom_line() + 
      geom_ribbon(aes(x = Index, ymin = Cases.q2.5,ymax = Cases.q97.5),alpha = 0.1)+
      labs(x = "Date", y = "Cases")+geom_point(data=data.cases,aes(x = Index, y = x))
    p2 <- ggplot(sol, aes(x = Index, y = Deaths.Mean)) + geom_line() + 
      geom_ribbon(aes(x = Index, ymin = Deaths.q2.5,ymax = Deaths.q97.5),alpha = 0.1)+
      labs(x = "Date", y = "Deaths")+geom_point(data=data.deaths,aes(x = Index, y = x))
    
  }
  grid.arrange(p1,p2,ncol=2,nrow=1)
}

plot.CI.all <- function(sol, diff = FALSE){
  data.cases <- window(xts(now.srag.zoo),end=max(sol$date))
  data.deaths <- window(xts(now.obito.srag.zoo),end=max(sol$date))
  if(diff){
    data.cases2 <- diff(data.cases)
    data.deaths2 <- diff(data.deaths)
    data.cases2[1,] <- data.cases[1,]
    data.deaths2[1,] <- data.deaths[1,]
    # sol.filter <- split(sol,sol$model)
    # df <- c()
    # for(i in 1:length(sol.filter)){
    #   values <- sol.filter[[i]]
    #   name_model <- names(sol.filter)[i]
    #   dates <- dplyr::select(values,"date")
    #   values <- dplyr::select(values,c("Cases.Mean","Cases.sd","Cases.q2.5",
    #                             "Cases.q50","Cases.q97.5","Deaths.Mean",
    #                             "Deaths.sd","Deaths.q2.5","Deaths.q50","Deaths.q97.5"))
    #   values <- xts(values, order.by = as.Date(dates$date))
    #         sol.diff <- apply(values,2,diff)
    #   sol.diff <- apply.weekly(sol.diff,function(x){c(sum(x$Cases.Mean),sum(x$Cases.sd),sum(x$Cases.q2.5),
    #                                                   sum(x$Cases.q50),sum(x$Cases.q97.5),sum(x$Deaths.Mean),
    #                                                   sum(x$Deaths.sd),sum(x$Deaths.q2.5),sum(x$Deaths.q50),sum(x$Deaths.q97.5))})
    #   dates <- rownames(sol.diff)
    #   sol.diff <- data.frame(coredata(sol.diff))
    #   colnames(sol.diff) <- c("Cases.Mean","Cases.sd","Cases.q2.5",
    #                     "Cases.q50","Cases.q97.5","Deaths.Mean",
    #                     "Deaths.sd","Deaths.q2.5","Deaths.q50","Deaths.q97.5")
    #   sol.diff$Date <- as.Date(dates)
    #   sol.diff$model <- name_model
    #   sol.diff <- sol.diff[1:dim(sol.diff)[1]-1,]
    #   rownames(sol.diff) <- NULL
    #   df <- rbind(df,sol.diff)
    # }
    
    print(sol)
    p1 <- ggplot(sol, aes(x = date, y = Cases.Mean, group = model, color = model, fill = model)) + geom_line() +
      geom_ribbon(aes(x = date, ymin = Cases.q2.5,ymax = Cases.q97.5),alpha = 0.3,color=NA)+
      labs(x = "Date", y = "New Cases") + ylim(0,NA)+
      geom_point(data=data.cases2,aes(x = Index, y = x, group = NULL,color = NULL, fill=NULL))
    p2 <- ggplot(sol, aes(x = date, y = Deaths.Mean, group = model,color = model, fill = model)) + geom_line() +
      geom_ribbon(aes(x = date, ymin = Deaths.q2.5,ymax = Deaths.q97.5),alpha = 0.3,color=NA)+
      labs(x = "Date", y = " New Deaths")+ylim(0,NA)+
      geom_point(data=data.deaths2,aes(x = Index, y = x,group = NULL,color=NULL,fill = NULL))
  }
  else{
    p1 <- ggplot(sol, aes(x = date, y = Cases.Mean, group = model, color = model)) + geom_line() +
      geom_ribbon(aes(x = date, ymin = Cases.q2.5,ymax = Cases.q97.5),alpha = 0.1)+
      labs(x = "Date", y = "Cases")+geom_point(data=data.cases,aes(x = Index, y = x))
    p2 <- ggplot(sol, aes(x = date, y = Deaths.Mean, group = model, color = model)) + geom_line() +
      geom_ribbon(aes(x = date, ymin = Deaths.q2.5,ymax = Deaths.q97.5),alpha = 0.1)+
      labs(x = "Date", y = "Deaths")+geom_point(data=data.deaths,aes(x = Index, y = x))

  }
  ggarrange(p1,p2,ncol=2,nrow=1,common.legend = TRUE)
}

#######TO RUN #########
#data<- read.csv("fits/new_model1_srag.1.CM.log0-2021_05_17_SP.csv") ###resultado percolacao
# data2<- read.csv("fits/new_model2_srag.1.CM.log0-2021_05_17_SP.csv") ###resultado sem perc (pequenos valores rodar novo)
#data2 <- read.csv("fits/model2_scn50_srag_2021_05_24_SP.csv")
#data3 <- read.csv("fits/model3_scn54_30npi_srag_2021_05_24_SP.csv")
# data2 <- read.csv("fits/SP/new3/new_model2_srag.1.CM.log0-2021_05_17_SP.csv")
#result1 <- generate_solutions(as.data.frame(generate_samples(data,n=100)),n.cores = 4)
#samples <- as.data.frame(generate_samples(data2,params=c("startdate","p"),n=100))
#samples$perc_threshold <- 0
#samples$home_steep <- 0
#samples$home_eff <- 0
#result2 <- generate_solutions(samples,n.cores = 4)
#samples <- as.data.frame(generate_samples(data3,params=c("startdate","p"),n=100))
#samples$perc_threshold <- 0
#samples$home_steep <- 0
#samples$home_eff <- 0
#result3 <- generate_solutions(samples,scenario = 58,n.cores = 4)
# result1 <- read.csv("fits/SP/resultados_CI/resultado_CI_model_percolation.csv")
# result2 <- read.csv("fits/SP/resultados_CI/resultado_CI_model_standard.csv")
# result3 <- read.csv("fits/SP/resultados_CI/resultado_CI_model_standard_2.csv")
# plot.CI(result1,diff = T) ## plota apenas um resultado
#result1.2 <- data.frame(date = as.Date(index(result1)),coredata(result1)) ###shenanigans pra transformar em dataframe
#result1.2$model <- "Percolation"
#result2.2 <- data.frame(date = as.Date(index(result2)),coredata(result2))
#result2.2$model <- "Standard"
#result3.2 <- data.frame(date = as.Date(index(result3)),coredata(result3))
#result3.2$model <- "Standard+30"
# result1$date <- as.Date(result1$date)
# result2$date <- as.Date(result2$date)
# result3$date <- as.Date(result3$date)
#result.both <- rbind(result1.1,result1.2)
#result.both <- result.both[1:(dim(result.both)[1]-1),]
# result.both.2 <- rbind(result1,result2,result3)
#plot.CI.all(result.both,diff=T) ### diff falso não verificado





