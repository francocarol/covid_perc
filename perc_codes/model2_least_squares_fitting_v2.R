if(!require(FME)){install.packages("FME"); library(FME)}
if(!require(zeallot)){install.packages("zeallot"); library(zeallot)} # provides %<-%/
if(!require(minpack.lm)){install.packages("minpack.lm"); library(minpack.lm)}
if(!require(doParallel)){install.packages("doParallel"); library(doParallel)}
if(!require(MASS)){install.packages("MASS"); library(MASS)}
if(!require(foreach)){install.packages("foreach"); library(foreach)}
if(!require(pse)){install.packages("pse"); library(pse)}
if(!require(xts)){install.packages("xts"); library(xts)}

# source files only if needed
if(!exists("run_covid19", mode="function")) source("Brazilian_Model_v2b.R")
if(!exists("update_parameters", mode="function")) source("import_functions_v2.R")
if(!exists("result_ts", mode="function")) source('parameter_calibration_functions.R')
if(!exists("test_modifying_parameters", mode="function")) source("test_modifying_parameters.R")

if(!exists("carregar_dados_cidades",mode="function")) source("dados_municipios.R")
# Load srag and covid data from specified city
city.zoo.list <- carregar_dados_cidades(estado = "SP", cidade = "Sao_Paulo")
for(i in 1:length(city.zoo.list)) assign(names(city.zoo.list)[i], city.zoo.list[[i]])

###
trim_intervention <- function(interventions, parameters){
  ceiling <- parameters["stopdate"] - parameters["startdate"] + 1
  
  for(i in 1:length(interventions)) {
    f <- interventions[[i]]
    if(!is.na(f["startdate",1])){
      if(f["startdate",1] <= ceiling) {
        val <- f[1:nrow(f), f["startdate",] <= ceiling]
        m <- matrix(val, dim(f))
        rownames(m) <- rownames(f)
        m["enddate",ncol(m)] <- ceiling
        interventions[[i]] <- m 
      } else {
        f[,] <- NA
        interventions[[i]] <- f
      }
    }
  }
  return(interventions)
}

###
# use this to avoid reloading parameters (in case of connection limit)
# REMOVE the file to force this to refresh
#if(! file.exists('parameters.R')){
#  parameters<-c()
#  parameters <- update_parameters(parameters, 9)
#  interventions <- update_scenario_parameters(scenario=18, parameters)
#  # save to local file
#  dump(c('parameters', 'date_parameters'), 'parameters.R')
#} else {
#  source('parameters.R')
#}

###############ELLIPTICAL CONFIDENCE INTERVAL ESTIMATION###############
confidence_interval_estimate <- function(v_matrix, par, dof, CI=0.95){
  # dof: degrees of freedom = # data points - # params
  S = matrix(v_matrix, nrow=length(par))
  Sinv = ginv(S)
  r=sqrt(qf(CI, 2, dof)*2)
  theta=seq(0,2*pi,length.out=100)
  z=cbind(r*cos(theta),r*sin(theta))
  # transform points of circle into points of ellipse using
  # svd of inverse covariance matrix
  Sinv_svd=svd(Sinv) # inverse of covariance matrix
  # transform from circle to ellispse
  xt=t(Sinv_svd$v)%*%diag(1/sqrt(Sinv_svd$d))%*%t(z)
  x=t(xt)
  # translate the ellipse so that center is the estimated parameter value
  x=x+matrix(rep(as.numeric(par),100),nrow=100,byrow=T)
  
  plot(x[,1],x[,2],type="l",xlab="p",ylab="perc threshold",lwd=2)
  points(par[1],par[2],pch=20,col="blue",cex=2)
}

########################RESIDUALS COMPUTATION FUNCTION######################
least_squares_solve <- function(params,date_params, flag = "covid", diff=0,
                                use=c('C','M'), logscale=TRUE, rollwindow=0,
                                interventions = interventions, weekly =TRUE){
  #############passa o parametro p para o conjunto de parametros############
  parameters2 <-parameters
  
  for (i in names(params)){parameters2[i] <- params[i]}
  for (i in names(date_params)){parameters2[i] <- as.numeric(date_params[i])}
  
  #####################DADOS DE MORTE DE SP#################################
  if(flag == "covid"){
    data_c <- now.covid.zoo
    data_cm <- now.obito.covid.zoo
  }
  if(flag == "srag"){
    data_c <- now.srag.zoo
    data_cm <- now.obito.srag.zoo
  }
  if (diff){
    data_c <- diff(data_c)
    data_cm <- diff(data_cm)
  }
  if (rollwindow > 1){
    data_c <- zoo::rollmean(data_c, rollwindow)
    data_cm <- zoo::rollmean(data_cm, rollwindow)
  }
  
  parameters2["stopdate"] <- max(time(data_c), time(data_cm))
  parameters2["reportc"] <- 0
  interventions <- trim_intervention(interventions, parameters2)
  ##########################################################################
  out <- run_covid19(parameters2,interventions,Y)
  res <- result_ts(out[[2]], parameters2, parameters2["startdate"], diff=F)
  
  if (diff) res <- diff(res)
  
  if (weekly){
    model_c_xts <- xts::apply.weekly(res$C, sum)
    model_cm_xts <- xts::apply.weekly(res$CM, sum)
    
    model_c <- zoo(model_c_xts, as.Date(index(model_c_xts)))
    model_cm <- zoo(model_cm_xts, as.Date(index(model_cm_xts)))
    
  } else {
    model_c <- res$C
    model_cm <- res$CM
  }
  
  c.res <- merge.zoo(model=model_c, obs=data_c, all=F)
  m.res <- merge.zoo(model=model_cm, obs=data_cm, all=F)
  # Manual correction of sub-notification of deaths
  # if(flag == "covid")
  #   m.res <- 0.4*m.res
  
  residuo <- c()
  if ('C' %in% use){
    if (logscale){
      residuo <- c(residuo, log(as.vector(c.res$model)+0.1)-log(as.vector(c.res$obs)+0.1))
    } else {
      residuo <- c(residuo, (as.vector(c.res$model)-as.vector(c.res$obs))/
                     mean(c.res$obs) / length(c.res$obs))
    }
  }
  if ('M' %in% use){
    if (logscale){
      residuo <- c(residuo, log(as.vector(m.res$model)+0.1)-log(as.vector(m.res$obs)+0.1))
    } else {
      residuo <- c(residuo, (as.vector(m.res$model)-as.vector(m.res$obs))/
                     mean(m.res$obs) / length(m.res$obs))
    }
  }
  return(residuo)
}
######################MAIN SOLVER FUNCTION#################################
solver_wrapper_fitval <- function(flag_ = "covid", diff=0, use=c('C', 'M'),
                                  logscale=TRUE, rollwindow=0, weekly = TRUE){
  params <- c(p = 0.015,perc_threshold = 0)
  lower <- c(p = 0, perc_threshold = 0)
  upper <- c(p = 0.2, perc_threshold = 0)
  data_start <- as.Date("2020-02-01")
  guesses <- list(c(p=0.025, perc_threshold=0),
                  c(p=0.02, perc_threshold=0),
                  c(p=0.025, perc_threshold=0),
                  c(p=0.02, perc_threshold=0),
                  c(p=0.015, perc_threshold=0),
                  c(p=0.03, perc_threshold=0),
                  c(p=0.02, perc_threshold=0),
                  c(p=0.025, perc_threshold=0),
                  c(p=0.02, perc_threshold=0),
                  c(p=0.015, perc_threshold=0)
  )
  # run in parallel in linux
  if(Sys.info()[["sysname"]] == "Linux"){
    registerDoParallel(cores=detectCores()-1)
  }
  # try *harder*
  #control <- nls.lm.control(ptol=1.490116e-12, ftol=1.490116e-12,
  #                          maxiter=200, maxfev=4*100*(length(par) + 1))
  control <- nls.lm.control() # default values
  
  results <-  foreach(j = 0:20,.combine = rbind) %:% foreach(i = 1:length(guesses)) %dopar% {
    print(c(i,j))
    params <- guesses[[i]]
    date_params <- c(startdate = as.Date(data_start + j))
    ###
    interventions_paral <- interventions
    
    ## To do: o ajuste da data de intervenção (lin 156) está em relação à data de 01-fev-2020. Por sua vez,
    ## este valor deve ser passado para o parâmetro startdate antes da tabela de intervenção (interventions) ser carregada no script.
    ## Por exemplo:
    ## parameters["startdate"] <- as.numeric(as.Date("2020-02-01"))
    ## interventions <- update_scenario_parameters(scenario=20, parameters)
    ## Uma outra possibilidade é implementar o cálculo da diferença em relação a qualquer data que tenha sido usada
    ## ao se carregar 'interventions':
    ## startdate_dif = parameters["startdate"] - as.numeric(as.Date(data_start + j))
    ## O valor de startdate_dif deve ser substituído por '- j' da linha 156, prestando atenção ao sinal (se subtração ou adição)
    ## dependendo se este valor será maior ou menor do que o primeiro valor de startdate a ser avaliado pelo fitting.
    
    for(int in 1:length(interventions_paral)) {
      if(!any(is.na(interventions_paral[[1]]))){
        interventions_paral[[int]][c("startdate", "enddate"),] = interventions_paral[[int]][c("startdate", "enddate"),] - j
        interventions_paral[[int]]["startdate", which(interventions_paral[[int]]["startdate",] < 0)] <- 0
      }
    }
    ###
    fitval <- nls.lm(par = params, fn = least_squares_solve, upper = upper,
                     lower=lower, control=control, date_params=date_params,
                     flag=flag_, diff=diff, use=use, logscale=logscale,
                     rollwindow=rollwindow, interventions = interventions_paral,
                     weekly = weekly)
    #print(fitval)
    
    # residuo = sum of squared residuals
    residuo<- sum(residuals(fitval)*residuals(fitval))
    # residuo_mean = (1/N) * sum of sqrd residuals
    residuo_mean <- residuo/length(residuals(fitval))
    # NEG LOG likelyhood = N * ln ( residuo_mean ) # ACTUALLY NEGATIVE LOG LIKELYHOOD !!!
    # it might be <0 because we are applying log to residuals which are <1
    # smaller residuals --> higher log likelyhood --> smaller negative log likelyhood
    likelihood <- length(residuals(fitval))*log(residuo_mean)
    # 
    cov_matrix <- tryCatch({vcov(fitval)},error = function(cond){return(as.numeric(rep(NA, length(params)^2)))})
    print(c(i,j,"ended"))
    c(startdate = as.numeric(date_params["startdate"]),coef(fitval),residuo = residuo,residuo_mean = residuo_mean,log_like = likelihood,cov_matrix = cov_matrix, info=fitval$info)
  }
  save(results, file = "results.Rdata")
  print("Parallel ok")
  print(results)
  results <- as.data.frame(do.call(bind_rows, results))
  print(results)
  return(results)
}

run_all_fits <- function(outfolder='./'){
  all_fits <- list()
  i <- 1
  # all fits
  for (flag in c("srag", "covid")){
    for (diff in 0:1){
      # using both cases and deaths
      for (use in list(c('C', 'M'), c('M'))){
        for (logscale in c(T, F)){
          fname <- paste0(outfolder, flag, '.', diff, '.',
                          paste0(use, collapse=''), '.','log',
                          as.integer(logscale))
          png(paste0(fname, '.png'))
          s <- solver_wrapper(flag, diff=diff, use=use, logscale=logscale)
          all_fits[[i]] <- list(flag=flag, diff=diff, use=use,
                                logscale=logscale, result=s$result,
                                opt=s$opt, fitval=s$fitval)
          i <- i+1
          dev.off()
          saveRDS(s, paste0(fname, '.rds'))
        }
      }
    }
  }
  return(all_fits)
}

analyse_fit_results <- function(results){
  results <- results %>%
    group_by(startdate) %>%
    filter(residuo == min(residuo)) %>%
    arrange(startdate) %>%
    as.data.frame()
  results <- results %>%
    mutate(data = zoo::as.Date.numeric(startdate),
           # likelyhood = exp( - neg_log_like )
           # max(log_like) = - min(neg_log_like))
           # rel_like = like / max(like)
           #          = exp(- neg_log_like) / exp(min(neg_log_like))
           #          = exp(- neg_log_like - min(like))
           # ?? why - ?
           # ?? neg_log_like = - log_like (which is >0, since log_like is <0)
           # ?? (it is already rel) neg_like = 
           # smaller residual --> higher log-likelyhood
           rel_like = exp(-log_like + min(log_like)),
           # prob of each startdate
           prob = rel_like/sum(rel_like))
}

## !! Importante: antes de rodar run_the_fitting(), certificar-se que as intervenções
## são carregadas com as datas em relação à data de 01 de fevereiro:
## parameters["startdate"] <- as.numeric(as.Date("2020-02-01"))
## interventions <- update_scenario_parameters(scenario = scenario, parameters)

# Recommended to run these line by line by hand
run_the_fitting <- function(flag='srag', diff=1, logscale=FALSE, filename){
  ## fitting
  print("Running!")
  s <- solver_wrapper_fitval(flag, diff=diff, use=c('C', 'M'), logscale=logscale, rollwindow = 1)
  print("s solved")
  results <- analyse_fit_results(s)
  print("saving resuts")
  # save result
  if (missing(filename))
    filename <- paste0('fits/', flag, '.', diff, '.CM.log',
                       as.integer(logscale), '-', data.base, '.csv')
  write.csv(results, filename, row.names=F)
  print("Saved")
  # optimal parameter set
  opt<-results[which.max(results$prob),]
  
  # ## plots
  # # log-likelihood
  # x11()
  # plot(results$data, results$log_like, type='l', xlab='start date', ylab='log-likelihood')
  # # probability
  # x11()
  # plot(results$data, results$prob, type='l', xlab='start date', ylab='probability')
  # x11()
  # # optimal and 95% CI of p and perc for optimal startdate
  # confidence_interval_estimate(unlist(opt[7:10]), unlist(opt[2:3]), dof=opt$residuo/opt$residuo_mean - 2)
  # # plot of the best fit compared to the data
  # x11()
  # test_modifying_parameters(parameters, interventions, Y, unlist(opt[1:3]), flag='srag')
  # 
  ## sensitivity analysis
  # change the fit.table file in the beginning of the script to point to the
  # one you just saved and then source the whole file
  
  ## updating the data
  # change data.base in the beginning of dados_municipio_SP.R with the last
  # value available in the site
  
  ## ALTERNATIVE methods
  # fit using more params
  # params <- c(p = 0.02, perc_threshold = 0.5, scale_ihr = 1)
  # lower <- c(p = 0.01, perc_threshold = 0.1, scale_ihr = 0.7)
  # upper <- c(p = 0.04, perc_threshold = 0.9, scale_ihr = 1.3)
  # datej <- 10:25 # startdate, in days sicne 2020-02-01
  # s <- solver_wrapper_n(params, lower, upper, datej, n=10, flag, diff=diff, use=c('C', 'M'), logscale=logscale)
  print("Done.")
}

######################################LEGACY FUNCTIONS##############################################
solver_wrapper <- function(flag_ = "covid", diff=0, use=c('C', 'M'),
                           logscale=TRUE, rollwindow=0){
  params <- c(p = 0.015,perc_threshold = 0)
  lower <- c(p = 0, perc_threshold = 0)
  upper <- c(p = 0.2, perc_threshold = 0)
  data_start <- as.Date("2020-02-01")
  # run in parallel in linux
  if(Sys.info()[["sysname"]] == "Linux"){
    registerDoParallel(cores=detectCores()-2)
  }
  results <-  foreach(j = 0:30,.combine = rbind) %:% foreach(i = 2:3) %dopar% {
    params["p"] <- 0.010 + 0.005*i
    date_params <- c(startdate = as.Date(data_start + j))
    fitval <- nls.lm(par = params, fn = least_squares_solve, upper = upper,
                     lower = lower, date_params = date_params,flag = flag_,
                     diff=diff, use=use, logscale=logscale,
                     rollwindow=rollwindow)
    print(fitval)
    #print(date_params)
    c(j=j, coef(fitval), sol=fitval$deviance, p_init = params["p"])
  }
  results <- as.data.frame(do.call(bind_rows, results))
  print(results)
  min_res <- min(results$sol)
  fit_params <- c()
  for(i in names(params))
    fit_params[i] <- results[results$sol == min_res,i]
  fit_params['startdate'] <- results[results$sol == min_res,'startdate']
  j_ <- results[results$sol == min_res, "j"]
  # cov_matrix <- results[results$sol == min_res,5:8]
  summary <- fit_params
  summary["startdate"] <- as.Date(data_start + j_)
  print(summary)
  # print(cov_matrix)
  
  test_modifying_parameters(parameters, interventions, Y, fit_params, flag=flag_)
  params["p"] <- results[results$sol == min_res,"p_init.p"]
  date_params <- c(startdate = as.Date(data_start + j_))
  fitval <- nls.lm(par=params, fn=least_squares_solve, upper=upper,
                   lower=lower, date_params=date_params, flag=flag_, diff=diff,
                   use=use)
  #print(vcov(fitval))
  #print(mean(abs(fitval$fvec)))
  return(list(results=results, opt=summary, fitval=fitval))
}

solver_wrapper_n <- function(params, lower, upper, datej, n=5, flag_ = "covid",
                             diff=0, use=c('C', 'M'), logscale=TRUE,
                             rollwindow=0){
  q <- rep('qunif', length(params))
  q.args <- list()
  for (p in names(params))
    q.args[[p]] <- list(min=unname(lower[p]), max=unname(upper[p]))
  n <- max(n, length(params)+2)
  
  # use LHS to sample the parameter space
  uncoupledLHS <- LHS(model=NULL, factors=names(params), N=n,
                      q=q, q.arg=q.args, method='random')
  starters <- uncoupledLHS$data
  
  data_start <- as.Date("2020-02-01")
  # run in parallel in linux
  if(Sys.info()[["sysname"]] == "Linux"){
    registerDoParallel(cores=detectCores()-1)
  }
  results <-  foreach(j = datej, .combine = rbind) %:% foreach(i = 1:n) %dopar% {
    #params["p"] <- 0.010 + 0.005*i
    params <- unlist(starters[i,])
    date_params <- c(startdate = as.Date(data_start + j))
    # Insert update from interventions here
    fitval <- nls.lm(par = params, fn = least_squares_solve, upper = upper,
                     lower = lower, date_params = date_params,flag = flag_,
                     diff=diff, use=use, logscale=logscale,
                     rollwindow=rollwindow)
    residuo<- sum(residuals(fitval)*residuals(fitval))
    residuo_mean <- residuo/length(residuals(fitval))
    likelihood <- length(residuals(fitval))*log(residuo_mean)
    cov_matrix <- tryCatch({vcov(fitval)},error = function(cond){return(as.numeric(rep(NA, length(params)^2)))})
    
    c(startdate = as.numeric(date_params["startdate"]),coef(fitval),residuo = residuo,residuo_mean = residuo_mean,log_like = likelihood,cov_matrix = cov_matrix, info=fitval$info)
  }
  results <- as.data.frame(do.call(bind_rows, results))
  return(results)
}