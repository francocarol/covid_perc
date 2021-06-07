if(!require(parallel)){install.packages("parallel"); library(parallel)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(zeallot)){install.packages("zeallot"); library(zeallot)} # provides %<-%/
if(!require(MASS)){install.packages("MASS"); library(MASS)}
if(!require(pse)){install.packages("pse"); library(pse)}
if(!require(extraDistr)){install.packages("extraDistr"); library(extraDistr)} # discrete uniform and truncated normal
if(!require(zoo)){install.packages("zoo"); library(zoo)}
if(!require(purrr)){install.packages("purrr"); library(purrr)} # partial
if(!require(abc)){install.packages("abc"); library(abc)}

# source before running this:
#source('Brazilian_Model_v1b.R')

default_parameters <- function(scenario){
  # base parameters and intervention
  date_parameters <- c(startdate=as.Date('2020-1-1'))
  parameters <- c()
  
  ########### take parameters from external tables ###########
  newps <- update_parameters(parameters, date_parameters)
  parameters <- newps[[1]]
  date_parameters <- newps[[2]]
  newps <- update_scenario_parameters(parameters, date_parameters,
                                      scenario=scenario)
  parameters <- newps[[1]]
  date_parameters <- newps[[2]]
  
  return(list(parameters, date_parameters))
}

samples_from_prior <- function(sensitive.params, N){
  # set up parameter's distributions
  distrs <- import_sensitivity_parameters_distr(sensitive.params)
  qs <- distrs[[1]]
  q.args <- distrs[[2]]
  
  uncoupledLHS <- LHS(model=NULL, factors=sensitive.params, N=N, q=qs,
                      q.arg=q.args, method='random')
  
  return(uncoupledLHS$data)
}

# sample from posterior: use dplyr:
# sample_n(df, N, replace=T)

sample_p_threshold <- function(v_matrix, parms, dof){
  # dof: degrees of freedom = # data points - # params
  S <- matrix(as.numeric(v_matrix), nrow=length(parms))
  Sinv <- ginv(S)
  r <- sqrt(rf(1, length(parms), dof)*2)
  theta <- runif(1, 0, 2*pi)
  # transform points of circle into points of ellipse using
  # svd of inverse covariance matrix
  Sinv_svd <- svd(Sinv) # inverse of covariance matrix
  # transform from circle to ellispse
  xt <- t(Sinv_svd$v)%*%diag(1/sqrt(Sinv_svd$d))%*%c(r*cos(theta), r*sin(theta))
  # translate the ellipse so that center is the estimated parameter value
  xt <- xt + parms
}

# fit.table <- read.csv('fits_srag.csv')
sample_parameter <- function(fit.table, size, params=c('p', 'perc_threshold')){
  rows <- sample(rownames(fit.table), size=size, prob=fit.table$prob,
                 replace=T)
  startdates <- fit.table[rows, 'startdate']
  pars <- fit.table[rows, c(params, paste0('cov_matrix', 1:(length(params)^2)))]
  dof <- fit.table[1,]$residuo/fit.table[1,]$residuo_mean - length(pars)
  
  # sample params
  sampled.pars <- apply(pars, 1, function(p){sample_p_threshold(
    v_matrix=p[paste0('cov_matrix', 1:(length(params)^2))],
    parms=p[params], dof=dof)})
  sampled.pars <- as.data.frame(t(sampled.pars))
  colnames(sampled.pars) <- params
  sampled.pars$startdate <- as.Date.numeric(startdates)
  return(sampled.pars)
}

result_ts <- function(out, parameters, date_parameters,
                      cols=list(C=Cindex+1, CM=CMindex+1), diff=TRUE){
  results <- list()
  for (v in names(cols)){
    results[[v]] <- rowSums(out[,cols[[v]]])
  }
  z <- zoo(bind_cols(results), as.Date(date_parameters["startdate"]) + out[,1])
  if (diff)
    return(diff(z))
  else
    return(z)
}

bind_ts_quantiles <- function(sim.out, quantiles=c(0.025, 0.5, 0.975)){
  QQ <- list()
  for (coln in colnames(sim.out[[1]])){
    listcol <- lapply(sim.out, function(x){x[,coln]})
    listcol[['fill']] <- 0.
    zoo.col <- do.call(merge.zoo, listcol)
    zoo.quantiles <- list()
    for (qi in quantiles)
      zoo.quantiles[[as.character(qi)]] <- 
      zoo(apply(zoo.col, 1, partial(quantile, prob=qi)),
          time(zoo.col))
    zoo.quantiles <- do.call(cbind.zoo, zoo.quantiles)
    colnames(zoo.quantiles) <- paste0(coln, '.', 100*quantiles)
    QQ[[coln]] <- zoo.quantiles
  }
  return(do.call(cbind.zoo, QQ))
}

# run simulations in parallel
run_simulations <- function(parameters, interventions, params.df, output.f,
                            post.process, fname, numCores){
  sensitive.params <- colnames(params.df)
  print(params.df)
  generate_modelRun <- function(sensitive.params, parameters,
                                interventions){
    f <- function(...){
      arg <- list(...)
      pars <- as.numeric(unlist(arg))
      print(pars[1])
      for (i in 1:length(sensitive.params)){
        parameters[sensitive.params[i]] <- pars[i]
        print(parameters[sensitive.params[i]])
      }
      if ('shift_start' %in% sensitive.params){
        # update all date parameters
        # only startdate is actually used in the function, the others are relative to it
        parameters['startdate'] <- parameters['startdate'] + parameters['shift_start']
        # pnames <- names(date_parameters)
        # pnames_dates <- pnames[substr(pnames, 1,5) == 'date_']
        # datepars <-  substr(pnames_dates, 6, nchar(pnames_dates))
        # some black magic with zeallot
        # notice the minus sign: startdate advances, relative dates are pushed back
        for(name in names(interventions)){
          interventions[[name]][,1]["startdate"] <- interventions[[name]][,1]["startdate"]-parameters['shift_start']
          interventions[[name]][,1]["enddate"] <- interventions[[name]][,1]["enddate"]-parameters['shift_start']
        }   
        
        # parameters[datepars] %<-% as.list(parameters[datepars] - parameters['shift_start'])
        parameters <- unlist(parameters)
      }
      # 
    print(parameters)
      # print(interventions)
      
      out <- run_covid19(parameters, interventions, Y)
      
      return(output.f(out[[2]], parameters, interventions))
    }
  }
  modelRun <- generate_modelRun(sensitive.params, parameters,
                                interventions)
  if(Sys.info()[["sysname"]] == "Linux"){
  if (missing(numCores))
    numCores <- detectCores() - 1
  sim.out <- mclapply(as.data.frame(t(params.df)), modelRun,
                      mc.cores=numCores)
  } else {
    sim.out <- lapply(as.data.frame(t(params.df)), modelRun)
  }
  if (! missing(post.process))
    sim.out <- post.process(sim.out)
  # save results to file
  # be careful to not overwrite working tables with broken ones!
  if (! missing(fname))
    write.csv(sim.out, file=fname, row.names=FALSE)
  
  return(sim.out)
}

### summary statistics
# cumulative incidence between days
ssSum <- function(zoo.obj) {
  return(sum(zoo.obj))
}

# growth rate during period
## copied from funcoes.R
## TODO: translate documentation
#' Ajusta um modelo exponencial para n de casos em função de n de dias
#' passados desde o início da série.
#' @details Ajusta um modelo generalizado misto para contagens
#'     (Poisson) dos valores de uma série temporal em função do número
#'     de dias transcorridos desde o início da série. Como o modelo
#'     Poisson tem função de ligação log, equivale a ajustar uma
#'     função de crescimento exponencial às contagens, mas com erro
#'     não gaussianos e sim Poisson.
ssGrowth <- function(zoo.obj){
  ndias <- as.vector(rev(max(time(zoo.obj)) - time(zoo.obj)))
  fit <- glm(zoo.obj ~ ndias, family = poisson)
  return(unname(coef(fit)[2]))
}

summary_stats <- function(zoo.obj, cases.start, cases.end,
                          deaths.start, deaths.end){
  cases <- window(zoo.obj[,1], start=cases.start, end=cases.end)
  deaths <- window(zoo.obj[,2], start=deaths.start, end=deaths.end)
  return(c(sum.cases=ssSum(cases),
           growth.cases=ssGrowth(cases),
           sum.deaths=ssSum(deaths),
           growth.deaths=ssGrowth(deaths)))
}

abc_filter <- function(params.df, sim.sstats, target, tol=0.05){
  test = abc(target = target, # observed summary statistics
             param = data.frame(params.df),
             sumstat = data.frame(sim.sstats),
             tol = tol, # controls the deviation of the estimates from the observed; large deviations are excluded
             method = 'rejection') # classic rejection method; simple, a bit restrictive
  
  print(summary(test))
  
  return(test)
}

