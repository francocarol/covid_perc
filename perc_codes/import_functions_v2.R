if(!require(googlesheets4)){install.packages("googlesheets4"); library(googlesheets4)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
# gsheet ORIGINAL:
#doc_link <- "https://docs.google.com/spreadsheets/d/1Xad_Qo5eu_DwEe42OGmPRVzbpLG20QnXOgFgXdyt3-4"
# gsheet Carol v1:
#doc_link <- "https://docs.google.com/spreadsheets/d/121IOvCPyvs1AbzHAxXT3rGqnIz9lRr06ZRSDBpE5hsA/edit?usp=sharing"
# gsheet Carol v2:
# doc_link <- "https://docs.google.com/spreadsheets/d/16OoSxuxizrw07-vl1c5jbZolv3T_WCECC5E5fVJA1Ig/edit?usp=sharing"
# outra gsheet:
#doc_link <- "https://docs.google.com/spreadsheets/d/1sUYRzBBadVZbq8Qchr3i4uLpPlG9NlF6SlBI25B61go/edit?usp=sharing"
# Brazilian Model v2
# doc_link <- "https://docs.google.com/spreadsheets/d/1MG6vy0TVxnWxideH4RimlF12cWAfk4EZTxLdkcy_33o/edit?usp=sharing"

# Brazilian Model v2 Nova (parâmetros são copiados de csv's em DATA)
doc_link <- "https://docs.google.com/spreadsheets/d/1uRlMRlO1CwUV6gqceGSm7buCxEgYxhlD2FZ_cCr2pyI/edit?usp=sharing"

#In case you need to set the wd to the current folder:
#setwd("/home/caroline/Documents/modelling_Covid19BR")

## TODO: replace hard-coded indexes with named vector
# sheets:
# CSV OK - 2. DistrEtaria
# 4. matrix_cont_all_loc
# 5. matrix_cont_all_home
# 6. matrix_cont_all_school
# 7. matrix_cont_all_work
# 8. matrix_cont_other_loc
# 9. baseline_parameters
# 10. pop birth
# 11. pop mort
# 12. DistrEtaria 2018
# 13. DistrEtaria SP
# 14. ihr
# 15. ifr
# 16. Scenario 1
# 17. Scenario 2: SP state currently
# 18. Scenario 3: SP state plan 1: total official lockdown

if(!exists("estado", mode="character")) {estado = "SP"}
if(!exists("cidade", mode="character")) {cidade = "Sao_Paulo"}

# public access
sheets_deauth()

# age distribution (2020 by state)

DistrEtaria <-  data.frame(read.csv("DATA/2. DistrEtaria2020.csv")[-1,])
DistrEtariaBR <- data.frame(agefloor=unlist(DistrEtaria[,1]),
                            pop=rowSums(DistrEtaria[,-1]))
DistrEtaria18 <- data.frame(read.csv("DATA/12. DistrEtaria2018.csv")[-1,])
DistrEtariaBR18 <- data.frame(agefloor=unlist(DistrEtaria18[,1]),
                              pop=rowSums(DistrEtaria18[,-1]))
# ---> CORRIGIDO! ANTES USAVA DistrEtaria (2020)
DistrEtaria_local <- data.frame(read.csv(paste0("DATA/", estado, "/", cidade, "/13. DistrEtaria_cidade2020.csv")))

# birth/mort
pop.birth.br <- data.frame(read.csv("DATA/10. pop_birth.csv"))
pop.mort.br <- data.frame(read.csv("DATA/11. pop_mort.csv"))

#########    SEVERITY AND MORTALITY
# age dependent hosp and mort
#ihr and ifr and prob_icu

ihr <- data.frame(read.csv("DATA/14. ihr.csv"))
ihr[,2] <- as.numeric(ihr[,2])
ifr <- data.frame(read.csv(paste0("DATA/", estado, "/", cidade, "/15. ifr.csv")))
ifr[,3] <- as.numeric(ifr[,3])
ifr[,4] <- as.numeric(ifr[,4])
ihr[,2]<-ihr[,2]/100   # csv data is in percentages
#ifr_original<-ifr/100   # csv data is in percentages
## TODO: is normalizing and multiplying by pdeath the best option?
ifr[,3] <- ifr[,3]/max(ifr[,3])
ifr[,4] <- ifr[,4]/max(ifr[,4])
prob_icu <- data.frame(read.csv("DATA/16. prob_icu.csv"))
prob_icu <- as.numeric(prob_icu[,2])/100
pclin <- data.frame((read.csv("DATA/symptomatics_per_age.csv")))
pclin[,2] <- as.numeric(pclin[,2])
# basic parameters

update_parameters <- function(parameters, estado = "SP", cidade = "Sao_Paulo"){
  date_params <- c("startdate", "stopdate")
  params <- as.data.frame(read.csv(paste0("DATA/", estado, "/", cidade, "/9. baseline_parameters.csv"), stringsAsFactors=FALSE)[,1:2])
  for(i in rownames(params)){
    if (is.na(params[i, "param"]))
      next
    if (params[i, "param"] %in% date_params){
      parameters[params[i, "param"]] <- as.numeric(as.Date(unlist(params[i, "value"]), format="%Y-%m-%d"))
    }
    else
      parameters[params[i, "param"]] <- as.numeric(params[i, "value"])
  }
  return(parameters)
}

# scenario parameters
###################### CHECK NUMBERING FOR SCENARIOS, in NEW TABLE !!!!!!!!
###################### CHECK NUMBERING FOR SCENARIOS, in NEW TABLE !!!!!!!!
###################### CHECK NUMBERING FOR SCENARIOS, in NEW TABLE !!!!!!!!
update_scenario_parameters <- function(parameters, scenario=3, estado = NA, cidade = NA){
  
  if(is.na(estado)|is.na(cidade)) {
  scenario.param <- data.frame(range_read(doc_link, scenario))[,1:5] } else {
    scenario.param <- read.csv(paste0("DATA/", estado, "/", cidade, "/base_scenario.csv"))
  }
  
  ##hard coded intervention names, but easier than the main code##
  intervention_names <- c("Self Isolation","Social Distancing","Handwashing","School Closing","Work From Home",
                          "Vaccination","Travel Ban","Cocoon Elderly","Cocoon Age","Screening","Quarantine","Quarantine Days",
                          "Window Time","Prob Test S","Prob Test E","Prob Test I",
                          "Prob Test CL","Prob Test HC","Prob Test ICUC","Prob Test X","Prob Test ICUH",
                          "Prob Test R","Prob Test ICU","Prob Test H","Home Tracing Eff","School Tracing Eff",
                          "Work Tracing Eff","Other Tracing Eff","School Closing 1","School Closing 2","School Closing 3","School Closing 4",
                          "Tests Per Day")
  interventions = list()
  j = 1
  for(i in intervention_names){
    f <- scenario.param[scenario.param$Intervention == i,]
    if(nrow(f) != 0){
      startdate <- c(as.numeric(as.Date(f$Start.date, format="%Y-%m-%d"))-parameters["startdate"])
      enddate <- c(as.numeric(as.Date(f$End.date, format="%Y-%m-%d"))-parameters["startdate"])
      value1 <- c(as.numeric(f$Value1))
      value2 <- c(as.numeric(f$Value2))
      f <- rbind(startdate,enddate,value1,value2)
      f["startdate", f["startdate",] < 0] <- 0
      colnames(f) <- NULL
    }
    else {
      startdate <-NA
      enddate <- NA
      value1 <- NA
      value2 <- NA
      f <- rbind(startdate,enddate,value1,value2)
      colnames(f) <- NULL
    }
    interventions[[j]] <- f
    j = j+1
  }
  names(interventions) <-  intervention_names
  return(interventions)
}

####### This function overwrite intervention parameters present in the spreadsheet
####### Other interventions should be previously loaded using update_scenarios_parameters() from a "base scenario"

reload_scenario_parameters <- function(parameters, scenario, interventions){
  scenario.param <- data.frame(range_read(doc_link, scenario))[,1:5]
  scenario.param <- scenario.param[complete.cases(scenario.param),]
  
  ##hard coded intervention names, but easier than the main code##
  intervention_names <- c("Self Isolation","Social Distancing","Handwashing","School Closing","Work From Home",
                          "Vaccination","Travel Ban","Cocoon Elderly","Cocoon Age","Screening","Quarantine","Quarantine Days",
                          "Window Time","Prob Test S","Prob Test E","Prob Test I",
                          "Prob Test CL","Prob Test HC","Prob Test ICUC","Prob Test X","Prob Test ICUH",
                          "Prob Test R","Prob Test ICU","Prob Test H","Home Tracing Eff","School Tracing Eff",
                          "Work Tracing Eff","Other Tracing Eff","School Closing 1","School Closing 2","School Closing 3","School Closing 4",
                          "Tests Per Day")
  
  scenarios_to_update =  unique(scenario.param$Intervention)[unique(scenario.param$Intervention) %in% intervention_names]
  
  contain.na <- function(x) return(any(grepl("NA", x)))
  for(i in scenarios_to_update){
    f <- scenario.param[scenario.param$Intervention == i,]
    if(nrow(f) != 0){
      if(!any(apply(f, 1, contain.na))) {
        startdate <- c(as.numeric(as.Date(f$Start.date, format="%Y-%m-%d"))-parameters["startdate"])
        enddate <- c(as.numeric(as.Date(f$End.date, format="%Y-%m-%d"))-parameters["startdate"])
        value1 <- c(as.numeric(f$Value1))
        value2 <- c(as.numeric(f$Value2))
        f <- rbind(startdate,enddate,value1,value2)
        f["startdate", f["startdate",] < 0] <- 0
        colnames(f) <- NULL } else {
          startdate = NA; enddate = NA; value1 = NA; value2 = NA
          f <- rbind(startdate,enddate,value1,value2)
          colnames(f) <- NULL
        }
    }
    j = which(names(interventions) == i)
    interventions[[j]] <- f
  }
  
  return(interventions)
}

switch_off_interventions <- function(parameters){
  #TODO: should update this?
  parameters["selfis_on"]<-10e5
  parameters["dist_on"]<-10e5
  parameters["hand_on"]<-10e5
  parameters["work_on"]<-10e5
  parameters["school_on"]<-10e5
  parameters["cocoon_on"]<-10e5
  parameters["vaccine_on"]<-10e5
  return(parameters)
}

import_sensitivity_parameters_distr <- function(parameters, estado = "SP", cidade = "Sao_Paulo"){
  # basic params
  basic.params <- as.data.frame(read.csv(paste0("DATA/", estado, "/", cidade, "/9. baseline_parameters.csv"), stringsAsFactors=FALSE)[,1:2])
  qs = list()
  q.args = list()
  for(i in rownames(basic.params)){
    if (length(pos <- which(parameters == basic.params[i, "param"]))){
      qs[[pos]] <- paste0("q", basic.params[i, "distribution"])
      # TODO: fix this to accept param number lower or higher than 2
      q.args[[pos]] <- list(as.numeric(basic.params[i, "param1"]),
                            as.numeric(basic.params[i, "param2"]))
      if (qs[[pos]] == 'qtnorm'){
        # truncated normal: lower value of zero
        q.args[[pos]][[3]] <- 0
      }
    }
  }
  # not formatted correctly yet
  #scenario.params <- data.frame(sheets_read(doc_link, scenario))
  #for (i in 1:9){
  #    if (length(pos = which(parameters == scenario.params[i, "param"]))){
  #        qs[[pos]] <- basic.params[i, "distribution"]
  #        q.args[[pos]] <- list(basic.params[i, "param1"], basic.params[i, "param1"])
  #    }
  
  #    parameters[scenario.param[i,1]] <-as.numeric(as.Date(scenario.param[[i,2]],format="%Y-%m-%d")-startdate)
  #}
  return(list(unlist(qs), q.args))
}

get_latest_fit <- function(flag = 'srag', diff = 1, logscale = FALSE,
                           fits_dir = 'fits/', estado = "SP", model_version= "model1"){
  files <- list.files(paste0(fits_dir, estado, "/"), pattern = paste0(model_version, '_', flag, '.', diff, '.CM.log',
                                                 as.integer(logscale),
                                                 '-.*csv'))
  datas <- str_extract(files, '\\d\\d\\d\\d_\\d\\d_\\d\\d')
  lastfile <- files[which(datas == max(datas, na.rm=T))]
  dir_lastfile <- paste0(estado, "/", lastfile)
}

update_params_from_fit <- function(parameters, file){
  d <- read.csv(file)
  dopt <- d[d$residuo == min(d$residuo),]
  # All parameters that match the ones at the fit result file will be replaced 
  pars_to_update <- names(parameters)[names(parameters) %in% names(dopt)]
  parameters[pars_to_update] <- unlist(dopt[pars_to_update])
  # parameters['p'] = dopt$p
  # parameters['perc_threshold'] = dopt$perc_threshold
  # parameters['startdate'] = dopt$startdate
  return(parameters)
}
