###########################################################################
## AGE-DEPENDANT SEIRS MODEL WITH 5-YEAR AGE CLASSES USING UN DEMOG DATA ##
###########################################################################
if(!require(deSolve)){install.packages("deSolve"); library(deSolve)}
if(!require(zeallot)){install.packages("zeallot"); library(zeallot)}

source("etaclass.R")
source("import_functions_v2.R")

#########  INCIDENCE DATA
#incdata_X<-read.csv("NYCcovidcases3.csv")

#########  DEMOGRAPHIC DATA
# population structure in 2020
popstruc <- DistrEtaria_local 
popstruc18 <- DistrEtariaBR18

A<-length(popstruc[,2])
initP<-sum(popstruc[,2]) # population size 
ageindcase<-20 # age of index case (years)
aci <- floor((ageindcase/5)+1) # age class of index case

# births by age of mother
popbirth <- pop.birth.br$Brasil
# convert from per year period to per person per day
popbirth<-popbirth/(popstruc18[,2]*365.25)

#natural mortality per person per year
popmort<-pop.mort.br
# convert from 1000s per year period to per person per day
mort<-popmort[,2]/(popstruc18[,2]*365.25)

##########   CONTACT DATA
# c_home <- as.matrix(sheets_read(doc_link, 5))
# c_school <- as.matrix(sheets_read(doc_link, 6))
# c_work <- as.matrix(sheets_read(doc_link, 7))
# c_other <- as.matrix(sheets_read(doc_link, 8))

# c_home <- data.frame(read.csv("DATA/5. matrix_cont_home.csv", stringsAsFactors = FALSE, dec=","))
# write.csv(c_home,"DATA/5. matrix_cont_home_PONTO.csv", row.names=FALSE)
# c_school <- data.frame(read.csv("DATA/6. matrix_cont_school.csv", stringsAsFactors = FALSE, dec=","))
# write.csv(c_school,"DATA/6. matrix_cont_school_PONTO.csv", row.names=FALSE)
# c_work <- data.frame(read.csv("DATA/7. matrix_cont_work.csv", stringsAsFactors = FALSE, dec=","))
# write.csv(c_work,"DATA/7. matrix_cont_work_PONTO.csv", row.names=FALSE)
# c_other <- data.frame(read.csv("DATA/8. matrix_cont_other_loc.csv", stringsAsFactors = FALSE, dec=","))
# write.csv(c_other,"DATA/8. matrix_cont_other_loc_PONTO.csv", row.names=FALSE)
Prem2020 <- TRUE

if(!exists("Prem2020")) {
  Prem2020 <- FALSE
  #outmsg <- paste("Prem2020 is", Prem2020)
  #print(outmsg)
}

if(!Prem2020) {
  print("Contact Matrix from Prem et al. 2017: former data")
  c_home <- data.frame(read.csv("DATA/5. matrix_cont_home_PONTO.csv", stringsAsFactors = FALSE))
  c_school <- data.frame(read.csv("DATA/6. matrix_cont_school_PONTO.csv", stringsAsFactors = FALSE))
  c_work <- data.frame(read.csv("DATA/7. matrix_cont_work_PONTO.csv", stringsAsFactors = FALSE))
  c_other <- data.frame(read.csv("DATA/8. matrix_cont_other_loc_PONTO.csv", stringsAsFactors = FALSE))
} else {
  print("Contact Matrix from Prem et al. 2020: New data")
  load("PREM/prem2020_bra.Rdata")
  c_home <- u_home
  c_school <- u_school
  c_work <- u_work
  c_other <- u_others
}

nce <-A-length(c_home[1,])

contact_home<-matrix(0,nrow=A,ncol=A)
contact_school<-matrix(0,nrow=A,ncol=A)
contact_work<-matrix(0,nrow=A,ncol=A)
contact_other<-matrix(0,nrow=A,ncol=A)

for (i in 1:(A-nce)){
  for (j in 1:(A-nce)){
    contact_home[i,j]<-c_home[i,j]
    contact_school[i,j]<-c_school[i,j]
    contact_work[i,j]<-c_work[i,j]
    contact_other[i,j]<-c_other[i,j]
  }
}

for (i in (A+1-nce):A){
  for (j in 1:(A-nce)){
    contact_home[i,j]<-c_home[(A-nce),j]
    contact_school[i,j]<-c_school[(A-nce),j]
    contact_work[i,j]<-c_work[(A-nce),j]
    contact_other[i,j]<-c_other[(A-nce),j]
  }
}
for (i in 1:(A-nce)){
  for (j in (A+1-nce):A){
    contact_home[i,j]<-c_home[i,(A-nce)]
    contact_school[i,j]<-c_school[i,(A-nce)]
    contact_work[i,j]<-c_work[i,(A-nce)]
    contact_other[i,j]<-c_other[i,(A-nce)]
  }
}
for (i in (A+1-nce):A){
  for (j in (A+1-nce):A){
    contact_home[i,j]<-c_home[(A-nce),(A-nce)]
    contact_school[i,j]<-c_school[(A-nce),(A-nce)]
    contact_work[i,j]<-c_work[(A-nce),(A-nce)]
    contact_other[i,j]<-c_other[(A-nce),(A-nce)]
  }
}


# average contacts per day from POLYMOD matrices 
avg_cont<-sum((contact_home+contact_other+contact_school+contact_work)%*%(popstruc[,2]/sum(popstruc[,2])))

#########    POP AGEING MATRIX ############################################
# per year ageing matrix
dd<-seq(1:A)/seq(1:A)
ageing <- t(diff(diag(dd),lag=1)/(5*365.25))
ageing<-cbind(ageing,0*seq(1:A)) # no ageing from last compartment

###########################################################################
# Define the indices for each variable and initial conditions
variables <- c('S', 'E', 'I', 'R', 'X',
               'H', 'HC', 'C', 'CM', 'V',
               'QS', 'QE', 'QI', 'QR', 'CL', 'QC',
               'ICU', 'ICUH', 'ICUC', 'CMC')
for (v in 1:length(variables)){
  eval(parse(text=paste0(variables[v], 'index <- seq(', 1+A*(v-1), ',', v*A, ')')))
  eval(parse(text=paste0('init', variables[v], ' <- 0*popstruc[,2]')))
}
initE[aci]<-1          # place random index case in E compartment
initS <- popstruc[,2] - (eval(parse(text=paste0(paste(paste0('init', variables[-1]), collapse='+')))))
Y <- eval(parse(text=paste0('c(', paste(paste0('init', variables), collapse=','), ')')))
######################################################################################################
vectorize_interventions <- function(parameters, interventions){
  ###This function converts every intervention to vectors with the length of the simulation,
  ###Sadly, every intervention must be hard coded since the parameter names are different,
  ###but the hard work has already done
  ######Self Isolation#########
  f = interventions$`Self Isolation`
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  isolation <- data.frame(vec, cov, eff)
  
  ######Social Distancing#########
  f = interventions$`Social Distancing`
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  distancing <- data.frame(vec, cov, eff)
  #########HandWashing###########
  f = interventions$Handwashing
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  handwashing <- data.frame(vec, cov, eff)
  
  #########School Closing###########Uniform closing##################
  f = interventions$`School Closing`
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  ##################################################################
  f = interventions$`School Closing 1`
  vec1 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov1 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff1 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec1[j] <- TRUE
        cov1[j] <- cov_
        eff1[j] <- eff_
      }
    }
  }
  ###################################################################
  f = interventions$`School Closing 2`
  vec2 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov2 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff2 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec2[j] <- TRUE
        cov2[j] <- cov_
        eff2[j] <- eff_
      }
    }
  }
  ##################################################################
  f = interventions$`School Closing 3`
  vec3 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov3 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff3 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec3[j] <- TRUE
        cov3[j] <- cov_
        eff3[j] <- eff_
      }
    }
  }
  f = interventions$`School Closing 4`
  vec4 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov4 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff4 <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec4[j] <- TRUE
        cov4[j] <- cov_
        eff4[j] <- eff_
      }
    }
  }
  schoolclosing <- data.frame(vec, cov, eff,vec1, cov1, eff1,vec2, cov2, eff2,vec3, cov3, eff3,vec4, cov4, eff4)
  
  #########Work From Home###########
  f = interventions$`Work From Home`
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  workhome <- data.frame(vec, cov, eff)
  #########Vaccination###########
  f = interventions$Vaccination
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  vaccination <- data.frame(vec, cov, eff)
  
  #########Travel Ban###########
  f = interventions$`Travel Ban`
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  mean_imports <- cov
  travelban <- data.frame(vec, mean_imports, eff)
  
  #########Cocoon Elderly###########
  f = interventions$`Cocoon Elderly`
  f2 = interventions$`Cocoon Age`
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cocoon_age <- c(rep(1,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  for( i in 1:ncol(f2)){
    if(!is.na(f2["startdate",i])){
      cocoon_age_ <- floor(as.numeric(f2["value1",i])/5 + 1)
      for(j in seq(f2["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f2["enddate",i])+1)){
        cocoon_age[j] <- cocoon_age_
      }
    }
  }
  cocoon <- data.frame(vec, cov, eff,cocoon_age)
  
  #########Screening###########
  f = interventions$Screening
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
      }
    }
  }
  screening <- data.frame(vec,cov)
  
  #########Quarantine###########
  f = interventions$Quarantine
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  cov <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      cov_ <-  as.numeric(f["value1",i])
      eff_ <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        cov[j] <- cov_
        eff[j] <- eff_
      }
    }
  }
  #########Quarantine Days###########
  f = interventions$`Quarantine Days`
  days <- c(rep(1,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      days_ <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        days[j] <- days_
      }
    }
  }
  f = interventions$`Window Time`
  window_time <- c(rep(1,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      time_ <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        window_time[j] <- time_
      }
    }
  }
  quarantine <- data.frame(vec,cov,eff,days, window_time)
  
  ########Testing#######
  f = interventions$`Tests Per Day`
  vec <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  tests_day <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  asymp_eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      tests <- as.numeric(f["value1",i])
      asymp <- as.numeric(f["value2",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        vec[j] <- TRUE
        tests_day[j] <- tests
        asymp_eff[j] <- asymp
      }
    }
  }
  f = interventions$`Prob Test S`
  test_s <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_s[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test E`
  test_e <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_e[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test I`
  test_i <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_i[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test CL`
  test_cl <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_cl[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test HC`
  test_hc <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_hc[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test ICUC`
  test_icuc <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_icuc[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test ICUH`
  test_icuh <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_icuh[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test R`
  test_r <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_r[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test H`
  test_h <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_h[j] <- prob
      }
    }
  }
  f = interventions$`Prob Test ICU`
  test_icu <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_icu[j] <- prob
      }
    }
  }
  
  f = interventions$`Prob Test H`
  test_h <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      prob <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        test_h[j] <- prob
      }
    }
  }
  testing = data.frame(vec,tests_day,asymp_eff,test_s,test_e,test_i,test_cl,test_icu,test_icuh,test_r,test_h,test_icuc,test_hc)
  
  f = interventions$`Home Tracing Eff`
  home_tracing_eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      tracing <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        home_tracing_eff[j] <- tracing
      }
    }
  }
  
  f = interventions$`School Tracing Eff`
  school_tracing_eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      tracing <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        school_tracing_eff[j] <- tracing
      }
    }
  }
  
  f = interventions$`Work Tracing Eff`
  work_tracing_eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      tracing <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        work_tracing_eff[j] <- tracing
      }
    }
  }
  
  f = interventions$`Other Tracing Eff`
  other_tracing_eff <- c(rep(0,parameters["stopdate"]-parameters["startdate"] + 1))
  for( i in 1:ncol(f)){
    if(!is.na(f["startdate",i]) & (f["startdate",i] <= parameters["stopdate"])){
      tracing <- as.numeric(f["value1",i])
      for(j in seq(f["startdate",i]+1,min(parameters["stopdate"]-parameters["startdate"],f["enddate",i])+1)){
        other_tracing_eff[j] <- tracing
      }
    }
  }
  tracing <- data.frame(home_tracing_eff,school_tracing_eff,work_tracing_eff,other_tracing_eff)
  
  return(list(isolation = isolation,distancing = distancing,handwashing = handwashing,schoolclosing = schoolclosing, workhome = workhome,
              vaccination = vaccination, travelban = travelban,cocoon = cocoon,screening = screening, quarantine = quarantine,testing = testing,tracing = tracing))
}

# interpola valores linearmente dentro de cada dia
interpola <- function(vec, t){
  t_ = floor(t)
  if (t_ > 1) {
    return((t - t_) * vec[t_] + (1-(t-t_)) * vec[t_+1])
  } else {
    return(vec[t_+1])
  }
}
 
run_covid19 <- function(parameters, interventions, Y){
  #########   INITIALISE SIMULATION/INTERVENTION START TIMES
  day_start <- 0
  day_stop <- as.numeric(parameters['stopdate']-parameters['startdate'])
  #print(parameters)
  times <- seq(day_start, day_stop)
  
  tin<-as.numeric(parameters['startdate']-as.numeric(as.Date("2020-01-01")))/365.25
  
  # Scale parameters to percentages/ rates
  parameters["rho"]<-parameters["rho"]/100
  parameters["omega"]<-(1/(parameters["omega"]*365))
  parameters["gamma"]<-1/parameters["gamma"]
  parameters["nui"]<-1/parameters["nui"]
  parameters["report"]<-parameters["report"]/100
  parameters["reportc"]<-parameters["reportc"]/100
  parameters["reporth"]<-parameters["reporth"]/100
  parameters["nus"]<-1/parameters["nus"]
  parameters["rhos"]<-parameters["rhos"]/100
  parameters["amp"]<-parameters["amp"]/100
  parameters["selfis_dur"]<-parameters["selfis_dur"]*7
  parameters["selfis_cov"]<-parameters["selfis_cov"]/100
  parameters["selfis_eff"]<-parameters["selfis_eff"]/100
  parameters["dist_dur"]<-parameters["dist_dur"]*7
  parameters["dist_cov"]<-parameters["dist_cov"]/100
  parameters["dist_eff"]<-parameters["dist_eff"]/100
  parameters["hand_dur"]<-parameters["hand_dur"]*7
  parameters["hand_eff"]<-parameters["hand_eff"]/100
  parameters["work_dur"]<-parameters["work_dur"]*7
  parameters["work_cov"]<-parameters["work_cov"]/100
  parameters["work_eff"]<-parameters["work_eff"]/100
  parameters["w2h"]<-parameters["w2h"]/100
  parameters["school_dur"]<-parameters["school_dur"]*7
  parameters["school_eff"]<-parameters["school_eff"]/100
  parameters["s2h"]<-parameters["s2h"]/100
  parameters["cocoon_dur"]<-parameters["cocoon_dur"]*7
  parameters["cocoon_cov"]<-parameters["cocoon_cov"]/100
  parameters["cocoon_eff"]<-parameters["cocoon_eff"]/100
  parameters["age_cocoon"]<-floor((parameters["age_cocoon"]/5)+1)
  parameters["travelban_eff"]<-parameters["travelban_eff"]/100
  parameters["vaccine_eff"]<-parameters["vaccine_eff"]/100
  parameters["vaccine_cov"]<-parameters["vaccine_cov"]/100
  parameters["vac_campaign"]<-parameters["vac_campaign"]*7
  parameters["travelban_dur"]<-parameters["travelban_dur"]*7
  parameters["screen_dur"]<-parameters["screen_dur"]*7
  parameters["screen_cov"]<-parameters["screen_cov"]/100
  parameters["quarantine_cov"]<-parameters["quarantine_cov"]/100
  parameters["quarantine_days"]<-parameters["quarantine_days"]
  #parameters["ratem"]<-1/parameters["ratem"]
  #parameters["ratemHC"]<-1/parameters["ratemHC"]
  #parameters["ratemICU"]<-1/parameters["ratemICU"]
  #parameters["ratemICUC"]<-1/parameters["ratemICUC"]
  #parameters["ratemVent"]<-1/parameters["ratemVent"]
  #parameters["ratemVentC"]<-1/parameters["ratemVentC"]
  parameters["give"]<-parameters["give"]/100
  parameters["pdeath_h"]<-parameters["pdeath_h"]/100
  parameters["pdeath_hc"]<-parameters["pdeath_hc"]/100
  parameters["pdeath_icu"]<-parameters["pdeath_icu"]/100
  parameters["pdeath_icuh"]<-parameters["pdeath_icuh"]/100
  parameters["pdeath_icuc"]<-parameters["pdeath_icuc"]/100
  #    parameters["pdeath_vent"]<-parameters["pdeath_vent"]/100
  #    parameters["pdeath_venticu"]<-parameters["pdeath_venticu"]/100
  #    parameters["pdeath_venth"]<-parameters["pdeath_venth"]/100
  #    parameters["pdeath_ventc"]<-parameters["pdeath_ventc"]/100
  parameters["nusc"]<-1/parameters["nusc"]
  parameters["nu_icu"]<-1/parameters["nu_icu"]
  parameters["nu_icuh"]<-1/parameters["nu_icuh"]
  parameters["nu_icuc"]<-1/parameters["nu_icuc"]
  parameters["nu_vent"]<-1/parameters["nu_vent"]
  # parameters["nu_venticu"]<-1/parameters["nu_venticu"]
  # parameters["nu_venth"]<-1/parameters["nu_venth"]
  # parameters["nu_ventc"]<-1/parameters["nu_ventc"]
  parameters["prob_vent"]<-parameters["prob_vent"]/100
  #parameters["beds_available"] = sum(popstruc[,2])*parameters["beds"]/1000 # maximum number of hospital beds - numeric 
  #parameters["icu_beds_available"] = sum(popstruc[,2])*parameters["ICU_beds"]/ #8000, # maximum number of hospital beds - numeric
  
  
  vector_interventions <- vectorize_interventions(parameters, interventions)
  params <- list(parameters,vector_interventions)
  
  
  # set up a function to solve the equations
  covid <- function(t, Y, params){
    parameters <- params[[1]]
    interventions <- params[[2]]
    with(as.list(c(Y, parameters)),
         {
           S <- Y[Sindex] #susceptibles
           E <- Y[Eindex] #exposed
           I <- Y[Iindex] #assymptomatic cases
           R <- Y[Rindex] #recovered cases
           X <- Y[Xindex] #symptomatic mild cases in self-isolation
           H <- Y[Hindex] #moderate cases in regular beds
           HC <- Y[HCindex] #moderate cases without beds
           C <- Y[Cindex]
           CM <- Y[CMindex]
           V <- Y[Vindex]
           QS <- Y[QSindex]
           QE <- Y[QEindex]
           QI <- Y[QIindex]
           QR <- Y[QRindex]
           CL <- Y[CLindex]
           QC <- Y[QCindex]
           ICU <- Y[ICUindex] #severe cases with icu beds
           ICUH <- Y[ICUHindex] #severe cases in regular beds
           ICUC <- Y[ICUCindex] #severe cases with no beds at all
           CMC <- Y[CMCindex]
           
           P <- (S+E+I+R+X+H+HC+V+QS+QE+QI+QR+CL+QC+ICU+ICUH+ICUC)
           # print(sum(QS+QE+QI+QR))
           
           # health system performance
           critH <- etaclass(sum(H, ICUH), beds_available, give)
           critICU <- etaclass(sum(ICU), icu_beds_available, give)
           
           # TODO: use parameter instead of hard-coded 2
           critICUH <- critH^2
           
           ###########get the intervention values from the vectors############
           t_ <- floor(t)+1

           isolation<-interventions$isolation$vec[t_]
           distancing<-interventions$distancing$vec[t_]
           handwash<-interventions$handwashing$vec[t_]
           workhome<-interventions$workhome$vec[t_]
           ######################################################
           schoolclose_all<-interventions$schoolclosing$vec[t_]
           schoolclose_1<-interventions$schoolclosing$vec1[t_]
           schoolclose_2<-interventions$schoolclosing$vec2[t_]
           schoolclose_3<-interventions$schoolclosing$vec3[t_]
           schoolclose_4<-interventions$schoolclosing$vec4[t_]
           ######################################################
           cocoon<-(interventions$cocoon$vec[t_])
           vaccine<-interventions$vaccination$vec[t_]
           travelban<-interventions$travelban$vec[t_]
           screen<-interventions$screening$vec[t_]
           ######
           selfis_cov <- interventions$isolation$cov[t_]
           selfis_eff <- interventions$isolation$eff[t_]
           if (isolation){
             # if(screen){selfis<-min(selfis_cov/(1-screen_cov),1)}
             # else{selfis<-selfis_cov*selfis_eff}
             selfis<-interventions$isolation$cov[t_]*interventions$isolation$eff[t_]
           }else{selfis<-0}
           #####
           dist_cov <- interventions$distancing$cov[t_]
           dist_eff <- interventions$distancing$eff[t_]
           #####
           hand_eff <- interventions$handwashing$eff[t_]
           #####
           work_cov <- interventions$workhome$cov[t_]
           work_eff <- interventions$workhome$eff[t_]
           ####################################################
           school_eff_all <- interventions$schoolclosing$eff[t_]
           school_eff_1 <- interventions$schoolclosing$eff1[t_]
           school_eff_2 <- interventions$schoolclosing$eff2[t_]
           school_eff_3 <- interventions$schoolclosing$eff3[t_]
           school_eff_4 <- interventions$schoolclosing$eff4[t_]
           school_cov_all <- interventions$schoolclosing$cov[t_]
           school_cov_1 <- interventions$schoolclosing$cov1[t_]
           school_cov_2 <- interventions$schoolclosing$cov2[t_]
           school_cov_3 <- interventions$schoolclosing$cov3[t_]
           school_cov_4 <- interventions$schoolclosing$cov4[t_]
           ####################################################
           cocoon_cov <- interventions$cocoon$cov[t_]
           cocoon_eff <- interventions$cocoon$eff[t_]
           age_cocoon <- interventions$cocoon$cocoon_age[t_]
           ####################################################
           quarantine_cov <- interpola(interventions$quarantine$cov, t)
           quarantine_days <- interpola(interventions$quarantine$days, t)
           tau_window <- interpola(interventions$quarantine$window_time, t)
           # ####################################################
           if(interventions$testing$vec[t_]){
             num_tests <- interpola(interventions$testing$tests_day, t)
             asymp_eff <- interventions$testing$asymp_eff[t_]
             flux_X <- sum(gamma*pclin[,2]*selfis*(1-scale_ihr*ihr[,2])*(E+QE))
             flux_CL <- sum(gamma*pclin*(1-selfis)*(1-scale_ihr*ihr[,2])*(E+QE))
             flux_H <- sum(gamma*scale_ihr*ihr[,2]*(1-prob_icu)*(1-critH)*(E+QE))
             flux_HC <- sum(gamma*scale_ihr*ihr[,2]*(1-prob_icu)*critH*(E+QE))
             flux_ICU <- sum(gamma*scale_ihr*ihr[,2]*prob_icu*(1-critICU)*(E+QE))
             flux_ICUH <- sum(gamma*scale_ihr*ihr[,2]*prob_icu*critICU*(1-critICUH)*(E+QE))
             flux_ICUC <- sum(gamma*scale_ihr*ihr[,2]*prob_icu*critICU*critICUH*(E+QE))
             
             cases <- c(flux_ICU, flux_ICUH, flux_H, flux_ICUC, flux_HC, flux_X+flux_CL)
             pis <- pmin(pmax((num_tests - cumsum(c(0, cases[-length(cases)]))) / (cases+1), 0), 1)
             c(pt_icu, pt_icuh, pt_h, pt_icuc, pt_hc, pt_cl) %<-% pis
             pt_x <- pt_cl
             
           } else{
             pt_s <- interventions$testing$test_s[t_]
             pt_e <- interventions$testing$test_e[t_]
             pt_i <- interventions$testing$test_i[t_]
             pt_cl <- interventions$testing$test_cl[t_]
             pt_x <- interventions$testing$test_cl[t_]
             pt_icu <- interventions$testing$test_icu[t_]
             pt_icuh <- interventions$testing$test_icuh[t_]
             pt_r <- interventions$testing$test_r[t_]
             pt_h <- interventions$testing$test_h[t_]
             pt_icuc <- interventions$testing$test_icuc[t_]
             pt_hc <- interventions$testing$test_hc[t_]
           }
           # ####################################################
           home_tracing_eff <- interpola(interventions$tracing$home_tracing_eff, t)
           school_tracing_eff <- interpola(interventions$tracing$school_tracing_eff, t)
           work_tracing_eff <- interpola(interventions$tracing$work_tracing_eff, t)
           other_tracing_eff <- interpola(interventions$tracing$other_tracing_eff, t)
           ####################################################
           if (workhome){
             work<-interventions$workhome$cov[t_]*interventions$workhome$eff[t_]
           }else{work<-1}

           # print(schoolclose_all)
           if (schoolclose_all){#####if uniform closing activated#####
             
             school_all_cov<-interventions$schoolclosing$cov[t_]
             school_all_eff<-interventions$schoolclosing$eff[t_]
             school_cov_1<-school_all_cov
             school_eff_1<-school_all_eff
             school_cov_2<-school_all_cov
             school_eff_2<-school_all_eff
             school_cov_3<-school_all_cov
             school_eff_3<-school_all_eff
             school_cov_4<-school_all_cov
             school_eff_4<-school_all_eff
             
           }
           else{
             if (schoolclose_1){##########SCHOOL CLOSE 1########
               school_cov_1<-interventions$schoolclosing$cov1[t_]
               school_eff_1<-interventions$schoolclosing$eff1[t_]
             }
             else{
               school_cov_1 <- 0
               school_eff_1 <- 0
             }
             if (schoolclose_2){#########SCHOOL CLOSE 2########
               school_cov_2<-interventions$schoolclosing$cov2[t_]
               school_eff_2<-interventions$schoolclosing$eff2[t_]
             }
             else{
               school_cov_2 <- 0
               school_eff_2 <- 0
             }
             if (schoolclose_3){#######SCHOOL CLOSE 3###########
               school_cov_3<-interventions$schoolclosing$cov3[t_]
               school_eff_3<-interventions$schoolclosing$eff3[t_]
             }
             else{
               school_cov_3 <- 0
               school_eff_3 <- 0
             }
             if (schoolclose_4){######SCHOOL CLOSE 4#############
               school_cov_4<-interventions$schoolclosing$cov4[t_]
               school_eff_4<-interventions$schoolclosing$eff4[t_]
             }
             else{
               school_cov_4 <- 0
               school_eff_4 <- 0
             }
           }
           if(distancing){
             dist<- interventions$distancing$cov[t_]*interventions$distancing$eff[t_]
           }else{dist<-1}
           if(handwash){
             hand<-interventions$handwashing$eff[t_]
           }else{hand<-0}
           if(vaccine){
             vac_rate <- (-log(1-vaccine_cov)/vac_campaign)
             vaccinate <- vac_rate
             vaccinate <- 0
           }else{vaccinate<-0}
           if(travelban){
             trvban_eff<- interventions$travelban$eff[t_]
           }else{trvban_eff<-0}
           ##############cocooning the elderly#####################
           cocoon_mat<-diag(1,nrow = length(popstruc$pop),ncol = length(popstruc$pop))
           for(jj in (age_cocoon-1):length(popstruc$pop)){cocoon_mat[jj,jj]<-(1-cocoon_cov*cocoon_eff)}
           
           
           ########################################################
           kappa<- (school_cov_1*school_eff_1*popstruc[1,2] + school_cov_2*school_eff_2*popstruc[2,2]+
                      school_cov_3*school_eff_3*popstruc[3,2] + school_cov_4*school_eff_4*popstruc[4,2])/sum(popstruc[1:4,2])
           beta<-c()
           beta[1]<-1-school_cov_1*school_eff_1
           beta[2]<-1-school_cov_2*school_eff_2
           beta[3]<-1-school_cov_3*school_eff_3
           beta[4]<-1-school_cov_4*school_eff_4
           school_mat<-diag(c(beta,rep(1-kappa,A-4)))      
           
           ##################PROPORTION OF CONTACTS OF EACH CATEGORY, COMPUTING THE AVERAGE CONTACTS#################
           P_school <- popstruc[,2]%*%(diag(1,A)-school_mat)%*%contact_school%*%(diag(1,A)-school_mat)%*%popstruc[,2]
           P_school_2<- popstruc[,2]%*%contact_school%*%popstruc[,2]
           P_work <- popstruc[,2] %*% contact_work %*% popstruc[,2]
           P_other <- popstruc[,2] %*% contact_other %*% popstruc[,2]
           
           p_work <- (P_work / (P_work + P_school_2 + P_other))[1,1]
           p_school <- (P_school/ (P_work + P_school_2 + P_other))[1,1]
           p_other <- (P_other / (P_work + P_school_2 + P_other))[1,1]
           
           weighted_intervention <- (workhome) * p_work * work_cov * work_eff +
             p_school+
             (distancing) * p_other * dist_cov * dist_eff
           # home_steep: take from parameters
           ###########percolation value computed from least squares, refer to overleaf to formula##########
           home_effective <- 1 - home_eff * (1+tanh(home_steep*(weighted_intervention - perc_threshold )))/2
           # contact matrices
           contact_home <- home_effective*contact_home
           cts<-(contact_home+distancing*(1-dist)*contact_other+(1-distancing)*contact_other
                 +school_mat%*%contact_school%*%school_mat # school close
                 +kappa*contact_home*s2h # inflating contacts at home when school closes
                 +(1-workhome)*contact_work  # normal work
                 +workhome*(1-work)*contact_work # people not working from home when homework is active
                 +contact_home*workhome*work*w2h # inflating contacts at home when working from home
           )
           
           tin<-(parameters["startdate"]-as.numeric(as.Date("2020-01-01")))/365.25
           
           # Final transmission related parameters
           contacts <- cocoon_mat%*%cts%*%cocoon_mat
           seas <- 1+amp*cos(2*pi*(t+tin-(phi*365.25/12))/365.25)
           importation <- mean_imports*(1-trvban_eff)
           HH<-H+ICUH + ICU
           HHC<-HC + ICUC 
           QQ<-QS + QE + QI+ QR + QC 
           QInf <- rho*QE + QI + QC + X + HHC
           ##these are the terms of not quarantined person
           lam <- (1-hand)*p*seas*(contacts%*%((rho*E+(I+CL+importation)+rhos*(HH))/P)) +
           ######^ contacts of quarantined people that infects non quarantined people
            (1-hand)*p*seas*(1-quarantine_eff_other)*(cocoon_mat%*%contact_other%*%cocoon_mat%*%QInf/P) +
           #self isolated can infect people in own household
            (1-hand)*p*seas*(1-quarantine_eff_home)*cocoon_mat%*%contact_home%*%cocoon_mat%*%QInf/P 
           
           #quarantined infecting quarantined
           lamq <- (1-hand)*p*seas*(1-quarantine_eff_home)*cocoon_mat%*%contact_home%*%cocoon_mat%*%QInf/P +
           ######^not quarantined infecting quarantined
            (1-hand)*p*seas*(1-quarantine_eff_other)*(cocoon_mat%*%contact_other%*%cocoon_mat%*%((rho*E+(I+CL+importation)+rhos*(HH))/P))
           # birth/death
           b1<-sum(popbirth*popstruc[,2])
           birth<-0*popbirth
           birth[1]<-b1
           c_home <- home_effective*contact_home #contacts in home
           c_home <- c_home*(1+kappa*s2h + workhome*(1-work)*w2h) #adding contacts from school and work to home
           c_school <- school_mat%*%contact_school%*%school_mat # school close
           c_work <- (1-workhome)*contact_work + workhome*(1-work)*contact_work
           c_others <- distancing*(1-dist)*contact_other+(1-distancing)*contact_other
           ######reconstructed separated matrix, begin to reconstruct the cocooned matrix############
           c_home<-cocoon_mat%*%c_home%*%cocoon_mat
           c_school <- cocoon_mat%*%c_school%*%cocoon_mat
           c_work <- cocoon_mat%*%c_work%*%cocoon_mat
           c_others <-cocoon_mat%*%c_others%*%cocoon_mat
           ##########################################################################################
           flux_X <- gamma*pclin[,2]*selfis*(1-scale_ihr*ihr[,2])*E
           flux_CL <- gamma*pclin[,2]*(1-selfis)*(1-scale_ihr*ihr[,2])*E
           flux_H <- gamma*scale_ihr*ihr[,2]*(1-prob_icu)*(1-critH)*E
           flux_HC <- gamma*scale_ihr*ihr[,2]*(1-prob_icu)*critH*E
           flux_ICU <- gamma*scale_ihr*ihr[,2]*prob_icu*(1-critICU)*E
           flux_ICUH <- gamma*scale_ihr*ihr[,2]*prob_icu*critICU*(1-critICUH)*E
           flux_ICUC <- gamma*scale_ihr*ihr[,2]*prob_icu*critICU*critICUH*E
           flux_I <- gamma*(1-pclin[,2])*(1-scale_ihr*ihr[,2])*E
           ##########################################################################################
           final_cts <- home_tracing_eff*c_home + school_tracing_eff*c_school + work_tracing_eff*c_work + others_tracing_eff*c_others
           Q_in <- quarantine_cov*tau_window/(P-QQ)*final_cts%*%(pt_x*flux_X + pt_cl*flux_CL + pt_h*flux_H + pt_hc*flux_HC + pt_icu*flux_ICU + pt_icuh*flux_ICUH + pt_icuc*flux_ICUC)

           # calculando pt_e e pt_i
           if(interventions$testing$vec[t_]){
             # TODO: definir test_pos_overdp na lista de variáveis
             test_pos_overdp <- 1 # (> 1, provavelmente 2)
             pt_e <- pt_i <- pmin(pmax(test_pos_overdp * (num_tests - sum(cases)) / (Q_in+1) * (I + E + S + R), 0), 1)
           }
           Q_in_2 <- quarantine_cov*tau_window/(P-QQ)*final_cts%*%Q_in * (pt_i*I + pt_e*E)
           #entrance rates in quarantine to each compartment
           ##################################################################################################################
           #####"standard" equations
           dSdt <- -S*lam - S*vaccinate + omega*R + ageing%*%S - mort*S + birth -
             (Q_in+Q_in_2)*S + (1/quarantine_days)*QS
           dEdt <- S*lam - gamma*E + ageing%*%E - mort*E + (1-vaccine_eff)*lam*V -
             (Q_in+Q_in_2)*E + (1/quarantine_days)*QE
           dIdt <- gamma*(1-pclin[,2])*(1-scale_ihr*ihr[,2])*E - nui*I + ageing%*%I - mort*I +
             (1/quarantine_days)*QI - (Q_in+Q_in_2)*I
           dRdt <- nui*I - omega*R + nui*X + nui*CL + ageing%*%R - mort*R +
             (1/quarantine_days)*QR - (Q_in+Q_in_2)*R + 
             nus*(1-pdeath_h*ifr[,4])*H + (1-pdeath_icu*ifr[,3])*nu_icu*ICU +
             (1-pdeath_icuc*ifr[,3])*nu_icuc*ICUC + (1-pdeath_hc*ifr[,4])*nusc*HC + (1-pdeath_icuh*ifr[,3])*nu_icuh*ICUH
           dCLdt <- (1 - quarantine_cov*pt_cl)*gamma*pclin[,2]*(1-selfis)*(1-scale_ihr*ihr[,2])*E - nui*CL + ageing%*%CL - mort*CL + (1/quarantine_days)*QC
           dXdt <- gamma*selfis*pclin[,2]*(1-scale_ihr*ihr[,2])*E - nui*X + ageing%*%X - mort*X
           #######hospitalization related equations
           #### H, ICU, ICUH ARE NOT QUARANTINED, BUT HC AND HICUC ARE#############
           dHdt <- gamma*scale_ihr*ihr[,2]*(1-prob_icu)*(1-critH)*(E+QE) - nus*H + ageing%*%H-mort*H
           dICUdt <- gamma*scale_ihr*ihr[,2]*prob_icu*(1-critICU)*(E+QE) - nu_icu*ICU + ageing%*%ICU - mort*ICU
           dICUHdt <- gamma*scale_ihr*ihr[,2]*prob_icu*critICU*(1-critICUH)*(E+QE) - nu_icuh*ICUH + ageing%*%ICUH - mort*ICUH
           dHCdt <- gamma*scale_ihr*ihr[,2]*(1-prob_icu)*critH*(E+QE) - nusc*HC + ageing%*%HC-mort*HC
           dICUCdt <- gamma*scale_ihr*ihr[,2]*prob_icu*critICU*critICUH*(E+QE) - nu_icuc*ICUC + ageing%*%ICUC - mort*ICUC
           ########quarantine related equations - WORK IN PROGRESS##################################################
           dQSdt <- (Q_in + Q_in_2)*S + ageing%*%QS - mort*QS - (1/quarantine_days)*QS - lamq*QS + omega*QR
           dQEdt <- (Q_in + Q_in_2)*E - gamma*QE +ageing%*%QE-mort*QE - (1/quarantine_days)*QE +lamq*QS
           dQIdt <- (Q_in + Q_in_2)*I + gamma*(1-scale_ihr*ihr[,2])*(1-pclin[,2])*QE-nui*QI+ageing%*%QI-mort*QI - (1/quarantine_days)*QI
           dQRdt <- (Q_in + Q_in_2)*R + nui*QI - omega*QR + nui*QC + ageing%*%QR - mort*QR - (1/quarantine_days)*QR
           dQCdt <- (1-selfis)*(quarantine_cov*pt_cl)*gamma*pclin[,2]*(1-scale_ihr*ihr[,2])*E+gamma*pclin[,2]*(1-scale_ihr*ihr[,2])*QE - nui*QC+ageing%*%QC-mort*QC - (1/quarantine_days)*QC
           #############VACCINATION - NOT IMPLEMENTED CORRECTLY YET#############
           dVdt <- vaccinate*S -(1-vaccine_eff)*lam*V +ageing%*%V - mort*V
           ######case and death reporting equations####
           dCdt <- report*gamma*(1-scale_ihr*ihr[,2])*(1-pclin[,2])*(E+QE)+reportc*gamma*(1-scale_ihr*ihr[,2])*pclin[,2]*(E+QE)+ reporth*gamma*scale_ihr*ihr[,2]*(E+QE)
           dCMdt <- nus*pdeath_h*ifr[,4]*H + nusc*pdeath_hc*ifr[,4]*HC +
             nu_icu*pdeath_icu*ifr[,3]*ICU + nu_icuc*pdeath_icuc*ifr[,3]*ICUC + pdeath_icuh*ifr[,3]*nu_icuh*ICUH + 
             mort*(H + HC + ICU + ICUC + ICUH + CL + X)
           dCMCdt <- nusc*pdeath_hc*ifr[,4]*HC +  nu_icuc*pdeath_icuc*ifr[,3]*ICUC + pdeath_icuh*ifr[,3]*nu_icuh*ICUH + mort*HC + mort*ICUC
           
           # return the rate of change
           list(c(dSdt,dEdt,dIdt,dRdt,dXdt,dHdt,dHCdt,dCdt,dCMdt,dVdt,dQSdt,dQEdt,dQIdt,dQRdt,dCLdt,dQCdt,dICUdt,dICUHdt,dICUCdt,dCMCdt))
         }
    ) 
  }
  
  out <- ode(y = Y, times = times, func = covid, parms = params)
  
  return(list(parameters, out))
}

