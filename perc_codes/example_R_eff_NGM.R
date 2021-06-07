if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(zeallot)){install.packages("zeallot"); library(zeallot)} # provides %<-%
if(!exists("sheets_read",mode="function")) library(googlesheets4)

source("import_functions_v2.R")
source("R_eff_NGM.R")
source("etaclass.R")


### Define state (estado), city (cidade) and model_version
estado =  "SP"; cidade = "Sao_Paulo"; model_version = "model1"

### Load demographic data
popstruc <- data.frame(read.csv(paste0("DATA/", estado, "/", cidade, "/13. DistrEtaria_cidade2020.csv")))
A <-  length(popstruc$pop)
popmort <- pop.mort.br # natural mortality per person per year
mort <- popmort[,2]/(popstruc[,2]*365.25) # convert from 1000s per year period to per person per day

## per year ageing matrix
dd <- seq(1:A)/seq(1:A)
ageing <- t(diff(diag(dd),lag = 1)/(5*365.25))
ageing <- cbind(ageing,0*seq(1:A)) 
demographic_parameters = list('mort' = mort,'ageing' = ageing)


### Carregar modelo e parâmetros de acordo com a cidade

base_parameters <- c()
base_parameters <- update_parameters(base_parameters)
base_parameters <- update_params_from_fit(base_parameters,paste0('fits/', get_latest_fit(estado = estado, model_version = model_version)[1]))
base_parameters <- set_scale(base_parameters)


##### Manually defining parameters, just for debugging
##base_parameters <- list('p' = 0.04, 'rho' = 0.5, 'omega' = 0, 'gamma' = 1/7, 'nui' = 1/11, 'nus' = 1/11, 'rhos' = 0.2, 'nusc' = 1/21, 'nu_icu' = 1/21, 'nu_icuh' = 1/21, 'nu_icuc' = 1/21)  

### Carregar matriz de contato e definir
load("PREM/prem2020_bra.Rdata")
c_home <- u_home
c_school <- u_school
c_work <- u_work
c_other <- u_others



nce <- A-length(c_home[1,])

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

#interventions <- update_scenario_parameters(base_parameters, scenario=4, estado=estado, cidade=cidade)
#####FROM HERE THE PARAMETERS ARE DYNAMIC IN MODEL SIMULATIONS#####
### Health system performance
H <- 0*popstruc[,2]
ICU <- 0*popstruc[,2]
ICUH <- 0*popstruc[,2]
critH <- etaclass(sum(H, ICUH), base_parameters['beds_available'], base_parameters['give'])
critICU <- etaclass(sum(ICU), base_parameters['icu_beds_available'], base_parameters['give'])
critICUH <- critH^2
hospitalization_parameters = list("critH" = critH, "critICU" = critICU, "critICUH" = critICUH, 'scale_ihr' = 1)

wnpi <- c()
R_eff_perc <- c()
R_eff_noperc <- c()
dscale <- 0.01
for(scale in seq(0,1,by=dscale)){
    index=round(scale/dscale)
    ## NPI scenario example ###
    ##handwashing
    hand <- 0.15
    ##work
    work <- scale * 1
    ##school (unstructured form - equivalent to setting all ages to the same covs and effs to all scholar ages in current model formulation)
    school  <- scale * 1
    ##distancing
    dist <- scale * 1
    ##self isolation
    selfis <- 0.3

    ## Building Correspondent Weighted Interventions
    P_school <- popstruc[,2] %*% contact_school %*% popstruc[,2]
    P_work <- popstruc[,2] %*% contact_work %*% popstruc[,2]
    P_other <- popstruc[,2] %*% contact_other %*% popstruc[,2]
    ##calculating p's       
    p_work <- (P_work / (P_work + P_school + P_other))[1,1]
    p_school <- (P_school/ (P_work + P_school + P_other))[1,1]
    p_other <- (P_other / (P_work + P_school + P_other))[1,1]       
    weighted_intervention <- p_work * work + p_school * school + p_other * dist
    wnpi[index] <- weighted_intervention
    #this is a list containg the names of the interventions, in dynamical tracking of R_NGM, building this relies on update_scenario_parameters and parsing the values in each point in time is required
    simplified_interventions<-list("hand" = hand, "work" = work, "school" = school, "dist" = dist, "selfis" = selfis)

    ## Percolation Function
    home_eff <- 1
    home_effective <- with(as.list(base_parameters),{ home_eff *(1+tanh(home_steep*(weighted_intervention - perc_threshold )))/2})
    cocoon_mat <- diag(1,nrow = A,ncol = A)
    for(jj in (base_parameters['age_cocoon'] - 1):length(popstruc$pop)){cocoon_mat[jj,jj] <- (1 - base_parameters['cocoon_cov'] * base_parameters['cocoon_eff'])}
    ## NPIs modified Contact Matrix with percolation
    C <- (1 - home_effective) * contact_home + (1 - school) * contact_school + (1 - work) * contact_work + (1 - dist) * contact_other
    C <- cocoon_mat %*% C %*% cocoon_mat
    C_home <-  (1 - home_effective) * (cocoon_mat %*% contact_home %*% cocoon_mat)
    C_other <- cocoon_mat %*% contact_other %*% cocoon_mat

    ## NPIs modified Contact Matrix without percolation
    C_np <- contact_home + (1 - school) * contact_school + (1 - work) * contact_work + (1 - dist) * contact_other
    C_np <- cocoon_mat %*% C_np %*% cocoon_mat
    C_home_np <- (cocoon_mat %*% contact_home %*% cocoon_mat)

    ##Assuming all susceptible
    S <- popstruc[,2]
    P <- popstruc[,2]
    ## Call for NGM
    R_eff_perc[index] <- R_NGM(S, C, C_home, C_other, P, base_parameters, demographic_parameters, simplified_interventions, hospitalization_parameters)
    R_eff_noperc[index] <- R_NGM(S, C_np, C_home_np, C_other, P, base_parameters, demographic_parameters, simplified_interventions, hospitalization_parameters)
}


df <- data.frame(weighted_npi = wnpi,
                 R_perc = R_eff_perc,
                 R_noperc = R_eff_noperc)

write.csv(df,'reffs_to_plot.csv')
