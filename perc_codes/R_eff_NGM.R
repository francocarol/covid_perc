######R_NGM#######
##Calculates Effective R via Next-Generation Matrix Method
##Receives a 1d-array, S, containing Susceptible age distribution
##Receives the full contact matrix,C, in which NPIs are already included.
##Receive the C_home (with percolation, hand and cocooning) and C_other (with only hand and cocooning?) contact matrices separatedly
##Receives the Total population, P
##Receives parameters as defined in example_R_eff_NGM.R
##Receives demographic parameters as defined in example_R_eff_NGM.R
##Receives simplified interventions as defined in example_R_eff_NGM.R
##Receives hospitalization_parameters as defined in example_R_eff_NGM.R
##For now it doesn't account for quarantining (nor vaccination)
##Still lacks a little bit of work to fit the current work-flow (mainly to remove simplified_interventions and use the interventions object created by  import_functions_V2.R update_scenario_parameters function) 
##Will try to run on a single simulation in order to track R_eff
if(!require(magic)){install.packages("magic"); library(magic)}
if(!require(rARPACK)){install.packages("rARPACK"); library(rARPACK)}
if(!require(matlib)){install.packages("matlib"); library(matlib)}






R_NGM <- function(S, C, C_home, C_other, P, parameters, demographic_parameters, simplified_interventions, hospitalization_parameters){
    with(as.list(c(parameters,demographic_parameters, simplified_interventions, hospitalization_parameters)),{
        n <- length(c('E', 'I', 'CL', 'X', 'H', 'HC', 'ICU', 'ICUH', 'ICUC')) # order of infected classes, important in weather the classes are isolated (use sml_q) or not (use sml)
        A <- length(S) 
        P <- 1 / P
        sml <-  diag(S) %*% C %*% diag(P) 
        sml_q <- diag(S) %*% (C_home + C_other) %*% diag(P)
        ##infected entries come only at E classes (the first A lines of F matrix are the only non-null ones)
        ##         E        , I  , CL , X    , H         , HC   , ICU       , ICUH      , ICUC
        ##         1        , 2  , 3  , 4    , 5         , 6    , 7         , 8         , 9 
        f <- cbind(rho * sml, sml, sml, sml_q, rhos * sml, sml_q, rhos * sml, rhos * sml, sml_q)         
        f_full <-  (1 - hand) * p * rbind(f, matrix(0, nrow = (n-1) * A, ncol = n * A))
        ##Define the exit matrix, V. It's composed of a block diagonal matrix, v_diag, and has it's first A collums and n*A rows as the movement from exposed to other infected classes, in v_1row
        v_diag <- adiag(diag((gamma + mort)) - ageing, 
                        diag((nui + mort)) - ageing, 
                        diag((nui + mort)) - ageing, 
                        diag((nui + mort)) - ageing, 
                        diag((nus + mort)) - ageing, 
                        diag((nusc + mort)) - ageing, 
                        diag((nu_icu + mort)) - ageing, 
                        diag((nu_icuh + mort)) - ageing, 
                        diag((nu_icuc + mort)) - ageing)
        v_1row <- rbind(matrix(0, nrow = A,ncol = A),
                        -diag(gamma * (1 - pclin[,2]) * (1 - scale_ihr * ihr[,2])), 
                        -diag((1 - quarantine_cov * pt_cl) * gamma * pclin[,2] * (1 - selfis) * (1 - scale_ihr * ihr[,2])),
                        -diag(gamma * selfis * pclin[,2] * (1 - scale_ihr * ihr[,2])),
                        -diag(gamma * scale_ihr * ihr[,2] * (1 - prob_icu) * (1 - critH)),
                        -diag(gamma * scale_ihr * ihr[,2] * (1 - prob_icu) * critH),
                        -diag(gamma * scale_ihr * ihr[,2] * prob_icu * (1 - critICU)),
                        -diag(gamma * scale_ihr * ihr[,2] * prob_icu * critICU * (1 - critICUH)),
                        -diag(gamma * scale_ihr * ihr[,2] * prob_icu * critICU * critICUH))
        v_1row <- cbind(v_1row, matrix(0, nrow = n*A, ncol = (n-1)*A)) # set all to n*A by n*A dimension
        v <- v_diag + v_1row
        m <- f_full %*% solve(v) #Proper NGM
        values <- eigs(m, 1, which = 'LM', sigma = NULL)
        return(as.numeric(values$values[1]))
    })
}


set_scale <- function(parameters){
    parameters["rho"]<-parameters["rho"]/100
    parameters["omega"]<-(1/(parameters["omega"]*365))
    parameters["gamma"]<-1/parameters["gamma"]
    parameters["nui"]<-1/parameters["nui"]
    parameters["report"]<-parameters["report"]/100
    parameters["reportc"]<-parameters["reportc"]/100
    parameters["reporth"]<-parameters["reporth"]/100
    parameters["nus"]<-1/parameters["nus"]
    parameters["rhos"]<-parameters["rhos"]/100
    ##parameters["amp"]<-parameters["amp"]/100
    ##parameters["selfis_dur"]<-parameters["selfis_dur"]*7
    ##parameters["selfis_cov"]<-parameters["selfis_cov"]/100
    ##parameters["selfis_eff"]<-parameters["selfis_eff"]/100
    ##parameters["dist_dur"]<-parameters["dist_dur"]*7
    ##parameters["dist_cov"]<-parameters["dist_cov"]/100
    ##parameters["dist_eff"]<-parameters["dist_eff"]/100
    ##parameters["hand_dur"]<-parameters["hand_dur"]*7
    ##parameters["hand_eff"]<-parameters["hand_eff"]/100
    ##parameters["work_dur"]<-parameters["work_dur"]*7
    ##parameters["work_cov"]<-parameters["work_cov"]/100
    ##parameters["work_eff"]<-parameters["work_eff"]/100
    ##parameters["w2h"]<-parameters["w2h"]/100
    ##parameters["school_dur"]<-parameters["school_dur"]*7
    ##parameters["school_eff"]<-parameters["school_eff"]/100
    ##parameters["s2h"]<-parameters["s2h"]/100
    parameters["cocoon_dur"]<-parameters["cocoon_dur"]*7
    parameters["cocoon_cov"]<-parameters["cocoon_cov"]/100
    parameters["cocoon_eff"]<-parameters["cocoon_eff"]/100
    parameters["age_cocoon"]<-floor((parameters["age_cocoon"]/5)+1)
    ##parameters["travelban_eff"]<-parameters["travelban_eff"]/100
    ##parameters["vaccine_eff"]<-parameters["vaccine_eff"]/100
    ##parameters["vaccine_cov"]<-parameters["vaccine_cov"]/100
    ##parameters["vac_campaign"]<-parameters["vac_campaign"]*7
    ##parameters["travelban_dur"]<-parameters["travelban_dur"]*7
    ##parameters["screen_dur"]<-parameters["screen_dur"]*7
    ##parameters["screen_cov"]<-parameters["screen_cov"]/100
    parameters["quarantine_cov"]<-parameters["quarantine_cov"]/100
    ##parameters["quarantine_days"]<-parameters["quarantine_days"]
    ##parameters["quarantine_eff_home"]<-parameters["quarantine_eff_home"]/100
    ##parameters["quarantine_eff_other"]<-parameters["quarantine_eff_other"]/100
    ##parameters["ratem"]<-1/parameters["ratem"]
    ##parameters["ratemHC"]<-1/parameters["ratemHC"]
    ##parameters["ratemICU"]<-1/parameters["ratemICU"]
    ##parameters["ratemICUC"]<-1/parameters["ratemICUC"]
    ##parameters["ratemVent"]<-1/parameters["ratemVent"]
    ##parameters["ratemVentC"]<-1/parameters["ratemVentC"]
    parameters["give"]<-parameters["give"]/100
    ##parameters["pdeath_h"]<-parameters["pdeath_h"]/100
    ##parameters["pdeath_hc"]<-parameters["pdeath_hc"]/100
    ##parameters["pdeath_icu"]<-parameters["pdeath_icu"]/100
    ##parameters["pdeath_icuh"]<-parameters["pdeath_icuh"]/100
    ##parameters["pdeath_icuc"]<-parameters["pdeath_icuc"]/100
    ##parameters["pdeath_vent"]<-parameters["pdeath_vent"]/100
    ##parameters["pdeath_venticu"]<-parameters["pdeath_venticu"]/100
    ##parameters["pdeath_venth"]<-parameters["pdeath_venth"]/100
    ##parameters["pdeath_ventc"]<-parameters["pdeath_ventc"]/100
    parameters["nusc"]<-1/parameters["nusc"]
    parameters["nu_icu"]<-1/parameters["nu_icu"]
    parameters["nu_icuh"]<-1/parameters["nu_icuh"]
    parameters["nu_icuc"]<-1/parameters["nu_icuc"]
    ##parameters["nu_vent"]<-1/parameters["nu_vent"]
    ##parameters["nu_venticu"]<-1/parameters["nu_venticu"]
    ##parameters["nu_venth"]<-1/parameters["nu_venth"]
    ##parameters["nu_ventc"]<-1/parameters["nu_ventc"]
    ##parameters["pclin"]<-parameters["pclin"]/100
    ##parameters["prob_vent"]<-parameters["prob_vent"]/100
    ##parameters["beds_available"] = sum(popstruc[,2])*parameters["beds"]/1000 # maximum number of hospital beds - numeric 
    ##parameters["icu_beds_available"] = sum(popstruc[,2])*parameters["ICU_beds"]/ #8000, # maximum number of hospital beds - numeric
    return(parameters)
}
