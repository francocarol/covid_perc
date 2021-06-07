################################################
################################################
################################################

### EVAL AIC from fitting results

################################################
################################################
################################################

calc_AIC <- function(file1, file2, file3) {
    files <- c(file1, file2, file3)
    data1 <- read.csv(file1)
    data2 <- read.csv(file2)
    data3 <- read.csv(file3)
    # n_par <- number of parameters in the fitted model
    # ?? should I actually consider all of them (sp. stardate) here?
    # k=2 for the usual AIC, or log(n) (n being the number of observations)
    k <- 2
    
    #AIC 1: full model (with percolation)
    log_likelyhood <- data1$log_like
    startdate1      <- as.character(data1$data[which.min(log_likelyhood)])
    p1              <- data1$p[which.min(log_likelyhood)]
    perc1           <- data1$perc_threshold[which.min(log_likelyhood)]
    h_steep1        <- data1$home_steep[which.min(log_likelyhood)]
    
    n_par_1 <- 4 # considering the full model ( p, T_perc, home_steep and startdate)
    AIC_1  <- +2*log_likelyhood[which.min(log_likelyhood)] + k*n_par_1
    
    #AIC 2: T_perc = home_steep = 0
    log_likelyhood <- data2$log_like
    startdate2      <- as.character(data2$data[which.min(log_likelyhood)])
    p2              <- data2$p[which.min(log_likelyhood)]
    perc2           <- NA
    h_steep2        <- NA
    
    n_par_2 <- 2 # considering the full model ( p and startdate)
    AIC_2  <- +2*log_likelyhood[which.min(log_likelyhood)] + k*n_par_2
    
    #AIC 3: T_perc = home_steep = 0 
    log_likelyhood <- data3$log_like
    startdate3      <- as.character(data3$data[which.min(log_likelyhood)])
    p3              <- data3$p[which.min(log_likelyhood)]
    perc3           <- NA
    h_steep3        <- NA
    
    n_par_3 <- 2 # considering the full model ( p and startdate)
    AIC_3  <- +2*log_likelyhood[which.min(log_likelyhood)] + k*n_par_3
    # + because it is the neg log_likelyhood
    
    # AICc is useful only when # observations / # parameters < 40
    # no need to use AICc because the sample size (number of data points?) is not so small (~300)
    AIC     <- c(AIC_1, AIC_2, AIC_3) 
    dAIC    <- AIC - min(AIC)
    #
    #weights <- exp(-1/(2*dAIC))/sum(exp(-1/(2*dAIC)))
    # this probably is not right
    p         <- c(p1, p2, p3)
    startdate <- c(startdate1, startdate2, startdate3)
    perc      <- c(perc1, perc2, perc3)
    h_steep   <- c(h_steep1, h_steep2, h_steep3)
    
    # AICc is useful only when # observations / # parameters < 40
    # no need to use AICc because the sample size (number of data points?) is not so small (~300)
    
    table <- cbind(Model=c("Version 1","Version 2", "Version 3"), p, startdate, h_steep, perc, AIC, dAIC, files)
}
