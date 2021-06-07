if(!require(xts)){install.packages("xts"); library(xts)}

source("plots.R")
# get observed data to compare

###################################

test_modifying_parameters <- function(parameters, interventions, Y, params = NA,
                                      flag = "covid", base_size=18, rollwindow = 1,
                                      rowmean_line_color = "grey",rowmean_line_size = 1.2,
                                      weekly = TRUE, max.date=NULL){
    plot.formatos <- theme_bw(base_size=base_size)
    for (i in names(params)){parameters[i] <- params[i]}
    
    out <- run_covid19(parameters, interventions, Y)
    res <- result_ts(out[[2]], parameters, parameters["startdate"], diff=F)
    
    
    ##res <- result_ts(out[[2]], parameters, date_parameters, diff=F)
    
    if(flag == "covid"){
        gC <- ggplot(res) +
            geom_line(aes(x=Index, y=C)) +
            geom_point(data=now.covid.zoo, aes(x=Index, y=now.covid.zoo)) +
            xlim(as.Date('2020-3-1'), max(time(now.covid.zoo))) +
            ylim(c(0, max(now.covid.zoo, window(res$C, end=max(time(now.covid.zoo)))))) +
            plot.formatos
        gM <- ggplot(res) +
            geom_line(aes(x=Index, y=CM)) +
            geom_point(data=now.obito.covid.zoo, aes(x=time(now.obito.covid.zoo),
                                                     y=now.obito.covid.zoo)) +
            xlim(as.Date('2020-3-1'), max(time(now.obito.covid.zoo))) +
            ylim(c(0, max(now.obito.covid.zoo, window(res$CM, end=max(time(now.obito.covid.zoo)))))) +
            plot.formatos
        
        data_c = diff(now.covid.zoo)
        data_m = diff(now.obito.covid.zoo)
        data_c_fit = diff(res$C)
        data_m_fit = diff(res$CM)
        
        if(!is.null(max.date)){
            xmax_c=as.Date(max.date)
            xmax_m=as.Date(max.date)
        }
        else{
            xmax_c=max(time(data_c))
            xmax_m=max(time(data_m))
        }
        
        if(weekly) {
            model_c_xts <- xts::apply.weekly(data_c_fit, sum)
            model_cm_xts <- xts::apply.weekly(data_m_fit, sum)
            data_c_fit <- zoo(model_c_xts, as.Date(index(model_c_xts)))
            data_m_fit <- zoo(model_cm_xts, as.Date(index(model_cm_xts)))
        }
        
        gCdiff <- ggplot(data_c_fit) +
            geom_line(aes(x=time(data_c_fit), y=data_c_fit)) +
            geom_point(data=data_c, aes(x=time(data_c), y=data_c)) +
            xlab("Time") + ylab("New Cases") + 
            xlim(as.Date('2020-3-1'), xmax_c) +
            ylim(c(0, max(data_c, window(data_c_fit, end=xmax_c)))) +
            plot.formatos
        gMdiff <- ggplot(data_m_fit) +
            geom_line(aes(x=time(data_m_fit), y=data_m_fit)) +
            geom_point(data=data_m, aes(x=time(data_m), y=data_m)) +
            xlab("Time") + ylab("New Deaths")+
            xlim(as.Date('2020-3-1'), xmax_m) +
            ylim(c(0, max(data_m, window(data_m_fit, end=xmax_m)))) +
            plot.formatos
        
        if(rollwindow > 1){
            
            zoo.rowmeans <- merge.zoo(cases = rollmean(now.covid.zoo, rollwindow), 
                                      cases_diff = rollmean(data_c, rollwindow), 
                                      deaths = rollmean(now.obito.covid.zoo, rollwindow),
                                      deaths_diff = rollmean(data_m, rollwindow))
            
            gCdiff  <- gCdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=cases_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            gMdiff  <- gMdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=deaths_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            
        }
    }
    else if(flag == "srag"){
        
        if(!is.null(max.date)){
            xmax_c=as.Date(max.date)
            xmax_m=as.Date(max.date)
            xmax_gC=as.Date(max.date)
            xmax_gM=as.Date(max.date)
        }
        else{
            xmax_c=max(time(data_c))
            xmax_m=max(time(data_m))
            xmax_gC=max(time(now.srag.zoo))
            xmax_gM=max(time(now.obito.srag.zoo))
        }
        gC <- ggplot(res) +
            geom_line(aes(x=Index, y=C)) +
            geom_point(data=now.srag.zoo, aes(x=Index, y=now.srag.zoo)) +
            xlim(c(as.Date('2020-3-1'), xmax_gC)) +
            ylim(c(0, max(now.srag.zoo, window(res$C, end=xmax_gC)))) +
            plot.formatos
        gM <- ggplot(res) +
            geom_line(aes(x=Index, y=CM)) +
            geom_point(data=now.obito.srag.zoo, aes(x=time(now.obito.srag.zoo),
                                                    y=now.obito.srag.zoo)) +
            xlim(as.Date('2020-3-1'), xmax_gM) +
            ylim(c(0, max(now.obito.srag.zoo, window(res$CM, end=xmax_gM)))) +
            plot.formatos
        
        data_c = diff(now.srag.zoo)
        data_m = diff(now.obito.srag.zoo)
        data_c_fit = diff(res$C)
        data_m_fit = diff(res$CM)
        
        if(weekly) {
            model_c_xts <- xts::apply.weekly(data_c_fit, sum)
            model_cm_xts <- xts::apply.weekly(data_m_fit, sum)
            data_c_fit <- zoo(model_c_xts, as.Date(index(model_c_xts)))
            data_m_fit <- zoo(model_cm_xts, as.Date(index(model_cm_xts)))
        }
        
        gCdiff <- ggplot(data_c_fit) +
            geom_line(aes(x=time(data_c_fit), y=data_c_fit)) +
            geom_point(data=data_c, aes(x=time(data_c), y=data_c)) +
            xlab("Time") + ylab("New Cases") + 
            xlim(as.Date('2020-3-1'), xmax_c) +
            ylim(c(0, max(data_c, window(data_c_fit, end=xmax_c)))) +
            plot.formatos
        gMdiff <- ggplot(data_m_fit) +
            geom_line(aes(x=time(data_m_fit), y=data_m_fit)) +
            geom_point(data=data_m, aes(x=time(data_m), y=data_m)) +
            xlab("Time") + ylab("New Deaths") + 
            xlim(as.Date('2020-3-1'), xmax_m) +
            ylim(c(0, max(data_m, window(data_m_fit, end=xmax_m)))) +
            plot.formatos
        
        if(rollwindow > 1){
            
            zoo.rowmeans <- merge.zoo(cases = rollmean(now.srag.zoo, rollwindow), 
                                      cases_diff = rollmean(data_c, rollwindow), 
                                      deaths = rollmean(now.obito.srag.zoo, rollwindow),
                                      deaths_diff = rollmean(data_m, rollwindow))
            
            gCdiff  <- gCdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=cases_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            gMdiff  <- gMdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=deaths_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            
        }    
        
    }
    grid.arrange(gC, gM, gCdiff, gMdiff)#,
    #top=paste0("casos ", flag))
}

###################################
################################### Only new cases and deaths:
###################################

test_modifying_parameters_new <- function(parameters, interventions, Y, params = NA,
                                      flag = "covid", base_size=18, rollwindow = 1,
                                      rowmean_line_color = "grey",rowmean_line_size = 1.2,
                                      weekly = TRUE, max.date=NULL){
    plot.formatos <- theme_bw(base_size=base_size)
    for (i in names(params)){parameters[i] <- params[i]}
    
    out <- run_covid19(parameters, interventions, Y)
    res <- result_ts(out[[2]], parameters, parameters["startdate"], diff=F)
    
    
    ##res <- result_ts(out[[2]], parameters, date_parameters, diff=F)
    
    if(flag == "covid"){
        gC <- ggplot(res) +
            geom_line(aes(x=Index, y=C)) +
            geom_point(data=now.covid.zoo, aes(x=Index, y=now.covid.zoo)) +
            xlim(as.Date('2020-3-1'), max(time(now.covid.zoo))) +
            ylim(c(0, max(now.covid.zoo, window(res$C, end=max(time(now.covid.zoo)))))) +
            plot.formatos
        gM <- ggplot(res) +
            geom_line(aes(x=Index, y=CM)) +
            geom_point(data=now.obito.covid.zoo, aes(x=time(now.obito.covid.zoo),
                                                     y=now.obito.covid.zoo)) +
            xlim(as.Date('2020-3-1'), max(time(now.obito.covid.zoo))) +
            ylim(c(0, max(now.obito.covid.zoo, window(res$CM, end=max(time(now.obito.covid.zoo)))))) +
            plot.formatos
        
        data_c = diff(now.covid.zoo)
        data_m = diff(now.obito.covid.zoo)
        data_c_fit = diff(res$C)
        data_m_fit = diff(res$CM)
        
        if(!is.null(max.date)){
            xmax_c=as.Date(max.date)
            xmax_m=as.Date(max.date)
        }
        else{
            xmax_c=max(time(data_c))
            xmax_m=max(time(data_m))
        }
        
        if(weekly) {
            model_c_xts <- xts::apply.weekly(data_c_fit, sum)
            model_cm_xts <- xts::apply.weekly(data_m_fit, sum)
            data_c_fit <- zoo(model_c_xts, as.Date(index(model_c_xts)))
            data_m_fit <- zoo(model_cm_xts, as.Date(index(model_cm_xts)))
        }
        
        gCdiff <- ggplot(data_c_fit) +
            geom_line(aes(x=time(data_c_fit), y=data_c_fit)) +
            geom_point(data=data_c, aes(x=time(data_c), y=data_c)) +
            xlab("Time") + ylab("New Cases") + 
            xlim(as.Date('2020-3-1'), xmax_c) +
            ylim(c(0, max(data_c, window(data_c_fit, end=xmax_c)))) +
            plot.formatos
        gMdiff <- ggplot(data_m_fit) +
            geom_line(aes(x=time(data_m_fit), y=data_m_fit)) +
            geom_point(data=data_m, aes(x=time(data_m), y=data_m)) +
            xlab("Time") + ylab("New Deaths")+
            xlim(as.Date('2020-3-1'), xmax_m) +
            ylim(c(0, max(data_m, window(data_m_fit, end=xmax_m)))) +
            plot.formatos
        
        if(rollwindow > 1){
            
            zoo.rowmeans <- merge.zoo(cases = rollmean(now.covid.zoo, rollwindow), 
                                      cases_diff = rollmean(data_c, rollwindow), 
                                      deaths = rollmean(now.obito.covid.zoo, rollwindow),
                                      deaths_diff = rollmean(data_m, rollwindow))
            
            gCdiff  <- gCdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=cases_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            gMdiff  <- gMdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=deaths_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            
        }
    }
    else if(flag == "srag"){
        
        if(!is.null(max.date)){
            xmax_c=as.Date(max.date)
            xmax_m=as.Date(max.date)
            xmax_gC=as.Date(max.date)
            xmax_gM=as.Date(max.date)
        }
        else{
            xmax_c=max(time(data_c))
            xmax_m=max(time(data_m))
            xmax_gC=max(time(now.srag.zoo))
            xmax_gM=max(time(now.obito.srag.zoo))
        }
        gC <- ggplot(res) +
            geom_line(aes(x=Index, y=C)) +
            geom_point(data=now.srag.zoo, aes(x=Index, y=now.srag.zoo)) +
            xlim(c(as.Date('2020-3-1'), xmax_gC)) +
            ylim(c(0, max(now.srag.zoo, window(res$C, end=xmax_gC)))) +
            plot.formatos
        gM <- ggplot(res) +
            geom_line(aes(x=Index, y=CM)) +
            geom_point(data=now.obito.srag.zoo, aes(x=time(now.obito.srag.zoo),
                                                    y=now.obito.srag.zoo)) +
            xlim(as.Date('2020-3-1'), xmax_gM) +
            ylim(c(0, max(now.obito.srag.zoo, window(res$CM, end=xmax_gM)))) +
            plot.formatos
        
        data_c = diff(now.srag.zoo)
        data_m = diff(now.obito.srag.zoo)
        data_c_fit = diff(res$C)
        data_m_fit = diff(res$CM)
        
        if(weekly) {
            model_c_xts <- xts::apply.weekly(data_c_fit, sum)
            model_cm_xts <- xts::apply.weekly(data_m_fit, sum)
            data_c_fit <- zoo(model_c_xts, as.Date(index(model_c_xts)))
            data_m_fit <- zoo(model_cm_xts, as.Date(index(model_cm_xts)))
        }
        
        gCdiff <- ggplot(data_c_fit) +
            geom_line(aes(x=time(data_c_fit), y=data_c_fit)) +
            geom_point(data=data_c, aes(x=time(data_c), y=data_c)) +
            xlab("Time") + ylab("New Cases") + 
            xlim(as.Date('2020-3-1'), xmax_c) +
            ylim(c(0, max(data_c, window(data_c_fit, end=xmax_c)))) +
            plot.formatos
        gMdiff <- ggplot(data_m_fit) +
            geom_line(aes(x=time(data_m_fit), y=data_m_fit)) +
            geom_point(data=data_m, aes(x=time(data_m), y=data_m)) +
            xlab("Time") + ylab("New Deaths") + 
            xlim(as.Date('2020-3-1'), xmax_m) +
            ylim(c(0, max(data_m, window(data_m_fit, end=xmax_m)))) +
            plot.formatos
        
        if(rollwindow > 1){
            
            zoo.rowmeans <- merge.zoo(cases = rollmean(now.srag.zoo, rollwindow), 
                                      cases_diff = rollmean(data_c, rollwindow), 
                                      deaths = rollmean(now.obito.srag.zoo, rollwindow),
                                      deaths_diff = rollmean(data_m, rollwindow))
            
            gCdiff  <- gCdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=cases_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            gMdiff  <- gMdiff + geom_line(data=zoo.rowmeans, aes(x=Index, y=deaths_diff), 
                                          colour = rowmean_line_color, size = rowmean_line_size) 
            
        }    
        
    }
    grid.arrange(gCdiff, gMdiff)#,
    #top=paste0("casos ", flag))
}