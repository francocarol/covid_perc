if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(ggcorrplot)){install.packages("ggcorrplot"); library(ggcorrplot)}
if(!require(RColorBrewer)){install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(gridExtra)){install.packages("gridExtra"); library(gridExtra)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}
if(!require(zoo)){install.packages("zoo"); library(zoo)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}

plot_scenarios_simple <- function(res, scn.labels, cols=c('C', 'H', 'U', 'M'),
                                  base_size=18){
    nscn <- ncol(res)/(1+length(cols))
    color_seq <- colorRampPalette(brewer.pal(8, "Accent"))(nscn)
    ylabels <- c("Symptomatic Cases", "Hospitalized Cases", "ICU Cases",
                 "Cumulative Deaths")
    plot.formatos <- theme_bw(base_size=base_size)

    plots <- list()
    for (label in 1:length(cols)){
        plots[[label]] <- ggplot(data=res)
        for (i in 1:nscn){
            plots[[label]] <- plots[[label]] +
                geom_line(aes_(x=as.name("times"),
                               y=as.name(paste0(cols[label], i)),
                               color=factor(scn.labels[i])),
                          size=1, show.legend=TRUE)
        }
        plots[[label]] <- plots[[label]] +
            #scale_y_log10(limits=c(10,NA)) +
            ylab(ylabels[label]) +
            scale_color_manual(values = color_seq, name = "Cenário", labels = scn.labels) +
            theme(legend.position = "right", legend.box = "vertical") +
            plot.formatos
    }
    return(plots)
}

plot_scenarios_simple_zoo <- function(res, scn.labels, logscale=F,
                                      cols=c('C', 'H', 'U', 'M'),
                                      ylabels=c("Total Observed Cases",
                                                "Hospitalized Cases", "ICU
                                                Cases", "Cumulative Deaths"),
                                      data=NULL, hlines=c(), base_size=18,
                                      lang='pt', legend.pos='bottom',
                                      label_months=F, line.size=1, color.palette="Accent",
                                      start.data=as.Date('2020-3-1'), end.data=NULL){
    nscn <- ncol(res)/length(cols)
    color_seq <- colorRampPalette(brewer.pal(8, color.palette))(nscn)
    plot.formatos <- theme_bw(base_size=base_size)
    if (lang=='pt') {
        legend.title <- "Cenários"
        x.label <- "Data"
    }
    if (lang=='en') {
        legend.title <- "Scenarios"
        x.label <- "Date"
    }
    plots <- list()
    for (label in 1:length(cols)){
        plots[[label]] <- ggplot(data=res)
        for (i in 1:nscn){
            plots[[label]] <- plots[[label]] +
                geom_line(aes_(x=~Index,
                               y=as.name(paste0(cols[label], i)),
                               color=as.character(i)),
                          size=line.size, show.legend=TRUE)
        }
        # add hlines
        if(cols[label] %in% names(hlines)){
            plots[[label]] <- plots[[label]] +
                geom_hline(yintercept=hlines[cols[label]], linetype="dashed",
                           color="black")
        }
        if (legend.pos=="right")
            plots[[label]] <- plots[[label]] +
                theme(legend.position = "right", legend.box = "vertical")
        if (legend.pos=="bottom")
            plots[[label]] <- plots[[label]] +
                theme(legend.position = "bottom", legend.box = "horizontal")

        plots[[label]] <- plots[[label]] +
            ylab(ylabels[label]) + xlab(x.label) +
            scale_color_manual(values = color_seq, name = legend.title, labels = scn.labels) +
            plot.formatos
        if (logscale)
            plots[[label]] <- plots[[label]] + scale_y_log10(limits=c(10,NA))
        if(label_months)
            plots[[label]] <- plots[[label]] +
                scale_x_date(date_breaks = "months" , date_labels = "%b") +
                theme(axis.text.x=element_text(angle=90, hjust=1))

        if(!is.null(data)){
            # TODO: include points in the legend (shape=? + scale_shape_manual))
            if(data=='srag'){
                if (!is.null(end.data)){
                    data_c  <- window(now.srag.zoo, start=start.data, end=end.data)
                    data_cm <- window(now.obito.srag.zoo, start=start.data, end=end.data)
                }
                else{
                    data_c  <- now.srag.zoo
                    data_cm <- now.obito.srag.zoo
                }
            } else if(data=='covid'){
                data_c = now.covid.zoo
                data_cm = now.obito.covid.zoo
            }
            if(cols[label] == 'C'){
                plots[[label]] <- plots[[label]] +
                    geom_point(data=fortify(data_c, melt=T),
                               aes(x=Index, y=Value), size=0.7)
            } else if(cols[label] == 'M') {
                plots[[label]] <- plots[[label]] +
                    geom_point(data=fortify(data_cm, melt=T),
                               aes(x=Index, y=Value), size=0.7)
            }
        }
    }
    return(plots)
}

plot_scenarios_sensitivity <- function(QQ, scn.labels, cols=c('C', 'H', 'U', 'M'),
                                       ylabels=c("Symptomatic Cases",
                                                 "Hospitalized Cases",
                                                 "ICU Cases",
                                                 "Cumulative Deaths"),
                                       quantiles=c(0.025, 0.5, 0.975),
                                       data=NULL, hlines=c(), logscale=FALSE,
                                       base_size=18, lang='pt',
                                       legend.pos='bottom', label_months=F){
    nscn <- length(QQ)
    if(missing(scn.labels))
        scn.labels <- 1:nscn
    QQ[["suffixes"]] <- 1:nscn
    QQ <- do.call(cbind.zoo, QQ)
    color_seq <- colorRampPalette(brewer.pal(8, "Accent"))(nscn)
    plot.formatos <- theme_bw(base_size=base_size)
    if (lang=='pt')
        legend.title <- "Cenários"
    if (lang=='en')
        legend.title <- "Scenarios"

    plots <- list()
    for (label in 1:length(cols)){
        plots[[label]] <- ggplot(data=QQ)
        for (i in 1:nscn){
            plots[[label]] <- plots[[label]] +
                geom_line(aes_(x=~Index,
                               y=as.name(paste(cols[label], 100*quantiles[2], i, sep='.')),
                               color=as.character(i)), show.legend=TRUE) +
                geom_ribbon(aes_(x=~Index,
                                 ymin=as.name(paste(cols[label], 100*quantiles[1], i, sep='.')),
                                 ymax=as.name(paste(cols[label], 100*quantiles[3], i, sep='.'))),
                            fill=color_seq[i], alpha=0.3)
        }
        # add hlines
        if(cols[label] %in% names(hlines)){
            plots[[label]] <- plots[[label]] +
                geom_hline(yintercept=hlines[cols[label]], linetype="dashed",
                           color="black")
        }
        if (legend.pos=="right")
            plots[[label]] <- plots[[label]] +
                theme(legend.position = "right", legend.box = "vertical")
        if (legend.pos=="bottom")
            plots[[label]] <- plots[[label]] +
                theme(legend.position = "bottom", legend.box = "horizontal")
        if (logscale)
            plots[[label]] <- plots[[label]] + scale_y_log10(limits=c(10,NA))
        if(label_months)
            plots[[label]] <- plots[[label]] +
                scale_x_date(date_breaks = "months" , date_labels = "%b") +
                theme(axis.text.x=element_text(angle=90, hjust=1))

        plots[[label]] <- plots[[label]] +
            ylab(ylabels[label]) + xlab('date') +
            scale_color_manual(values = color_seq, name = legend.title, labels = scn.labels) +
            scale_x_date(date_labels = "%d/%b", name="") +
            plot.formatos

        if(!is.null(data)){
            # TODO: include points in the legend (shape=? + scale_shape_manual))
            if(data=='srag'){
                data_c <- now.srag.zoo
                data_cm <- now.obito.srag.zoo
            } else if(data=='covid'){
                data_c = diff(now.covid.zoo)
                data_cm = diff(now.obito.covid.zoo)
            }
            if(cols[label] == 'C'){
                plots[[label]] <- plots[[label]] +
                    geom_point(data=fortify(data_c, melt=T),
                               aes(x=Index, y=Value))
            } else if(cols[label] == 'M') {
                plots[[label]] <- plots[[label]] +
                    geom_point(data=fortify(data_cm, melt=T),
                               aes(x=Index, y=Value))
            }
        }
    }
    return(plots)
}

plot.hist.posteriors <- function(posterior){
    post.gathered <- as.data.frame(abc.out$unadj.values) %>% gather()
    p <- ggplot(post.gathered, aes(value)) +
        geom_histogram(bins=8) +
        facet_wrap(~key, scales='free_x')
}

plot.cor.posteriors <- function(posterior){
    corr<-cor(posterior)
    p <- ggcorrplot(corr, hc.order = TRUE, type = "lower", lab=T,
                    outline.col = "white")
}

# arrange plots in a grid using POG
combine.plots <- function(plots, legend.pos='bottom', margin=0.3){
    if(legend.pos == 'bottom')
        return((((plots[[1]] + theme(legend.position = "bottom", legend.direction="vertical") | plots[[2]] + theme(legend.position = "none")| plots[[3]] + theme(legend.position = "none")) / (plots[[4]]  + theme(legend.position = "none")| plots[[5]] + theme(legend.position = "none") | plots[[6]] + theme(legend.position = "none"))) / guide_area()) + plot_layout(guide = "collect", heights = c(1, 1, margin)))
    if(legend.pos == 'right')
        return(((plots[[1]] + theme(legend.position = "right") | plots[[2]] + theme(legend.position = "none")| plots[[3]] + theme(legend.position = "none")) / (plots[[4]]  + theme(legend.position = "none")| plots[[5]] + theme(legend.position = "none") | plots[[6]] + theme(legend.position = "none"))) + plot_layout(guide = "collect"))
}

combine.plots.R_eff<-function(plots,legend.pos='bottom',margin=0.3){
    if(legend.pos == 'bottom')
        return((((plots[[7]] + theme(legend.position = "bottom", legend.direction="vertical")| plots[[8]] + theme(legend.position = "none")) / guide_area()) + plot_layout(guide = "collect", heights = c(1, 1, margin))))
    if(legend.pos == 'right')
        return(((plots[[7]] + theme(legend.position = "right") | plots[[8]] + theme(legend.position = "none"))))
}


### Plot results by age group

plot.ages <- function(data, index, age_group = 1:19, plot.label = TRUE, palette = viridis(19), n.breaks = 18,
                      margins = c(4,4,2,5), popstruc = NA, control_pop = FALSE, title = "", y_lim = NA,
                      start_date = "2020-02-19", lang = "pt", ylabel_age = NA, cex.lab = 1, warn = TRUE,
                      diff = FALSE, x_lim = NA) { 
    #To-do: start_date, automatizar 
    date_parameters = c(startdate = as.Date(start_date))
    date_sequence = as.Date(date_parameters["startdate"]) + data[,1]
    
    nindex = length(index)/19
    
    # Somar grupos de idades caso mais de um compartimento seja escolhido
    if(nindex > 1) {
        sumdata = matrix(NA, nrow(data), 19)
        for(j in 1:19) sumdata[,j] = rowSums(data[,index[j+(0:(nindex-1))*19]+1])
        if(diff) subdata <- apply(subdata, 2, diff)
        z <- zoo(sumdata, date_sequence)
        colnames(z) = ((1:19)*5)-5
    } else {
    w = index+1
    subdata <- data[,w]
    if(diff) subdata <- apply(subdata, 2, diff)
    z <- zoo(subdata, date_sequence)
    colnames(z) = ((1:19)*5)-5
    }
    
    if(lang == "pt") {
        if(is.na(ylabel_age)){
        if(control_pop == FALSE) ylabel_age = "Número acumulado de mortes"
        if(control_pop == TRUE)  ylabel_age = "Proporção de mortes acumuladas por faixa etária"}
        xlabel_age = "Data"}
    
    if(lang == "en") {
        if(is.na(ylabel_age)){
        if(control_pop == FALSE) ylabel_age = "Cumulative number of deaths"
        if(control_pop == TRUE)  ylabel_age = "Proportion of cumulative deaths according to age group"}
        xlabel_age = "Date"}
    
    
    if(all(control_pop == FALSE & is.na(y_lim))) y_lim = c(0,max(z[,age_group]))
    if(all(control_pop == TRUE & is.na(y_lim))) y_lim = c(0,max( t(z[,age_group])/popstruc$pop[age_group] ))
    
    # To do: adjust so it calculates the date ranges
    month_seq = seq(as.Date("2020-03-01") ,as.Date("2021-06-01") , by = "months")
    if(is.na(x_lim)) {x_lim = date_sequence[c(1,length(date_sequence))]} else {x_lim = as.Date(x_lim)}
    
    par(mar = margins)
    plot(z[,1], type = "n", ylim = y_lim, 
         xlim = x_lim,
         xlab = xlabel_age, ylab = ylabel_age, main = title,
         cex.lab = cex.lab, xaxt = "n")
    axis(1, at = month_seq, labels = format(month_seq, "%b"), las = 2)
    
    if(control_pop == FALSE) for(i in age_group) lines(z[,i], col = palette[i])
    if(control_pop == TRUE)  for(i in age_group) lines(z[,i]/popstruc$pop[i], col = palette[i])
    # plot(z[,i]/popstruc$pop[i])
    # i = 10
    # popstruc$pop[i]
    # max(z[,i])
    # max(z[,i])/popstruc$pop[i]
    # plot(z[,i])
    
    if (plot.label) { 
        # Set up and plot labels
        if(control_pop == FALSE) finalcount <- t(z)[,nrow(z)]# t(z[nrow(z),])
        if(control_pop == TRUE) finalcount <- t(z)[age_group,nrow(z)]/popstruc$pop[age_group]
        
        if(n.breaks == 1) {
        labeldf <- data.frame(label = "0 ~ 90*", ypos = mean(finalcount))
        } else {
        dist <- hist(finalcount, breaks = seq(0,max(finalcount),length.out = n.breaks), plot = F)
        
        labeldf <- data.frame()
        for(i in 1:length(dist$breaks)){
            
            if(i < length(dist$breaks)) {
                invalue <- finalcount >= dist$breaks[i] & finalcount < dist$breaks[i+1] } else {
                    invalue <- finalcount >= dist$breaks[i]
                }
            
            if(any(invalue)==TRUE){
                
                if (sum(invalue) > 3){
                    if(all(invalue[c(min(which(invalue)):max(which(invalue)))])) {
                        age_label = paste0(names(invalue)[min(which(invalue))], " - ",  names(invalue)[max(which(invalue))])
                    } else {
                        age_label = paste0(names(invalue)[min(which(invalue))], " ~ ",  names(invalue)[max(which(invalue))], "*")
                        if(warn == TRUE) warning("* Not all age groups are represented within this range")
                    }
                    
                } else {  
                    age_label = paste(names(invalue)[invalue], collapse = ", ")
                }
                
                tempdf <- data.frame(label = age_label, ypos = mean(finalcount[invalue]))
                labeldf <- rbind(labeldf, tempdf)
            }}}
        
        for(i in 1:nrow(labeldf)) mtext(labeldf$label[i], side = 4, line = 0.5, at = labeldf$ypos[i], cex = 0.7, las = 1)
    }
    
}

# Plot number of cases according to age groups similar to SP prevalence studies

# data = out2
# index = Cindex
# ylabel_age = "Proporção de casos acumulados por faixa etária"; warn = F; control_pop = T; y_lim = c(0,1); cex.lab = 0.9
# margins = c(4,4,2,3); palette = viridis(6); title = ""; start_date = "2020-02-19"; lin.col = "red"; lin.wd = 1 ;by.age = T

plot.prev <- function(data, index, palette = viridis(6), 
                      margins = c(4,4,2,5), popstruc = NA, control_pop = FALSE, title = "",
                      start_date = "2020-02-19", lang = "pt", ylabel_age = NA, diff = FALSE,
                      cex.lab = 1, warn = TRUE, y_lim = NA, lin.col = "red", lin.wd = 1,
                      by.age = TRUE, return.max = FALSE, x_lim = NA) { 
    #To-do: start_date, automatizar
    date_parameters = c(startdate = as.Date(start_date))
    date_sequence = as.Date(date_parameters["startdate"]) + data[,1]
    
    nindex = length(index)/19
    
    # Somar grupos de idades caso mais de um compartimento seja escolhido
    if(nindex > 1) {
        sumdata = matrix(NA, nrow(data), 19)
        for(j in 1:19) sumdata[,j] = rowSums(data[,index[j+(0:(nindex-1))*19]+1])
        sumdata2 <- data.frame(a0_19 = rowSums(sumdata[,1:4]),
                               a20_29 = rowSums(sumdata[,5:6]),
                               a30_39 = rowSums(sumdata[,7:8]),
                               a40_49 = rowSums(sumdata[,9:10]),
                               a50_59 = rowSums(sumdata[,11:12]),
                               a60_99 = rowSums(sumdata[,13:19]))
        if(diff == TRUE) subdata2 <- apply(subdata2, 2, diff)
        z <- zoo(sumdata2, date_sequence)
    } else {
        w = index+1
        subdata <- data[,w]
        subdata2 <- data.frame(a0_19 = rowSums(subdata[,1:4]),
                               a20_29 = rowSums(subdata[,5:6]),
                               a30_39 = rowSums(subdata[,7:8]),
                               a40_49 = rowSums(subdata[,9:10]),
                               a50_59 = rowSums(subdata[,11:12]),
                               a60_99 = rowSums(subdata[,13:19]))
        if(diff == TRUE) subdata2 <- apply(subdata2, 2, diff)
        z <- zoo(subdata2, date_sequence)
    }
    
    if(lang == "pt") {
        if(is.na(ylabel_age)){
            if(control_pop == FALSE) ylabel_age = "Número acumulado de mortes"
            if(control_pop == TRUE)  ylabel_age = "Proporção de mortes acumuladas por faixa etária"}
        xlabel_age = "Data"}
    
    if(lang == "en") {
        if(is.na(ylabel_age)){
            if(control_pop == FALSE) ylabel_age = "Cumulative number of deaths"
            if(control_pop == TRUE)  ylabel_age = "Proportion of cumulative deaths according to age group"}
        xlabel_age = "Date"}
    
    popstruc2 <- c(sum(popstruc$pop[1:4]), sum(popstruc$pop[5:6]), sum(popstruc$pop[7:8]), sum(popstruc$pop[9:10]),
                   sum(popstruc$pop[11:12]),sum(popstruc$pop[13:19]))
    #  
    if(all(control_pop == FALSE & is.na(y_lim))) y_lim = c(0,max(z[,]))
    if(all(control_pop == TRUE & is.na(y_lim))) y_lim = c(0,max( t(z[,])/popstruc2))
    
    # if(control_pop == TRUE) {
    #     z$all <- rowSums(z)/sum(popstruc)
    # }

    ## To do: adjust to adapt according to data provided        
    month_seq = seq(as.Date("2020-03-01") ,as.Date("2021-06-01") , by = "months")
    if(is.na(x_lim)) {x_lim = date_sequence[c(1,length(date_sequence))]} else {x_lim = as.Date(x_lim)}
    
    par(mar = margins)
    plot(z[,1], type = "n", ylim = y_lim, 
         xlim = x_lim,
         xlab = xlabel_age, ylab = ylabel_age, main = title,
         cex.lab = cex.lab, xaxt = "n")
    axis(1, at = month_seq, labels = format(month_seq, "%b"), las = 2)
    if(control_pop == FALSE) for(i in 1:6) lines(z[,i], col = palette[i])
    if(control_pop == TRUE)  {
        if(by.age == TRUE) {for(i in 1:6) lines(z[,i]/popstruc2[i], col = palette[i]) } else {
        lines(z$all, col = lin.col, lwd = lin.wd)}
    }
    if(return.max) {
        whichmax <- function(column){
            which(column==max(column))
        }
        
        peak = index(z)[apply(z, 2, whichmax)]
        names(peak) = c("a0_19",  "a20_29", "a30_39", "a40_49", "a50_59", "a60_99")
       return(peak)
    }
}


## Plot Reff for different scenarios

plot_scenarios_reff <- function(result, scenarios, color.palette = "Dark2", 
                                base_size = 18, lang = "pt", scn.labels){ 
    
    Re_cori_scn <- data.frame(result$reffs[[1]]$Reff_cori_sample, 
                              Index = index(result$reffs[[1]]$Reff_cori_sample), 
                              scenario = scenarios[1])
    
    for(i in 2:length(result)) {
        Re_cori_scn <- rbind(Re_cori_scn,
                             data.frame(result$reffs[[i]]$Reff_cori_sample,
                                        Index = index(result$reffs[[i]]$Reff_cori_sample),
                                        scenario = scenarios[i]))    
    }
    
    nscn <- length(result)
    color_seq <- colorRampPalette(brewer.pal(8, color.palette))(nscn)
    
    if(missing(scn.labels)) scn.labels <- 1:nscn
    
    plot.formatos <- theme_bw(base_size=base_size)
    if (lang=='pt') {
        scn.labels <- names(scn.ALL.PT)[ match(scenarios,scn.ALL.PT) ]
        legend.title <- "Cenários"
        x.label <- "Data"
        y.label <- "Número reprodutivo efetivo"
    }
    if (lang=='en') {
        scn.labels <- names(scn.ALL)[ match(scenarios,scn.ALL) ]
        legend.title <- "Scenarios"
        x.label <- "Date"
        y.label <- "Effective reproduction number"
    }
    
    ggplot.Rt <- 
        ggplot(data = Re_cori_scn, aes(x=Index,y=Mean, group = scenario)) +
        geom_ribbon(aes(ymin = as.numeric(Lower), ymax = as.numeric(Upper), fill=factor(scenario)), 
                    alpha=0.3) +
        geom_line(aes(color = factor(scenario))) +
        scale_fill_manual(values = color_seq,  name = legend.title, labels = scn.labels) +
        scale_color_manual(values = color_seq) +
        geom_hline(yintercept=1, linetype="dashed", col="black") +
        scale_x_date( date_labels = "%d/%b", name=x.label) +
        ylab(y.label) + guides(color=FALSE)
    
    ggplot.Rt <- ggplot.Rt +
        theme(legend.position = "bottom", 
              legend.box = "horizontal",
              legend.direction = "vertical")
    
    return(ggplot.Rt)
}

