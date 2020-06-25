library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#FUnction to Obtain Betas for the Model 
combbeta <- function(scen, time, tstart, t_dur, R0Dec) {
  gamma <- 1/GenTime(3.3, 2.8)
  if(scen == 0) {
    output <-  1.7*gamma
  }
  if(scen == 1) {
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
           R0Dec*gamma, 1.7*gamma)
  }
  if(scen == 2) {
    R0lin <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(R0Dec, 1.7), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
            R0lin(time)*gamma, 1.7*gamma)
  }
  if(scen == 3) {
    R0lin <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(1.7, R0Dec), method="linear", rule =2) 
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
           R0lin(time)*gamma, 1.7*gamma)
  }
  if(scen == 4) {
    R0incdec <-approxfun(x=c(tstart, tstart+(t_dur/2), tstart+t_dur), y= c(1.7, R0Dec, 1.7), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
           R0incdec(time)*gamma, 1.7*gamma)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart+(t_dur*0.0833) & time <= tstart+(t_dur*0.25)) | (time >= tstart+(t_dur*0.4167) & time <= tstart+(t_dur*0.5833)) |
             (time >= tstart+(t_dur*0.75 ) & time <= tstart+(t_dur*0.9167)), 
           R0Dec*gamma, 1.7*gamma)
  }
  return(output)
}

plot(seq(0,365),combbeta(5, seq(0,365), 71, 12*7, 0.8))

#SEIR set of ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, R0Dec)
    
    dS = - beta*S*(E+I)
    dE = beta*S*(E+I) - theta*E
    dI = theta*E - gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dE, dI, dR, dC)))
  })
} 

#### Baseline Model #### 

init <- c(S = 0.9999, E = 0,I = 0.0001, R = 0, C = 0)
times <- seq(0,300,by = 1)
parms = c(gamma = 1/GenTime(3.3, 2.8),
          theta = 1/5,
          scen = 1,
          tstart = 51,
          t_dur = 12*7,
          R0Dec = 0.5)
outbase <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

datalist <- list()

for(j in 1:length(seq(1,5))) {
  datalist[[j]] <- local({
    j=j
    data <- data.frame(matrix(nrow = 9, ncol = 0))
    explor_scen <- c(0, seq(1,5)[j])
    for(i in 1:2) {
      parms["scen"] <- explor_scen[i]
      out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], "r0" = combbeta(explor_scen[i], times, parms[["tstart"]], parms[["t_dur"]], parms[["R0Dec"]])/parms[["gamma"]])
      out$re <- out$r0*out$S
      data <- rbind(data, out)
    }
    
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    plotr0 <- melt(data, id.vars = c("time", "group"), measure.vars = ("r0"))
    plotre <- melt(data, id.vars = c("time", "group"), measure.vars = ("re"))
    
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      labs(x ="Time (Days)", y = "Prevalence", col = "") + scale_y_continuous(limits = c(0 , 0.25),expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) + scale_x_continuous( expand = c(0, 0)) + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkred", "darkred"))

    shade <- data.frame(xmin =  parms[["tstart"]], 
                        xmax = parms[["tstart"]]+(parms[["t_dur"]]),
                        ymin = 0, ymax = Inf)
    
    if(parms[["scen"]]  == 1) {
     p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                          fill = "darkblue") + geom_line(size = 1.1, stat = "identity")
    }
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity")
    }
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity")
    }
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity")
    }
    if(parms[["scen"]]  == 5) {
      
      shade <- data.frame(xmin =  c(parms[["tstart"]]+(parms[["t_dur"]]*0.0833), parms[["tstart"]]+(parms[["t_dur"]]*0.4167), 
                             parms[["tstart"]]+(parms[["t_dur"]]*0.75)),
                   xmax = c(parms[["tstart"]]+(parms[["t_dur"]]*0.25), parms[["tstart"]]+(parms[["t_dur"]]*0.5833), 
                            parms[["tstart"]]+(parms[["t_dur"]]*0.9167)),
                   ymin = 0, ymax = Inf)
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin,  ymax = ymax, 
                     xmin = xmin, xmax = xmax), alpha = 0.2, fill = "darkblue") + geom_line(size = 1.1, stat = "identity")
    }
    
    p2 <- ggplot(plotr0, aes(x = time, y = value, col = group, alpha= group))  + theme_bw() +
      labs(x ="Time (Days)", y = expression(R[0]), col = "") + scale_y_continuous(limits = c(0 , 1.75), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_line(size = 1.1, stat = "identity") + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    
    p3 <- ggplot(plotre, aes(x = time, y = value, col = group, alpha= group)) + theme_bw() +
      labs(x ="Time (Days)", y = expression(R[e]), col = "") + scale_y_continuous(limits = c(0 , 1.75), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_hline(yintercept = 1, size = 1.1, lty = 2, col = "black") + geom_line(size = 1.02, stat = "identity") + 
     scale_alpha_manual(values = c(0.35,1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.5, 0.5))
    return(combplot)
  })
}

combplot <- ggarrange(datalist[[1]], datalist[[2]], datalist[[3]], datalist[[4]], datalist[[5]],nrow = 1, ncol = 5, 
          legend = "none", align = "hv", labels =  c("A","B","C","D","E") ,font.label = c(size = 20))

ggsave(combplot, filename = "5_scenarios.png", dpi = 300, type = "cairo", width = 20, height = 10, units = "in")


