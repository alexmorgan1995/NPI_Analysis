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
    output <- ifelse((time >= tstart & time <= tstart+(t_dur*0.16667)) | (time >= tstart+(t_dur*0.3333) & time <= tstart+(t_dur*0.5)) |
             (time >= tstart+(t_dur*0.6667) & time <= tstart+(t_dur*0.83333)), 
           R0Dec*gamma, 1.7*gamma)
  }
  return(output)
}

plot(seq(0,365),combbeta(5, seq(0,365), 71, 12*7, 0.8))

#SEIR set of ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, R0Dec)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 


# Baseline Model Analysis - 5 Scenarios -----------------------------------

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 0,
          tstart = 71,
          t_dur = 12*7,
          R0Dec = 0.8)

datalist <- list()

for(j in 1:length(seq(1,5))) {
  datalist[[j]] <- local({
    j=j
    data <- data.frame(matrix(nrow = 9, ncol = 0))
    explor_scen <- c(0, seq(1,5)[j])
    
 
    for(i in 1:2) {
      
      parms["scen"] <- explor_scen[i]
      
      if(parms["scen"] != 0 && parms["scen"] != 1) {
        parms["t_dur"] = 24*7
      }
      
      out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], "r0" = combbeta(explor_scen[i], times, parms[["tstart"]], parms[["t_dur"]], parms[["R0Dec"]])/parms[["gamma"]])
      out$re <- out$r0*out$S
      data <- rbind(data, out)
    }
 
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    plotr0 <- melt(data, id.vars = c("time", "group"), measure.vars = ("r0"))
    plotre <- melt(data, id.vars = c("time", "group"), measure.vars = ("re"))
    
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.125),expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=18), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
      scale_x_continuous( expand = c(0, 0)) + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkred", "darkred"))

    p2 <- ggplot(plotr0, aes(x = time, y = value, col = group, alpha= group))  + theme_bw() +
      scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=18),axis.title.x = element_blank(), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_line(size = 1.1, stat = "identity") + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    p3 <- ggplot(plotre, aes(x = time, y = value, col = group, alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=18),axis.title.x = element_text(size=18), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_hline(yintercept = 1, size = 1.1, lty = 2, col = "black") + geom_line(size = 1.02, stat = "identity") + 
      scale_alpha_manual(values = c(0.35,1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    shade <- data.frame(xmin =  parms[["tstart"]], 
                        xmax = parms[["tstart"]]+(parms[["t_dur"]]),
                        ymin = 0, ymax = Inf)
    
    if(parms[["scen"]]  == 1) {
     p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                          fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", col = "")
     p2 <- p2 + labs(x ="", y = expression(R[0]), col = "")
     p3 <- p3 + labs(x ="Time (Days)", y = expression(R[e]), col = "")
    }
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 5) {
      
      shade <- data.frame(xmin =  c(parms[["tstart"]], parms[["tstart"]]+(parms[["t_dur"]]*0.333), 
                             parms[["tstart"]]+(parms[["t_dur"]]*0.667)),
                   xmax = c(parms[["tstart"]]+(parms[["t_dur"]]*0.1667), parms[["tstart"]]+(parms[["t_dur"]]*0.53), 
                            parms[["tstart"]]+(parms[["t_dur"]]*0.833)),
                   ymin = 0, ymax = Inf)
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin,  ymax = ymax, 
                     xmin = xmin, xmax = xmax), alpha = 0.2, fill = "darkblue") + geom_line(size = 1.1, stat = "identity") +
      labs(x ="", y = "", col = "")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.45, 0.5))
    return(combplot)
  })
}

combplot <- ggarrange(datalist[[1]], datalist[[2]], datalist[[3]], datalist[[4]], datalist[[5]],nrow = 1, ncol = 5, 
          legend = "none", align = "h", labels =  c("A","B","C","D","E") ,font.label = c(size = 25))

ggsave(combplot, filename = "5_scenarios.png", dpi = 300, type = "cairo", width = 18, height = 10, units = "in")

# Sensitivity Analysis ----------------------------------------------------

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,365,by = 1)
parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 0,
          tstart = 71,
          t_dur = 12*7,
          R0Dec = 0.8)

sens <- list("tstart" = 1:150,
             "R0Dec" = seq(0.1,1.7, by = 0.01),
             "t_dur" = seq(1,365, by = 1))

senslist <- list()

for(i in 1:3) {
  senslist[[i]] <- local({
    sensitivity <- sens[[i]]
    datasens <- data.frame(matrix(nrow = 0, ncol = 4))
    
    for(k in 1:5) {
      data <- cbind(data.frame(matrix(ncol = 4, nrow = length(sensitivity))), "group" = as.factor(seq(1,5)[k]))
      parms["scen"] <- seq(1,5)[k]
      
      if(parms["scen"] != 0 && parms["scen"] != 1 && i != 3) {
        parms["t_dur"] = 24*7
      }
      
      for(j in 1:length(sensitivity)) {
        parms[names(sens)[i]] <- sensitivity[j]
        #print(parms[names(sens)[i]])
        out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
        data[j,1] <- max(out$I)
        data[j,2] <- max(out$C)
        data[j,3] <- sensitivity[j]
        data[j,4] <- names(sens)[i]
      }
      datasens <- rbind(datasens, data)
      print(k)
    }
    
    colnames(datasens) <- c("peak", "cum", "sensval", "sens","group")
    
    p1 <- ggplot(datasens, aes(x = sensval, y = peak, col = group)) + geom_line(size = 1.02) + theme_bw()  + 
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + labs(x = as.character(names(sens)[i]), y = "Peak I(t)", col = "Scenario")
    
    p2 <- ggplot(datasens, aes(x = sensval, y = cum, col = group)) + geom_line(size = 1.02) + theme_bw() +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + labs(x = names(sens)[i], y = "Cumulative Incidence", col = "Scenario")
    
    if(names(sens)[i] == "tstart") {
      p1 <- p1 + labs(x = "Intervention Trigger Date (Days)", y = "Peak I(t)", col = "Scenario") +
        scale_y_continuous(limits = c(0.02,0.125), expand = c(0,0)) + scale_x_continuous(limits = c(0,150),expand = c(0, 0)) 
      p2 <- p2 + labs(x = "Intervention Trigger Date (Days)", y = "Cumulative Incidence", col = "Scenario") +
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) + scale_x_continuous(limits = c(0,150), expand = c(0, 0)) 
    }
    
    if(names(sens)[i] == "R0Dec") {
      p1 <- p1 + labs(x = expression(Intervention ~ R[0]), y = "Peak I(t)", col = "Scenario") +
        scale_y_continuous(limits = c(0.02,0.125),expand = c(0,0)) + scale_x_continuous(limits = c(0,1.7),expand = c(0, 0)) 
      p2 <- p2 + labs(x = expression(Intervention ~ R[0]), y = "Cumulative Incidence", col = "Scenario") +
        scale_y_continuous(limits = c(0,1),expand = c(0,0)) + scale_x_continuous(limits = c(0,1.7),expand = c(0, 0)) 
    }
    
    if(names(sens)[i] == "t_dur") {
      p1 <- p1 + labs(x = "Length of Intervention (Days)", y = "Peak I(t)", col = "Scenario") +
        scale_y_continuous(limits = c(0.02,0.125), expand = c(0,0)) + scale_x_continuous(limits = c(0,365),expand = c(0, 0)) 
      p2 <- p2 + labs(x = "Length of Intervention (Days)", y = "Cumulative Incidence", col = "Scenario") +
        scale_y_continuous(limits = c(0,1),expand = c(0,0)) + scale_x_continuous(limits = c(0,365),expand = c(0, 0)) 
    }
    
    sensplot <- ggarrange(p1,p2, ncol =  2, nrow = 1, common.legend =  TRUE, legend = "none", 
                          align = "h")
    if(names(sens)[i] == "t_dur") {
      sensplot <- ggarrange(p1,p2, ncol =  2, nrow = 1, common.legend =  TRUE, legend = "bottom", 
                            align = "h")
    }
    
    return(sensplot)
  })
}

combsensplot <- ggarrange(senslist[[1]], senslist[[2]], senslist[[3]],nrow = 3, ncol = 1, common.legend = TRUE,
                          legend = "bottom", align = "v", labels =  c("A","B","C") ,font.label = c(size = 20),
                          heights = c(1,1,1.1))

ggsave(combsensplot, filename = "5_scenarios_sensitivity.png", dpi = 300, type = "cairo", width = 9, height = 11, units = "in")



# Multi-Parameter Optimisation --------------------------------------------

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,365,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 0,
          tstart = 71,
          t_dur = 12*7,
          R0Dec = 0.8)

parameterspace <- expand.grid("trigday" = seq(0,175), "length" = seq(0,200))

#For this use sapply - will need to create a function first 

