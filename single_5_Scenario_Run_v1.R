library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Need to rescale the model for c(t) and not beta(t) 

#FUnction to Obtain Betas for the Model 
combbeta <- function(scen, time, tstart, t_dur, cmin) {
  gamma <- 1/GenTime(3, 2.8)
  betascale <- (2.8*gamma)*0.7
  if(scen == 0) {
    output <- betascale
  }
  if(scen == 1) {
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     betascale*cmin, betascale)
  }
  if(scen == 2) {
    cminfun <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(cmin, 1), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     cminfun(time)*betascale, 
                     betascale)
  }
  if(scen == 3) {
    cminfun <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(1, cmin), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     cminfun(time)*betascale, 
                     betascale)
  }
  if(scen == 4) {
    cminfun <-approxfun(x=c(tstart, tstart+(t_dur/2), tstart+t_dur), y= c(1, cmin, 1), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     cminfun(time)*betascale, betascale)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart & time <= tstart+(t_dur*0.16667)) | (time >= tstart+(t_dur*0.3333) & time <= tstart+(t_dur*0.5)) |
                       (time >= tstart+(t_dur*0.6667) & time <= tstart+(t_dur*0.83333)), 
                     betascale*cmin, betascale)
  }
  return(output)
}

plot(seq(0,400),combbeta(5, seq(0,400), 12, 24*7, 0.4))

#SEIR set of ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, cmin)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 


# Baseline Model Analysis - 5 Scenarios -----------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

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
                   "group" =  c("baseline", "scenario")[i], 
                   "beta" = combbeta(explor_scen[i], times, parms[["tstart"]], parms[["t_dur"]], parms[["cmin"]]))
      out$re <- out$beta/parms[["gamma"]]*out$S
      data <- rbind(data, out)
    }
 
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    plotbeta <- melt(data, id.vars = c("time", "group"), measure.vars = ("beta"))
    plotre <- melt(data, id.vars = c("time", "group"), measure.vars = ("re"))
    
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.150),expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
            axis.title.y=element_text(size=18), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
      scale_x_continuous( expand = c(0, 0)) + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkred", "darkred"))

    p2 <- ggplot(plotbeta, aes(x = time, y = value, col = group, alpha= group))  + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.3), expand = c(0,0)) +
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
    
    shade <- data.frame(xmin =  parms[["tstart"]], xmax = parms[["tstart"]]+(parms[["t_dur"]]), ymin = 0, ymax = Inf)
    
    if(parms[["scen"]]  == 1) {
     p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                          fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", col = "", title = "Scenario 1")
     p2 <- p2 + labs(x ="", y = expression(beta[(t)]), col = "")
     p3 <- p3 + labs(x ="Time (Days)", y = expression(R[e]), col = "")
    }
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 2")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 3")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 4")
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
      labs(x ="", y = "", col = "",  title = "Scenario 5")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.45, 0.5))
    dump <- list(combplot, plotdata, max(out$C))
    return(dump)
  })
}

combplot <- ggarrange(datalist[[1]][[1]], datalist[[2]][[1]], datalist[[3]][[1]], datalist[[4]][[1]], datalist[[5]][[1]],nrow = 1, ncol = 5, 
          legend = "none", align = "h")

ggsave(combplot, filename = "5_scenarios.png", dpi = 300, type = "cairo", width = 18, height = 10, units = "in")

for(i in 1:5) {
print(datalist[[i]][[2]][datalist[[i]][[2]]$group == "scenario",]
      [which.max(datalist[[i]][[2]]$value[datalist[[i]][[2]]$group == "scenario"]),])
print(datalist[[i]][[3]])
}


# Sensitivity Analysis ----------------------------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

sens <- list("tstart" = seq(1,150, by = 1),
             "cmin" = seq(0.01,1, by = 0.01),
             "t_dur" = seq(1,400, by = 1))

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
      print(paste0("Sensitivity Test: ", i,", Scenario:", k))
    }
    
    colnames(datasens) <- c("peak", "cum", "sensval", "sens","group")
    
    p1 <- ggplot(datasens, aes(x = sensval, y = peak, col = group)) + geom_line(size = 1.02) + theme_bw()  + 
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.8,0.6,0.4,0.6),"cm")) + labs(x = as.character(names(sens)[i]), y = "Peak I(t)", col = "Scenario")
    
    p2 <- ggplot(datasens, aes(x = sensval, y = cum, col = group)) + geom_line(size = 1.02) + theme_bw() +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.8,0.6,0.4,0.6),"cm")) + labs(x = names(sens)[i], y = "Cumulative Incidence", col = "Scenario")
    
    if(names(sens)[i] == "tstart") {
      p1 <- p1 + labs(x = "Intervention Trigger (Days)", y = "I(t) Peak", col = "Scenario") +
        scale_y_continuous(limits = c(0.02,0.150), expand = c(0,0)) + scale_x_continuous(limits = c(0,150),expand = c(0, 0)) 
      p2 <- p2 + labs(x = "Intervention Trigger (Days)", y = "Cumulative Incidence", col = "Scenario") +
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) + scale_x_continuous(limits = c(0,150), expand = c(0, 0)) 
    }
    
    if(names(sens)[i] == "cmin") {
      p1 <- p1 + labs(x = expression(Intervention ~ c[min]), y = "I(t) Peak", col = "Scenario") +
        scale_y_continuous(limits = c(0.02,0.150),expand = c(0,0)) + scale_x_continuous(limits = c(0,1),expand = c(0, 0)) 
      p2 <- p2 + labs(x = expression(Intervention ~ c[min]), y = "Cumulative Incidence", col = "Scenario") +
        scale_y_continuous(limits = c(0,1),expand = c(0,0)) + scale_x_continuous(limits = c(0,1),expand = c(0, 0)) 
    }
    
    if(names(sens)[i] == "t_dur") {
      p1 <- p1 + labs(x = "Intervention Duration (Days)", y = "I(t) Peak", col = "Scenario") +
        scale_y_continuous(limits = c(0.001,0.150), expand = c(0,0)) + scale_x_continuous(limits = c(0,400),expand = c(0, 0)) 
      p2 <- p2 + labs(x = "Intervention Duration (Days)", y = "Cumulative Incidence", col = "Scenario") +
        scale_y_continuous(limits = c(0,1),expand = c(0,0)) + scale_x_continuous(limits = c(0,400),expand = c(0, 0)) 
    }
    
    sensplot <- ggarrange(p1,p2, ncol =  2, nrow = 1, common.legend =  TRUE, legend = "none", 
                          align = "h")
    if(names(sens)[i] == "t_dur") {
      sensplot <- ggarrange(p1,p2, ncol =  2, nrow = 1, common.legend =  TRUE, legend = "bottom", 
                            align = "h")
    }
    dump <- list(sensplot, datasens)
    return(dump)
  })
}

combsensplot <- ggarrange(senslist[[1]][[1]], senslist[[2]][[1]], senslist[[3]][[1]],nrow = 3, ncol = 1, common.legend = TRUE,
                          legend = "bottom", align = "v", labels =  c("A","B","C") ,font.label = c(size = 25),
                          heights = c(1,1,1.1), vjust = 0.9, hjust = -0.1)

ggsave(combsensplot, filename = "5_scenarios_sensitivity.png", dpi = 300, type = "cairo", width = 9, height = 11, units = "in")

# Multi-Parameter Optimisation --------------------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)

parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

#parameterspace <- expand.grid("trigday" = seq(0,175, by =2), "length" = seq(0,200, by =2), "scen" = seq(1,5))

parameterspace <- expand.grid("trigday" = seq(0,100, by =5), "length" = seq(1,250, by =5))

scensens <- list()

for(j in 1:5) {
  
  scensens[[j]] = local({
    i = 0
    scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 5))
    
    for(i in 1:nrow(parameterspace)) {
      
      print(paste0("Scenario ", j," - ", round(i/nrow(parameterspace), digits = 2)))
      parms["tstart"] <- parameterspace[i,1]
      parms["t_dur"] <- parameterspace[i,2] 
      parms["scen"] <- j
      
      out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
      scendata[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                        "tstart" = parms[["tstart"]], "t_dur" = parms[["t_dur"]])
    }
    
    colnames(scendata) <- c("peak", "cum", "scen", "tstart", "t_dur")
    
    p1 <- ggplot(scendata, aes(x = tstart, y = t_dur, fill= peak))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = -0.2, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + 
      labs(x = "Intervention Trigger", y = "Intervention Duration", fill = "I(t) Peak", title = paste("Scenario", j)) + 
      scale_fill_viridis_c(direction = -1)
    
    p2<- ggplot(scendata, aes(x = tstart, y = t_dur, fill = cum))  + geom_tile() +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = -0.2),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + 
      labs(x = "Intervention Trigger", y = "Intervention Duration", fill = "Cumulative\nIncidence", title = "") + 
      scale_fill_viridis_c(direction = -1, option = "magma") 
    
    combplot <- ggarrange(p1,p2, ncol = 2, nrow = 1, widths = c(1,1.05), align = "h")
    print(combplot)
    dump <- list(combplot, scendata)
    return(dump)
  })
  
}

combplot <- ggarrange(scensens[[1]][[1]],scensens[[2]][[1]],scensens[[3]][[1]],scensens[[4]][[1]],scensens[[5]][[1]],
                      nrow = 5, ncol = 1)

ggsave(combplot, filename = "Heat_5_scenarios_sensitivity.png", dpi = 300, type = "cairo", width = 10, height = 16, units = "in")
