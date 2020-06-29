library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")
start_time <- Sys.time()

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

combbetamult <- function(scen, time, tstart1, t_dur1, tstart2, t_dur2, R0Dec1, R0Dec2) {
  gamma <- 1/GenTime(3.3, 2.8)
  if(scen == 0) {
    output <-  1.7*gamma
  }
  if(scen == 1) {
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse((time >= (tstart1 + t_dur1 + tstart2)),
                            R0Dec2*gamma,
                            R0Dec1*gamma),
                     1.7*gamma)
  }
  if(scen == 2) {
    R0lin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(R0Dec1, 1.7), method="linear", rule =2)
    R0lin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(R0Dec2, 1.7), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0lin1(time)*gamma,
                            R0lin2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 3) {
    R0lin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(1.7, R0Dec1), method="linear", rule =2)
    R0lin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(1.7, R0Dec2), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0lin1(time)*gamma,
                            R0lin2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 4) {
    R0incdec1 <-approxfun(x=c(tstart1, tstart1+(t_dur1/2), tstart1+t_dur1), y= c(1.7, R0Dec1, 1.7), method="linear", rule =2)
    R0incdec2 <-approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 +(t_dur2/2)), (tstart1 + t_dur1 + tstart2 +t_dur2)), y= c(1.7, R0Dec2, 1.7), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0incdec1(time)*gamma,
                            R0incdec2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart1 & time <= tstart1+(t_dur1*0.16667)) | (time >= tstart1+(t_dur1*0.3333) & time <= tstart1+(t_dur1*0.5)) |
                       (time >= tstart1+(t_dur1*0.6667) & time <= tstart1+(t_dur1*0.83333)) |
                       (time >= tstart1 + t_dur1 + tstart2 & time <= tstart1 + t_dur1 + tstart2+(t_dur2*0.16667)) | 
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.3333) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.5)) |
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.6667) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.83333)),
                     ifelse((time >= tstart1 + t_dur1 + tstart2),
                            R0Dec2*gamma,
                            R0Dec1*gamma),
                     1.7*gamma)
  }
  return(output)
}

plot(seq(0,365),combbetamult(5, seq(0,365), 71, 12*7, 20, 12*7, 0.8, 0.6))

#SEIR set of ODEs
SIRmulti <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbetamult(scen, time, tstart1, t_dur1, tstart2, t_dur2, R0Dec1, R0Dec2)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

# Multiple Intervention Case Study - Trajectory ----------------------------------------

# Run the Model


init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,1000,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 0,
          tstart1 = 71,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          R0Dec1 = 0.8,
          R0Dec2 = 0.8)

datalist <- list()

for(j in 1:length(seq(1,5))) {
  datalist[[j]] <- local({
    j=j
    data <- data.frame(matrix(nrow = 9, ncol = 0))
    explor_scen <- c(0, seq(1,5)[j])
    
    
    for(i in 1:2) {
      
      parms["scen"] <- explor_scen[i]
      
      if(parms["scen"] != 0 && parms["scen"] != 1) {
        parms["t_dur1"] = 12*7
        parms["t_dur2"] = 12*7
      }
      
      out <- cbind(data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], 
                   "r0" = combbetamult(explor_scen[i], times, parms[["tstart1"]], parms[["t_dur1"]], parms[["tstart2"]], parms[["t_dur2"]],
                                       parms[["R0Dec1"]], parms[["R0Dec2"]])/parms[["gamma"]])
      out$re <- out$r0*out$S
      data <- rbind(data, out)
    }
    
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    plotr0 <- melt(data, id.vars = c("time", "group"), measure.vars = ("r0"))
    plotre <- melt(data, id.vars = c("time", "group"), measure.vars = ("re"))
    
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.125),expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
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
    
    shade <- data.frame(xmin =  c(parms[["tstart1"]], parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]]), 
                        xmax = c(parms[["tstart1"]]+parms[["t_dur1"]], parms[["tstart1"]]+parms[["t_dur1"]]+parms[["tstart2"]]+parms[["t_dur2"]]), 
                        ymin = 0, ymax = Inf)
    
    if(parms[["scen"]]  == 1) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", col = "", title = "Scenario 1")
      p2 <- p2 + labs(x ="", y = expression(R[0]), col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = expression(R[e]), col = "")
    }
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 2")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 3")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 4")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 5) {
      
      shade <- data.frame(xmin =  c(parms[["tstart1"]], 
                                    parms[["tstart1"]]+(parms[["t_dur1"]]*0.333), 
                                    parms[["tstart1"]]+(parms[["t_dur1"]]*0.667),
                                    parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]],
                                    parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]] + (parms[["t_dur2"]]*0.3333),
                                    parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]] + (parms[["t_dur2"]]*0.667)),
                          
                          
                          xmax = c(parms[["tstart1"]]+(parms[["t_dur1"]]*0.1667), 
                                   parms[["tstart1"]]+(parms[["t_dur1"]]*0.53), 
                                   parms[["tstart1"]]+(parms[["t_dur1"]]*0.833),
                                   parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]] + (parms[["t_dur2"]]*0.16667),
                                   parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]] + (parms[["t_dur2"]]*0.53),
                                   parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]] + (parms[["t_dur2"]]*0.833)),
                          
                          ymin = 0, ymax = Inf)
      
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin,  ymax = ymax, 
                                                              xmin = xmin, xmax = xmax), alpha = 0.2, fill = "darkblue") + geom_line(size = 1.1, stat = "identity") +
        labs(x ="", y = "", col = "",  title = "Scenario 5")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.45, 0.5))
    return(combplot)
  })
}

combplotmulti <- ggarrange(datalist[[1]], datalist[[2]], datalist[[3]], datalist[[4]], datalist[[5]],nrow = 1, ncol = 5, 
                           legend = "none", align = "h")


ggsave(combplotmulti, filename = "5_scenarios_multi.png", dpi = 300, type = "cairo", width = 18, height = 10, units = "in")

# Multiple Optimisations - Trigger Date --------------------------------------------------

optimdata <- expand.grid("tstart1" = seq(0,125, by = 5), "tstart2" = seq(0,125, by = 5))

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,1000,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 1,
          tstart1 = 71,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          R0Dec1 = 0.8,
          R0Dec2 = 0.8)

outcomelist <- list()

for(j in 1:5) {
  j = j
  
  parms["scen"] <- j
  
  if(parms["scen"] != 0 && parms["scen"] != 1) {
    parms["t_dur1"] = 12*7
    parms["t_dur2"] = 12*7
  }
  
  outcomelist[[j]] <- local({
    
    optim <- data.frame(matrix(ncol = 6, nrow = nrow(optimdata)))

    for(i in 1:nrow(optimdata)) {
      parms["tstart1"] <- optimdata[i,1]
      parms["tstart2"] <- optimdata[i,2]
 
      out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
      optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                     "tstart1" = parms[["tstart1"]], "tstart2" = parms[["tstart2"]],
                     "realstart2" = parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]])
    }
    
    colnames(optim) <- c("peak", "cum", "scen", "tstart1", "tstart2", "realstart2")
    
    p1 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= peak))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = -0.2, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + labs(x = "Intervention 1 Trigger", y = "Intervention 2 Trigger", fill = "Peak I(t)", 
                                                        title = paste("Scenario", j)) + 
      scale_fill_viridis_c(direction = -1)
    
    p2 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= cum))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + labs(x = "Intervention 1 Trigger", y = "Intervention 2 Trigger", fill = "Cumulative\nIncidence", title = "") + 
      scale_fill_viridis_c(direction = -1, option = "magma")
    
    combplot <- ggarrange(p1,p2, ncol = 2, nrow = 1, widths = c(1,1.05), align = "h")
    
    print(paste0("Scenario: ", j, " Complete"))
    
   return(combplot)
  })   
}

multicombplot <- ggarrange(outcomelist[[1]],outcomelist[[2]],outcomelist[[3]],outcomelist[[4]],outcomelist[[5]],
                      nrow = 5, ncol = 1)

ggsave(multicombplot, filename = "Heat_5_multi_sensitivity_trig.png", dpi = 300, type = "cairo", width = 10, height = 16, units = "in")

# Multiple Optimisations - R0 --------------------------------------------------

r0optim <- expand.grid("r01" = seq(0, 1.7, by = 0.1), "r01" = seq(0, 1.7, by = 0.1))

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,1000,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 1,
          tstart1 = 71,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          R0Dec1 = 0.8,
          R0Dec2 = 0.8)

outcomelistr0 <- list()

for(j in 1:5) {
  j = j
  
  parms["scen"] <- j
  
  if(parms["scen"] != 0 && parms["scen"] != 1) {
    parms["t_dur1"] = 12*7
    parms["t_dur2"] = 12*7
  }
  
  outcomelistr0[[j]] <- local({
    
    optim <- data.frame(matrix(ncol = 5, nrow = nrow(r0optim)))
    
    for(i in 1:nrow(r0optim)) {
      parms["R0Dec1"] <- r0optim[i,1]
      parms["R0Dec2"] <- r0optim[i,2]
      
      out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
      optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                     "R0Dec1" = parms[["R0Dec1"]], "R0Dec2" = parms[["R0Dec2"]])
    }
    
    colnames(optim) <- c("peak", "cum", "scen", "R0Dec1", "R0Dec2")
    
    p1 <- ggplot(optim, aes(x = R0Dec1, y = R0Dec2, fill= peak))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = -0.2, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + labs(x = "Intervention 1 R0", y = "Intervention 2 R0", fill = "Peak I(t)", 
                                                        title = paste("Scenario", j)) + 
      scale_fill_viridis_c(direction = -1)
    
    p2 <- ggplot(optim, aes(x = R0Dec1, y = R0Dec2, fill= cum))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + 
      scale_fill_viridis_c(direction = -1, option = "magma") + labs(x = "Intervention 1 R0", y = "Intervention 2 R0", fill = "Cumulative\nIncidence", 
                                                                    title = "")
    
    combplot <- ggarrange(p1,p2, ncol = 2, nrow = 1, widths = c(1,1.05), align = "h")
    
    print(paste0("Scenario: ", j, " Complete"))
    
    return(combplot)
  })   
}

multicombplotr0 <- ggarrange(outcomelistr0[[1]],outcomelistr0[[2]],outcomelistr0[[3]],outcomelistr0[[4]],outcomelistr0[[5]],
                           nrow = 5, ncol = 1)

ggsave(multicombplotr0, filename = "Heat_5_multi_sensitivity_r0.png", dpi = 300, type = "cairo", width = 10, height = 16, units = "in")


end_time <- Sys.time()
end_time - start_time
