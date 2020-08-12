library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis/supplementary")
start_time <- Sys.time()

# Model Functions ----------------------------------------------------------
#Generation Time Function
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Beta Function for Single Intervention
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

#Multiple Interventions
combbetamult <- function(scen, time, tstart1, t_dur1, tstart2, t_dur2, cmin1, cmin2) {
  gamma <- 1/GenTime(3, 2.8)
  betascale <- (2.8*gamma)*0.7
  if(scen == 0) {
    output <-  betascale
  }
  if(scen == 1) {
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse((time >= (tstart1 + t_dur1 + tstart2)),
                            betascale*cmin2,
                            betascale*cmin1),
                     betascale)
  }
  if(scen == 2) {
    betalin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(cmin1, 1), method="linear", rule =2)
    betalin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(cmin2, 1), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            betalin1(time)*betascale,
                            betalin2(time)*betascale),
                     betascale)
  }
  if(scen == 3) {
    betalin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(1, cmin1), method="linear", rule =2)
    betalin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(1, cmin2), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            betalin1(time)*betascale,
                            betalin2(time)*betascale),
                     betascale)
  }
  if(scen == 4) {
    betaincdec1 <-approxfun(x=c(tstart1, tstart1+(t_dur1/2), tstart1+t_dur1), y= c(1, cmin1, 1), method="linear", rule =2)
    betaincdec2 <-approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 +(t_dur2/2)), (tstart1 + t_dur1 + tstart2 +t_dur2)), y= c(1, cmin2, 1), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            betaincdec1(time)*betascale,
                            betaincdec2(time)*betascale),
                     betascale)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart1 & time <= tstart1+(t_dur1*0.16667)) | (time >= tstart1+(t_dur1*0.3333) & time <= tstart1+(t_dur1*0.5)) |
                       (time >= tstart1+(t_dur1*0.6667) & time <= tstart1+(t_dur1*0.83333)) |
                       (time >= tstart1 + t_dur1 + tstart2 & time <= tstart1 + t_dur1 + tstart2+(t_dur2*0.16667)) | 
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.3333) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.5)) |
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.6667) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.83333)),
                     ifelse((time >= tstart1 + t_dur1 + tstart2),
                            cmin2*betascale,
                            cmin1*betascale),
                     betascale)
  }
  return(output)
}

plot(seq(0,365),combbetamult(3, seq(0,365), 52, 6*7, 20, 6*7, 0.4, 1))

#SEIR set of ODEs - Single
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

#SEIR set of ODEs - Repeated
SIRmulti <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbetamult(scen, time, tstart1, t_dur1, tstart2, t_dur2, cmin1, cmin2)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

# Single Intervention Sensitivity Analysis --------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,365,by = 1)

parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

parameterspace <- expand.grid("trigday" = seq(0,100, by =5), "length" = seq(1,200, by =5))
cminrange <- c(0.25, 0.5, 0.75)

scensens <- list()

for(j in 1:5) {
  
  scensens[[j]] = local({
    i = 0
    scendataframe <- data.frame(matrix(nrow = 0, ncol = 6))
    parms["scen"] <- j
    
    for(z in 1:length(cminrange)) {
      scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 6))
      
      parms["cmin"] <- cminrange[z]
      
      for(i in 1:nrow(parameterspace)) {
        parms["tstart"] <- parameterspace[i,1]
        parms["t_dur"] <- parameterspace[i,2] 
        
        out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
        scendata[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                          "tstart" = parms[["tstart"]], "t_dur" = parms[["t_dur"]], "cmin" = parms[["cmin"]])
      }
      
      print(paste0("Scenario ", j," | cmin: ", cminrange[z]))
      scendataframe <- rbind(scendataframe, scendata)
    }
    
    colnames(scendataframe) <- c("peak", "cum", "scen", "tstart", "t_dur","cmin")
    
    datalist <- list("p1data" = scendataframe[scendataframe$cmin == 0.25,],
                     "p2data" = scendataframe[scendataframe$cmin == 0.5,],
                     "p3data" = scendataframe[scendataframe$cmin == 0.75,])
    
    for(t in 1:3) { 
      data <- datalist[[t]]
      
      p1 <- ggplot(data, aes(x = tstart, y = t_dur, fill= peak))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15), axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = - 0.2, face = "bold"),
              plot.subtitle = element_text(size = 15, vjust = 2, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm"))
      
      p2 <- ggplot(data, aes(x = tstart, y = t_dur, fill = cum))  + geom_tile() +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
              plot.subtitle = element_text(size = 15, vjust = 2, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) 
      
      #COmmon Legends across Cmins
      if(j == 1) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.01, 0.15, by = (0.15-0.01)/4), limits = c(0.01, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks= seq(0, 0.8, by = (0.8)/4), limits = c(0, 0.8) )  
      }
      if(j == 2) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.025, 0.15, by = (0.15-0.025)/4), limits = c(0.025, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      if(j == 3) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.03, 0.15, by = (0.15-0.03)/4), limits = c(0.03, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      if(j == 4) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.03, 0.15, by = (0.15-0.03)/4), limits = c(0.03, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      if(j == 5) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.04, 0.15, by = (0.15-0.04)/4), limits = c(0.04, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      
      #Differentiating cmin titles
      if(t == 1) {
        p1 <- p1 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = paste0("Scenario ", j),
                        subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.25))
        p2 <- p2 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "Total\nCumulative\nIncidence", title = paste0("Scenario ", j),
                        subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.25))
      }
      if(t == 2) {
        p1 <- p1 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = "",
                        subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.5))
        p2 <- p2 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "Total\nCumulative\nIncidence", title = "",
                        subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.5))
      }
      if(t == 3) {
        p1 <- p1 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = "",
                        subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.75))
        p2 <- p2 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "Total\nCumulative\nIncidence", title = "",
                        subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.75))
      }
      
      assign(paste0("p",1,"peak",c(0, 5, 1)[t]), p1)
      assign(paste0("p",1,"cum",c(0, 5, 1)[t]), p2)
      
    }
    
    combplotpeak <- ggarrange(p1peak0, p1peak5, p1peak1, ncol = 3, nrow = 1, widths = c(1,1,1), align = "h")
    combplotcum <- ggarrange(p1cum0, p1cum5, p1cum1, ncol = 3, nrow = 1, widths = c(1,1,1), align = "h")
    
    return(list(combplotpeak,combplotcum))
  })
  
}

combplotsenspeak <- ggarrange(scensens[[1]][[1]],scensens[[2]][[1]],scensens[[3]][[1]],scensens[[4]][[1]],scensens[[5]][[1]],
                              nrow = 5, ncol = 1)
combplotsenscum <- ggarrange(scensens[[1]][[2]],scensens[[2]][[2]],scensens[[3]][[2]],scensens[[4]][[2]],scensens[[5]][[2]],
                             nrow = 5, ncol = 1)

ggsave(combplotsenspeak, filename = "cmin_Heat_5_sensitivity_peak.png", dpi = 300, type = "cairo", width = 13, height = 16, units = "in")
ggsave(combplotsenscum, filename = "cmin_Heat_5_sensitivity_cum.png", dpi = 300, type = "cairo", width = 13, height = 16, units = "in")

# Multiple Optimisations - Trigger + Length --------------------------------------------------

lengthdata <- expand.grid("int1length" = seq(3,9, by = 3), "int2length" = seq(3,9, by = 3))
optimdata <- expand.grid("tstart1" = seq(0,100, by = 5), "tstart2" = seq(0,100, by = 5))
lengthdata2345 <- expand.grid("int1length" = seq(6,18, by = 6), "int2length" = seq(6,18, by = 6))

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,420,by = 1)

parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 2,
          tstart1 = 52,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          cmin1 = 0.4,
          cmin2 = 0.4)

outcomelist <- list()

for (z in 1:5) {
  
  parms["scen"] = z
  
  
  for(j in 1:nrow(lengthdata)) {
    j = j
    
    parms["t_dur1"] <- lengthdata[j,1] * 7
    parms["t_dur2"] <- lengthdata[j,2] * 7 
    
    if(parms["scen"] != 0 && parms["scen"] != 1) {
      parms["t_dur1"] = lengthdata2345[j,1] * 7
      parms["t_dur2"] = lengthdata2345[j,2] * 7
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
      bquote("Trigger ("*italic(t[p])*")")
      colnames(optim) <- c("peak", "cum", "scen", "tstart1", "tstart2", "realstart2")
      
      p1 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= peak))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + labs(x = bquote("Trigger 1 ("*italic(t[p1])*")"), y = bquote("Trigger 2 ("*italic(t[p2])*")"), fill = "I(t) Peak", 
                                                          title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " wks"))
      
      p2 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= cum))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + labs(x = bquote("Trigger 1 ("*italic(t[p1])*")"), y = bquote("Trigger 2 ("*italic(t[p2])*")"), fill = "Total\nCumulative\nIncidence", 
                                                          title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " wks"))
      
      if(z == 1) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.044, 0.15, by = (0.15-0.044)/4), limits = c(0.044, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks= seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      if(z == 2) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.042, 0.15, by = (0.15-0.042)/4), limits = c(0.042, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.45, 0.8, by = (0.8-0.45)/4), limits = c(0.45, 0.8) )  
      }
      if(z == 3) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.035, 0.15, by = (0.15-0.035)/4), limits = c(0.035, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      if(z == 4) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.042, 0.15, by = (0.15-0.042)/4), limits = c(0.042, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      if(z == 5) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.05, 0.15, by = (0.15-0.05)/4), limits = c(0.05, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      
      print(paste0("Scenario: ", z, " | Length Analysis: ", (j/nrow(lengthdata))*100, "%"))
      
      return(list(p1,p2))
    })   
  }
  
  peak_multilength <- ggarrange(outcomelist[[1]][[1]], outcomelist[[2]][[1]], outcomelist[[3]][[1]],
                                outcomelist[[4]][[1]], outcomelist[[5]][[1]], outcomelist[[6]][[1]],
                                outcomelist[[7]][[1]], outcomelist[[8]][[1]], outcomelist[[9]][[1]],
                                nrow = 3, ncol =3, align = "hv")
  
  ggsave(peak_multilength, filename = paste0("heatpeak_multilength",z,".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
  
  cum_multilength <- ggarrange(outcomelist[[1]][[2]], outcomelist[[2]][[2]], outcomelist[[3]][[2]],
                               outcomelist[[4]][[2]], outcomelist[[5]][[2]], outcomelist[[6]][[2]],
                               outcomelist[[7]][[2]], outcomelist[[8]][[2]], outcomelist[[9]][[2]],
                               nrow = 3, ncol =3, align = "hv")
  
  ggsave(cum_multilength, filename = paste0("heatcum_multilength", z, ".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
}

# Multiple Optimisations - CMIN + Length --------------------------------------------------

lengthdata <- expand.grid("int1length" = seq(3,9, by = 3), "int2length" = seq(3,9, by = 3))
lengthdata2345 <- expand.grid("int1length" = seq(6,18, by = 6), "int2length" = seq(6,18, by = 6))
cminoptim <- expand.grid("cmin1" = seq(0, 1, by = 0.05), "cmin2" = seq(0, 1, by = 0.05))

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,1000,by = 1)

parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 2,
          tstart1 = 52,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          cmin1 = 0.4,
          cmin2 = 0.4)

outcomelistcmin <- list()

for (z in 1:5) {
  parms["scen"] = z
  
  
  for(j in 1:nrow(lengthdata)) {
    j = j
    
    parms["t_dur1"] <- lengthdata[j,1] * 7
    parms["t_dur2"] <- lengthdata[j,2] * 7 
    
    if(parms["scen"] != 0 && parms["scen"] != 1) {
      parms["t_dur1"] = lengthdata2345[j,1] * 7
      parms["t_dur2"] = lengthdata2345[j,2] * 7
    }
    
    outcomelistcmin[[j]] <- local({
      
      optim <- data.frame(matrix(ncol = 5, nrow = nrow(cminoptim)))
      
      for(i in 1:nrow(cminoptim)) {
        parms["cmin1"] <- cminoptim[i,1]
        parms["cmin2"] <- cminoptim[i,2]
        
        out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
        optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                       "cmin1" = parms[["cmin1"]], "cmin2" = parms[["cmin2"]])
      }
      
      colnames(optim) <- c("peak", "cum", "scen", "cmin1", "cmin2")
      
      p1 <- ggplot(optim, aes(x = cmin1, y = cmin2, fill= peak))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + labs(x = bquote(.(Intervention ~ 1 ~ italic(c[min]))), y = bquote(.(Intervention ~ 2 ~ italic(c[min]))), fill = "I(t) Peak", 
                                                          title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " weeks"))
      
      p2 <- ggplot(optim, aes(x = cmin1, y = cmin2, fill= cum))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + 
        labs(x = bquote(.(Intervention ~ 1 ~ italic(c[min]))), y = bquote(.(Intervention ~ 2 ~ italic(c[min]))), fill = "Total\nCumulative\nIncidence", 
                                                                      title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " wks"))
      
      
      if(z == 1) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.04, 0.15, by = (0.15-0.04)/4), limits = c(0.04, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks= seq(0.46, 0.8, by = (0.8-0.46)/4), limits = c(0.46, 0.8) )  
      }
      if(z == 2) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.05, 0.15, by = (0.15-0.05)/4), limits = c(0.05, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.48, 0.8, by = (0.8-0.48)/4), limits = c(0.48, 0.8) )  
      }
      if(z == 3) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.04, 0.15, by = (0.15-0.04)/4), limits = c(0.04, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.43, 0.8, by = (0.8-0.43)/4), limits = c(0.43, 0.8) )  
      }
      if(z == 4) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.03, 0.15, by = (0.15-0.03)/4), limits = c(0.03, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.45, 0.8, by = (0.8-0.45)/4), limits = c(0.45, 0.8) )  
      }
      if(z == 5) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.05, 0.15, by = (0.15-0.05)/4), limits = c(0.05, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.48, 0.8, by = (0.8-0.48)/4), limits = c(0.48, 0.8) )  
      }
      
      print(paste0("Scenario: ", z, " | Length Analysis: "))
      
      return(list(p1,p2))
    })   
  }
  
  peak_multilength_cmin <- ggarrange(outcomelistcmin[[1]][[1]], outcomelistcmin[[2]][[1]], outcomelistcmin[[3]][[1]],
                                     outcomelistcmin[[4]][[1]], outcomelistcmin[[5]][[1]], outcomelistcmin[[6]][[1]],
                                     outcomelistcmin[[7]][[1]], outcomelistcmin[[8]][[1]], outcomelistcmin[[9]][[1]],
                                nrow = 3, ncol =3, align = "hv")
  
  ggsave(peak_multilength_cmin, filename = paste0("heatpeakcmin_multilength",z,".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
  
  cum_multilength_cmin <- ggarrange(outcomelistcmin[[1]][[2]], outcomelistcmin[[2]][[2]], outcomelistcmin[[3]][[2]],
                                    outcomelistcmin[[4]][[2]], outcomelistcmin[[5]][[2]], outcomelistcmin[[6]][[2]],
                                    outcomelistcmin[[7]][[2]], outcomelistcmin[[8]][[2]], outcomelistcmin[[9]][[2]],
                               nrow = 3, ncol =3, align = "hv")
  
  ggsave(cum_multilength_cmin, filename = paste0("heatcumcmin_multilength", z, ".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
}

# Multiple Intervention Case Study - cmin Trajectory - All Scenarios ----------------------------------------

cminoptim <- expand.grid("cmin1" = c(0.25, 0.5, 0.75), "cmin2" = c(0.25, 0.5, 0.75))

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)

parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 1,
          tstart1 = 52,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          cmin1 = 0.4,
          cmin2 = 0.4)

datalist <- list()

for(z in 1:5){
  
  for(j in 1:nrow(cminoptim)) {
    
    datalist[[j]] <- local({
      
      j=j
      
      parms["cmin1"] = cminoptim[j,1]
      parms["cmin2"] = cminoptim[j,2]
      parms["scen"] = z
      
      if(parms["scen"] != 0 && parms["scen"] != 1) {
        parms["t_dur1"] = 12*7
        parms["t_dur2"] = 12*7
      }

      out <- cbind(data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms)), 
                   "beta" = combbetamult(z, times, parms[["tstart1"]], parms[["t_dur1"]], parms[["tstart2"]], parms[["t_dur2"]],
                                         parms[["cmin1"]], parms[["cmin2"]]))
      out$re <- (out$beta/parms["gamma"])*out$S
      
      shade <- data.frame(xmin =  c(parms[["tstart1"]], parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]]), 
                          xmax = c(parms[["tstart1"]]+parms[["t_dur1"]], parms[["tstart1"]]+parms[["t_dur1"]]+parms[["tstart2"]]+parms[["t_dur2"]]), 
                          ymin = 0, ymax = Inf)
      
      p1 <- ggplot(data = out, aes(x = time, y = I)) + theme_bw() +
        scale_y_continuous(limits = c(0 , 0.15),expand = c(0,0)) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
              plot.title = element_text(size = 18, vjust = 2, hjust = 0.5, face = "bold"),
              axis.title.y=element_text(size=15), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        scale_x_continuous(expand = c(0, 0))  +
        geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                  fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", 
                                                                                       title = paste0(cminoptim[j,1], " / ", cminoptim[j,2]))
      p2 <- ggplot(out, aes(x = time, y = beta))  + theme_bw() +
        scale_y_continuous(limits = c(0 , 0.4), expand = c(0,0)) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
              axis.title.y=element_text(size=15),axis.title.x = element_blank(), 
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
        geom_line(size = 1.1, stat = "identity")  + labs(x ="", y = bquote(italic(beta["(t)"])))
      
      p3 <- ggplot(out, aes(x = time, y = re)) + theme_bw() +
        scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) +
        theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=12),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
        geom_hline(yintercept = 1, size = 1.1, lty = 2, col = "black") + geom_line(size = 1.02, stat = "identity") + labs(x ="Time (Days)", y = expression(R[e]))
      
      combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.45, 0.55))
      return(combplot)
    })
  }
  
  
  combplotmulti_TRAJ <- ggarrange(datalist[[1]], datalist[[2]], datalist[[3]], 
                                  datalist[[4]], datalist[[5]], datalist[[6]],
                                  datalist[[7]], datalist[[8]], datalist[[9]],
                                  nrow = 3, ncol = 3, legend = "none", align = "h")
  
  ggsave(combplotmulti_TRAJ, filename = paste0("1_scenarios_multi_TRAJ_", z,".png"), dpi = 300, type = "cairo", width = 10, height = 15, units = "in")
  
}
