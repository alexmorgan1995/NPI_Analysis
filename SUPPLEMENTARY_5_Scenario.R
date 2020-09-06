library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis/supplementary")
start_time <- Sys.time()

#Generation Time as a function of doubling time and the R0 
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

# Beta Model Functions ----------------------------------------------------------
#Function to describe beta(t) for single intervention
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

#Function to describe beta(t) for multi-intervention
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


# ODE Functions --------------------------------------------------------------------

#SIR set of ODEs - Single
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

#SIR set of ODEs - Repeated
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

# Plots describing Beta(t) - For Supplementary Table --------------------------------------------------------------

#Initial Conditions
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 3),
          scen = 3,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

betaplots <- list()

#Plotting Loop - will save to working directory
for(i in 1:5) {
  betaplots[[i]] <- local({
    parms["scen"] = seq(1:5)[i]
    out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)),
                 "beta" = combbeta(parms["scen"], times, parms[["tstart"]], parms[["t_dur"]], parms[["cmin"]]))
    out$re <- out$beta/parms[["gamma"]]*out$S    
    
    p1 <- ggplot(data = out, aes(x = time, y = beta)) + theme_bw() + geom_line(size = 1.1, stat = "identity") + 
      scale_y_continuous(limits = c(0 , 0.3), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=18),axis.title.x = element_text(size=15), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = expression(beta[(t)]))
    
    return(p1)
    
  })
  ggsave(betaplots[[i]], filename = paste0("Beta_Scenario_", seq(1:5)[i],".png"), dpi = 300, type = "cairo", width = 6, height = 3, units = "in")
}

# Exploring c_min for the single intervention sensitivity analysis - Figure S3-4 --------------------------------

#Initial Conditions and parameters
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,365,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

#Parameter range explored and the cmin values explored for each analysis 
parameterspaceOG <- expand.grid("trigday" = seq(0,100, by =5), "length" = seq(1,250, by =5))
parameterspacescen1 <- expand.grid("trigday" = seq(0,100, by =5), "length" = seq(1,125, by = 2.5))
cminrange <- c(0.25, 0.5, 0.75)

scensens <- list()

#Plotting loop for each of the 5 scenarios 
for(j in 1:5) {
  
  scensens[[j]] = local({
    i = 0
    
    scendataframe <- data.frame(matrix(nrow = 0, ncol = 6)) #Initialise empty dataframe
    
    #Setting the explored scenario - if scenario 1 then use a limited range for dt
    #This is part of the dt range compensations - dt is then scaled for scenario 1 to provide a relative axis 
    parameterspace <- parameterspaceOG
    parms["scen"] <- j
    if(parms["scen"] == 1) {
      parameterspace <- parameterspacescen1
    }
    
    #Run a plotting loop - for each of the considered C_min values - (length of 3)
    for(z in 1:length(cminrange)) {
      
      scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 6)) #Initialise empty dataframe 
      
      parms["cmin"] <- cminrange[z] #for the explored c_min values - set this in the parameters
      
      for(i in 1:nrow(parameterspace)) { #For the actual sensitivity analysis - iterate through the explored dt and tp 
        parms["tstart"] <- parameterspace[i,1]
        parms["t_dur"] <- parameterspace[i,2] 
        
        #Run the model 
        out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
        scendata[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                          "tstart" = parms[["tstart"]], "t_dur" = parms[["t_dur"]], "cmin" = parms[["cmin"]])
      }
      
      print(paste0("Scenario ", j," | cmin: ", cminrange[z])) #Progress notification
      scendataframe <- rbind(scendataframe, scendata) #Bind to the empty dataframe
    }
    
    colnames(scendataframe) <- c("peak", "cum", "scen", "tstart", "t_dur","cmin")
    
    #create a list of the three dataframes - stratified by the cmin value explored
    datalist <- list("p1data" = scendataframe[scendataframe$cmin == 0.25,],
                     "p2data" = scendataframe[scendataframe$cmin == 0.5,],
                     "p3data" = scendataframe[scendataframe$cmin == 0.75,])
    
    for(t in 1:3) { #Identical to z - iterate through 3 values (the cmin values explored)
      data <- datalist[[t]]
      
      #Generic plotting for the peak prevalence (p1) and attack rate (p2)
      p1 <- ggplot(data, aes(x = tstart, y = t_dur, fill= peak))  + geom_tile()  +
        scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15), axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = - 0.2, face = "bold"),
              plot.subtitle = element_text(size = 15, vjust = 2, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm"))
      
      p2 <- ggplot(data, aes(x = tstart, y = t_dur, fill = cum))  + geom_tile() +
        scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
              plot.subtitle = element_text(size = 15, vjust = 2, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) 
      
      #Create common legends for each scenario for the peak prevalence and the attack rate 
      if(j == 1) {
        
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.035, 0.15, by = (0.15-0.035)/4), limits = c(0.035, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks= seq(0.25, 0.8, by = (0.8-0.25)/4), limits = c(0.25, 0.8) )  
      }
      if(j == 2) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.022, 0.15, by = (0.15-0.022)/4), limits = c(0.022, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.32, 0.8, by = (0.8-0.32)/4), limits = c(0.32, 0.8) )  
      }
      if(j == 3) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.03, 0.15, by = (0.15-0.03)/4), limits = c(0.03, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.36, 0.8, by = (0.8-0.36)/4), limits = c(0.36, 0.8) )  
      }
      if(j == 4) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.028, 0.15, by = (0.15-0.028)/4), limits = c(0.028, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.34, 0.8, by = (0.8-0.34)/4), limits = c(0.34, 0.8) )  
      }
      if(j == 5) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.04, 0.15, by = (0.15-0.04)/4), limits = c(0.04, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.36, 0.8, by = (0.8-0.36)/4), limits = c(0.36, 0.8) )  
      }
      
      #Axis scaling for scenario 1
      formatter2 <- function(x){ 
        x*2
      }
      
      #Formatting according to the cmin value explored - the axis titles and labels have to be different for each plot
      
      if(t == 1) {
        
        if(parms["scen"] == 1){
          p1 <- p1 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Relative Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = paste0("Scenario ", j),
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.25)) + scale_y_continuous(expand = c(0,0), labels = formatter2)
          p2 <- p2 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Relative Duration ("*italic(d[t])*")"), fill = "Attack\nRate", title = paste0("Scenario ", j),
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.25)) + scale_y_continuous(expand = c(0,0), labels = formatter2)
          
        } else { 
          p1 <- p1 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = paste("Scenario", j),
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.25)) + 
          scale_y_continuous(expand = c(0,0))
          p2 <- p2 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "Attack\nRate", title = paste("Scenario", j),
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.25)) + 
            scale_y_continuous(expand = c(0,0))
          }
      }
      if(t == 2) {
        
        if(parms["scen"] == 1){
          p1 <- p1 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Relative Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.5)) + scale_y_continuous(expand = c(0,0), labels = formatter2)
          p2 <- p2 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Relative Duration ("*italic(d[t])*")"), fill = "Attack\nRate", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.5)) + scale_y_continuous(expand = c(0,0), labels = formatter2)
          
        } else { 
          p1 <- p1 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.5)) + 
            scale_y_continuous(expand = c(0,0))
          p2 <- p2 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "Attack\nRate", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.5)) + 
            scale_y_continuous(expand = c(0,0))
        }
      }
      if(t == 3) {
        if(parms["scen"] == 1){
          p1 <- p1 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Relative Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.75)) + scale_y_continuous(expand = c(0,0), labels = formatter2)
          p2 <- p2 + labs(x = bquote("Trigger ("*italic(t[p])*")"), y = bquote("Relative Duration ("*italic(d[t])*")"), fill = "Attack\nRate", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.75)) + scale_y_continuous(expand = c(0,0), labels = formatter2)
          
        } else { 
          p1 <- p1 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.75)) + 
            scale_y_continuous(expand = c(0,0))
          p2 <- p2 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "Attack\nRate", title = "",
                          subtitle = bquote("Intervention" ~ italic(c[min])~"="~0.75)) + 
            scale_y_continuous(expand = c(0,0))
        }
      }
      
      #change the name of plot 1 and 2 according to what cmin value has been explored 
      assign(paste0("p",1,"peak",c(0, 5, 1)[t]), p1)
      assign(paste0("p",1,"cum",c(0, 5, 1)[t]), p2)
      
    }
    
    #Combine the different cmin runs for each scenario together - so that each scenario has 3 panels for each cmin
    #This is done for peak prevalence and the attack rate 
    combplotpeak <- ggarrange(p1peak0, p1peak5, p1peak1, ncol = 3, nrow = 1, widths = c(1,1,1), align = "h")
    combplotcum <- ggarrange(p1cum0, p1cum5, p1cum1, ncol = 3, nrow = 1, widths = c(1,1,1), align = "h")
    
    return(list(combplotpeak,combplotcum))
  })
  
}

#Create 2 figures one for peak prevalence and the attack rate 
combplotsenspeak <- ggarrange(scensens[[1]][[1]],scensens[[2]][[1]],scensens[[3]][[1]],scensens[[4]][[1]],scensens[[5]][[1]],
                              nrow = 5, ncol = 1)
combplotsenscum <- ggarrange(scensens[[1]][[2]],scensens[[2]][[2]],scensens[[3]][[2]],scensens[[4]][[2]],scensens[[5]][[2]],
                             nrow = 5, ncol = 1)

#Save the two plots to the working directory 
ggsave(combplotsenspeak, filename = "cmin_Heat_5_sensitivity_peak.png", dpi = 300, type = "cairo", width = 13, height = 16, units = "in")
ggsave(combplotsenscum, filename = "cmin_Heat_5_sensitivity_cum.png", dpi = 300, type = "cairo", width = 13, height = 16, units = "in")

# Multiple Optimisations - Trigger + Length --------------------------------------------------

#Parameter range for the trigger and length parameter sensitivity analysis
lengthdata <- expand.grid("int1length" = seq(3,9, by = 3), "int2length" = seq(3,9, by = 3)) #exploring dt as a parameter for scenario 1
optimdata <- expand.grid("tstart1" = seq(0,100, by = 5), "tstart2" = seq(0,100, by = 10)) 
lengthdata2345 <- expand.grid("int1length" = seq(6,18, by = 6), "int2length" = seq(6,18, by = 6)) #exploring dt as a parameter for all other scenarios
#Differences in scenario as dt compensations are necessary for scenario 1

#Initial conditions and parameter values 
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

outcomelist <- list() #Initialising empty list for data and plot storage

for (z in 1:5) { #Iterate for every considered scenario
  
  parms["scen"] = z

  for(j in 1:nrow(lengthdata)) {#For the dt values explored (only 3x3)
    j = j
    
    if(parms["scen"] != 0 && parms["scen"] != 1) {
      parms["t_dur1"] = lengthdata2345[j,1] * 7
      parms["t_dur2"] = lengthdata2345[j,2] * 7
    } else {
      parms["t_dur1"] <- lengthdata[j,1] * 7
      parms["t_dur2"] <- lengthdata[j,2] * 7 
    }
    
    outcomelist[[j]] <- local({#For the chosen dt values run the sensitivity analysis
      
      optim <- data.frame(matrix(ncol = 6, nrow = nrow(optimdata))) #Initialise the empty dataframe 
      
      for(i in 1:nrow(optimdata)) { #iterate for all combinations of trigger date 1 and 2 
        parms["tstart1"] <- optimdata[i,1] #Specify both parameters in the model
        parms["tstart2"] <- optimdata[i,2]
        
        #Run the model
        out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
        optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                       "tstart1" = parms[["tstart1"]], "tstart2" = parms[["tstart2"]],
                       "realstart2" = parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]])
      }
      bquote("Trigger ("*italic(t[p])*")")
      colnames(optim) <- c("peak", "cum", "scen", "tstart1", "tstart2", "realstart2")
      
      #Generic plotting for the peak prevalence and the attack rate
      p1 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= peak))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=18),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.8, "cm"),
              legend.key.width =  unit(2, "cm")) + labs(x = bquote("Trigger 1 ("*italic(t[p1])*")"), y = bquote("Trigger 2 ("*italic(t[p2])*")"), fill = "I(t) Peak", 
                                                          title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " wks"))
      
      p2 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= cum))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=18),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.8, "cm"),
              legend.key.width =  unit(2, "cm")) + labs(x = bquote("Trigger 1 ("*italic(t[p1])*")"), y = bquote("Trigger 2 ("*italic(t[p2])*")"), fill = "Attack\nRate", 
                                                          title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " wks"))
      
      #Scenario specific plotting for the peak prevalence and atatck rate - have made the legends for each scenario common
      if(z == 1) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.044, 0.15, by = (0.15-0.044)/4), limits = c(0.044, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks= seq(0.4, 0.8, by = (0.8-0.4)/4), limits = c(0.4, 0.8) )  
      }
      if(z == 2) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.038, 0.15, by = (0.15-0.038)/4), limits = c(0.038, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.45, 0.8, by = (0.8-0.45)/4), limits = c(0.45, 0.8) )  
      }
      if(z == 3) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.03, 0.15, by = (0.15-0.03)/4), limits = c(0.03, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.42, 0.8, by = (0.8-0.42)/4), limits = c(0.42, 0.8) )  
      }
      if(z == 4) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.034, 0.15, by = (0.15-0.034)/4), limits = c(0.034, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.42, 0.8, by = (0.8-0.42)/4), limits = c(0.42, 0.8) )  
      }
      if(z == 5) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.05, 0.15, by = (0.15-0.05)/4), limits = c(0.05, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.44, 0.8, by = (0.8-0.44)/4), limits = c(0.44, 0.8) )  
      }
      
      print(paste0("Scenario: ", z, " | Length Analysis: ", (j/nrow(lengthdata))*100, "%")) #Progress notification for the iterator
      
      return(list(p1,p2))
    })   
  }
  
  #Once we have obtained the sensitivity analysis data for each scenario across all the dt values we concatenate all of them together
  #This is done for the peak prevalence and the attack rate
  peak_multilength <- ggarrange(outcomelist[[1]][[1]], outcomelist[[2]][[1]], outcomelist[[3]][[1]],
                                outcomelist[[4]][[1]], outcomelist[[5]][[1]], outcomelist[[6]][[1]],
                                outcomelist[[7]][[1]], outcomelist[[8]][[1]], outcomelist[[9]][[1]],
                                nrow = 3, ncol =3, align = "hv", common.legend = TRUE, legend = "bottom")
  cum_multilength <- ggarrange(outcomelist[[1]][[2]], outcomelist[[2]][[2]], outcomelist[[3]][[2]],
                               outcomelist[[4]][[2]], outcomelist[[5]][[2]], outcomelist[[6]][[2]],
                               outcomelist[[7]][[2]], outcomelist[[8]][[2]], outcomelist[[9]][[2]],
                               nrow = 3, ncol =3, align = "hv", common.legend = TRUE, legend = "bottom")
  
  #Then combine the two outcome measure plots together 
  combplottstart <- ggarrange(NULL,NULL, 
                              NULL,peak_multilength, 
                              NULL,NULL,
                              NULL,cum_multilength,
                                nrow = 4, ncol = 2, align = "v", heights = c(0.1,1,0.1,1), widths = c(0.05,1),
                              labels = c("","",
                                         "","A",
                                         "","", 
                                         "","B"), font.label = c(size = 45), vjust = -1.1, hjust = 0.8)
  
  #Save the sensitivity analysis according to the scenario explored (z)
  ggsave(combplottstart, filename = paste0("heatcum_multilength_comb",z,".png"), dpi = 300, type = "cairo", width = 15, height = 22, units = "in")
  
}

# Multiple Optimisations - CMIN + Length --------------------------------------------------

#Parameter range for the trigger and length parameter sensitivity analysis
lengthdata <- expand.grid("int1length" = seq(3,9, by = 3), "int2length" = seq(3,9, by = 3)) #exploring dt as a parameter for scenario 1
lengthdata2345 <- expand.grid("int1length" = seq(6,18, by = 6), "int2length" = seq(6,18, by = 6)) #exploring dt as a parameter for all other scenarios
cminoptim <- expand.grid("cmin1" = seq(0, 1, by = 0.05), "cmin2" = seq(0, 1, by = 0.05))
#Differences in scenario as dt compensations are necessary for scenario 1

#Initial conditions and parameters
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

outcomelistcmin <- list() #Initialise the empty list to store data

for (z in 1:5) { #Iterate through the for loop for every scenario explored
  parms["scen"] = z

  for(j in 1:nrow(lengthdata)) { #Iterate through the 3x3 different dt values explored for intervention 1 and 2
    j = j
    
    #The bottom if function defines whether compensations are necessary for the dt values for scenario 1 compared to all other scenarios
    #Scenario 1 is explored with half the dt length relative to all other scenarios to compensate for a different beta(t) profile
    if(parms["scen"] != 0 && parms["scen"] != 1) {
      parms["t_dur1"] = lengthdata2345[j,1] * 7
      parms["t_dur2"] = lengthdata2345[j,2] * 7
    } else{
      parms["t_dur1"] <- lengthdata[j,1] * 7
      parms["t_dur2"] <- lengthdata[j,2] * 7 
    }
    
    outcomelistcmin[[j]] <- local({
      
      optim <- data.frame(matrix(ncol = 5, nrow = nrow(cminoptim))) #Initialise the empty dataframe to store outcome measures for explore cmin values
      
      for(i in 1:nrow(cminoptim)) { #Specify the cmin1 and 2 values explored, run the model and collect outcome measures
        parms["cmin1"] <- cminoptim[i,1]
        parms["cmin2"] <- cminoptim[i,2]
        
        out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
        optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                       "cmin1" = parms[["cmin1"]], "cmin2" = parms[["cmin2"]])
      }
      
      colnames(optim) <- c("peak", "cum", "scen", "cmin1", "cmin2")
      
      #Generic plotting for the peak prevalence and the attack rate 
      p1 <- ggplot(optim, aes(x = cmin1, y = cmin2, fill= peak))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=18),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.8, "cm"),
              legend.key.width =  unit(2, "cm")) + labs(x = bquote(.(Intervention ~ 1 ~ italic(c[min]))), y = bquote(.(Intervention ~ 2 ~ italic(c[min]))), fill = "I(t) Peak", 
                                                          title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " weeks"))
      
      p2 <- ggplot(optim, aes(x = cmin1, y = cmin2, fill= cum))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=18),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.8, "cm"),
              legend.key.width =  unit(2, "cm")) + 
        labs(x = bquote(.(Intervention ~ 1 ~ italic(c[min]))), y = bquote(.(Intervention ~ 2 ~ italic(c[min]))), fill = "Attack\nRate", 
                                                                      title = paste0("Int 1 & 2 Duration: ", parms["t_dur1"]/7, "/", parms["t_dur2"]/7, " wks"))
      
      
      #Scenario specific plotting for the peak prevalence and the attack rate - making sure that the legend is common for panels for each respective scenario
      if(z == 1) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.04, 0.15, by = (0.15-0.04)/4), limits = c(0.04, 0.15) )  
        p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma", breaks= seq(0.46, 0.8, by = (0.8-0.46)/4), limits = c(0.46, 0.8) )  
      }
      if(z == 2) {
        p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks=seq(0.048, 0.15, by = (0.15-0.048)/4), limits = c(0.048, 0.15) )  
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
      
      print(paste0("Scenario: ", z, " | Length Analysis: ", (j/nrow(lengthdata))*100, "%")) #Progress notification
      
      return(list(p1,p2))
    })   
  }
  
  #For each scenario plot the different dt values explored for the peak prevalence and the attack rate
  peak_multilength_cmin <- ggarrange(outcomelistcmin[[1]][[1]], outcomelistcmin[[2]][[1]], outcomelistcmin[[3]][[1]],
                                     outcomelistcmin[[4]][[1]], outcomelistcmin[[5]][[1]], outcomelistcmin[[6]][[1]],
                                     outcomelistcmin[[7]][[1]], outcomelistcmin[[8]][[1]], outcomelistcmin[[9]][[1]],
                                nrow = 3, ncol =3, align = "hv", common.legend = TRUE, legend = "bottom")
  
  cum_multilength_cmin <- ggarrange(outcomelistcmin[[1]][[2]], outcomelistcmin[[2]][[2]], outcomelistcmin[[3]][[2]],
                                    outcomelistcmin[[4]][[2]], outcomelistcmin[[5]][[2]], outcomelistcmin[[6]][[2]],
                                    outcomelistcmin[[7]][[2]], outcomelistcmin[[8]][[2]], outcomelistcmin[[9]][[2]],
                               nrow = 3, ncol =3, align = "hv", common.legend = TRUE, legend = "bottom")
  
  #Combine the peak prevalence and the attack rate plots together
  combplotcmin <- ggarrange(NULL,NULL, 
                              NULL,peak_multilength_cmin, 
                              NULL,NULL,
                              NULL,cum_multilength_cmin,
                              nrow = 4, ncol = 2, align = "v", heights = c(0.1,1,0.1,1), widths = c(0.05,1),
                              labels = c("","",
                                         "","A",
                                         "","", 
                                         "","B"), font.label = c(size = 45), vjust = -1.1, hjust = 0.8)
  
  #Save the combined peak prevalence and attack rate plots for each scenario
  ggsave(combplotcmin, filename = paste0("heatcum_multicmin_comb",z,".png"), dpi = 300, type = "cairo", width = 15, height = 22, units = "in")
}
