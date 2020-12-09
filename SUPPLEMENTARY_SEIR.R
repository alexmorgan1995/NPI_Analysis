library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis/supplementary")
start_time <- Sys.time()

# Common Functions --------------------------------------------------------
#Generation Time Function
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

# SEIR MODEL - FUNCTIONS --------------------------------------------------------------
#Function defining beta(t) for the SEIR model 
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
SEIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, cmin)
    
    dS = - beta*S*(I)
    dE = beta*S*(I)- sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dE, dI, dR, dC)))
  })
} 

# SEIR MODEL - Trajectory plot with baseline parameters  -----------------------------------

#Initial conditions and parameters
init <- c(S = 0.99999, E = 0, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          sigma = 1/3,
          t_dur = 12*7,
          cmin = 0.4)

#Initiate an empty dataframe
datalist <- list()

for(j in 1:length(seq(1,5))) { #For each of the 5 scenarios conduct a model run
  datalist[[j]] <- local({
    j=j
    data <- data.frame(matrix(nrow = 9, ncol = 0)) #Create empty dataframe 
    explor_scen <- c(0, seq(1,5)[j]) #Define the scenario based on the current model iteration
    
    for(i in 1:2) {
      
      parms["scen"] <- explor_scen[i]
      if(parms["scen"] != 0 && parms["scen"] != 1) { #if the scenario is not 1 or 0 then double dt
        parms["t_dur"] = 24*7
      }
      
      #Run the model 
      out <- cbind(data.frame(ode(y = init, func = SEIR, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], 
                   "beta" = combbeta(explor_scen[i], times, parms[["tstart"]], parms[["t_dur"]], parms[["cmin"]]))
      out$re <- out$beta/parms[["gamma"]]*out$S
      data <- rbind(data, out)
    }
    
    #Manipulate the data into a suitable format for plotting 
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    #Identify the peak prevalence and the attack rate for each of the scenarios + text for annotation on the plot
    peak <- round(max(out$I), 3)
    cum <- round(max(out$C), 3)
    datatext <- data.frame(x = c(200, 200), y = c(0.15, 0.125), label = c( paste0("italic(I)[italic(max)]", " ==", peak), 
                                                                           paste0("italic(I)[italic(c)](italic(t)[italic(max)])", " ==", cum)))
    
    #Generic plotting (common across scenarios) for the peak prevalence and the attack rate
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.16),expand = c(0,0)) + 
      scale_x_continuous( expand = c(0, 0)) + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkred", "darkred"))
    
    #Identifying the region on the plot where the intervention is occuring
    shade <- data.frame(xmin =  parms[["tstart"]], xmax = parms[["tstart"]]+(parms[["t_dur"]]), ymin = 0, ymax = Inf)
    
    #Scenario specific plotting (removing the axis labels and titles for certain panels)
    if(parms[["scen"]]  == 1) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", col = "", title = "Scenario 1") +
        theme(legend.position = "none", 
              plot.title = element_text(size = 24, vjust = 3, hjust = 0.5, face = "bold"), axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
              axis.title.y=element_text(size=15), axis.title.x = element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 4.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
    }
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "", col = "", title = "Scenario 2") +
        theme(legend.position = "none", 
              plot.title = element_text(size = 24, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x=element_text(size=15),axis.text.y=element_blank(),
              axis.title.y=element_text(size=15), axis.title.x = element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 4.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")

    }
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "", col = "", title = "Scenario 3") +
        theme(legend.position = "none", 
              plot.title = element_text(size = 24, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x=element_text(size=15),axis.text.y=element_blank(),
              axis.title.y=element_text(size=15), axis.title.x = element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 4.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")

    }
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "", col = "", title = "Scenario 4")  +
        theme(legend.position = "none",  
              plot.title = element_text(size = 24, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x =element_text(size=15),axis.text.y=element_blank(),
              axis.title.y=element_text(size=15), axis.title.x = element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 4.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")

    }
    if(parms[["scen"]]  == 5) {
      
      shade <- data.frame(xmin =  c(parms[["tstart"]], parms[["tstart"]]+(parms[["t_dur"]]*0.333), 
                                    parms[["tstart"]]+(parms[["t_dur"]]*0.667)),
                          xmax = c(parms[["tstart"]]+(parms[["t_dur"]]*0.1667), parms[["tstart"]]+(parms[["t_dur"]]*0.53), 
                                   parms[["tstart"]]+(parms[["t_dur"]]*0.833)),
                          ymin = 0, ymax = Inf)
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin,  ymax = ymax, 
                                                              xmin = xmin, xmax = xmax), alpha = 0.2, fill = "darkblue") + geom_line(size = 1.1, stat = "identity") +
        labs(x ="Time (Days)", y = "", col = "",  title = "Scenario 5") +
        theme(legend.position = "none",  
              plot.title = element_text(size = 24, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x=element_text(size=15),axis.text.y=element_blank(),
              axis.title.y=element_text(size=15), axis.title.x = element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 4.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")

    }

    dump <- list(p1, plotdata, max(out$C)) #Retaining the plots and dataframes needed outside of the local() for-loop
    return(dump)
  })
}

#Combine the plots together for each scenario
traj <- ggarrange(NULL,datalist[[1]][[1]], datalist[[2]][[1]], datalist[[3]][[1]], datalist[[4]][[1]], datalist[[5]][[1]], NULL,nrow = 1, ncol = 7, 
                      legend = "none", align = "h", widths = c(0.1,1,0.9,0.9,0.9,0.9,0.1))

# Single Intervention - MultiParm Sensitivity Analysis ------------------------

init <- c(S = 0.99999, E = 0, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          sigma = 1/3,
          t_dur = 12*7,
          cmin = 0.4)

#Creating two parameter spaces to explore - With the bottom one used for scenario 1
#This is to create a relative scale for scenario 1 in terms of dt
parameterspaceOG <- expand.grid("trigday" = seq(0,125, by =5), "length" = seq(1,252, by =5))
parameterspacescen1 <- expand.grid("trigday" = seq(0,125, by =5), "length" = seq(1,126, by = 2.5))

#Create a empty list to store data
scensens <- list()

#Run for each of the 5 scenarios
for(j in 1:5) {
  
  scensens[[j]] = local({
    
    i = 0
    parms["scen"] <- j
    
    #If scenario 1 then create a vector "parameterspace" using the relative scale tp/dt parameter space
    if(parms["scen"] == 1) {
      parameterspace <- parameterspacescen1
    } else{parameterspace <- parameterspaceOG }
    
    scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 5)) #Create an empty dataframe
    
    for(i in 1:nrow(parameterspace)) { # for each combination of dt and tp parameters identify the peak and attack rate
      
      print(paste0("Scenario ", j," - ", round(i/nrow(parameterspace), digits = 2))) #Progress notification 
      
      #Implement the explored parameter set in the model
      parms["tstart"] <- parameterspace[i,1]
      parms["t_dur"] <- parameterspace[i,2] 
      
      out <- data.frame(ode(y = init, func = SEIR, times = times, parms = parms))
      scendata[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                        "tstart" = parms[["tstart"]], "t_dur" = parms[["t_dur"]])
    }
    
    colnames(scendata) <- c("peak", "cum", "scen", "tstart", "t_dur")
    
    #This function is used to convert the scenario 1 dt axis into a relative one (after having used a more limited dt range)
    formatter2 <- function(x){ 
      x*2
    }
    
    #Generic Plotting for Peak Prevalence
    p1 <- ggplot(scendata, aes(x = tstart, y = t_dur, fill= peak)) + geom_tile() +
      scale_x_continuous(expand = c(0, 0)) + theme_bw()
    
    #Changing the axis for scenario 1
    if(parms["scen"] == 1){
      p1 <- p1 + scale_fill_viridis_c(direction = -1, breaks = c(0.07, 0.08, 0.09, 0.1)) +
        labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "") +
        scale_y_continuous(expand = c(0,0), labels = formatter2) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),legend.key.width =  unit(1, "cm")) 
    } else {
      p1 <- p1 + scale_fill_viridis_c(direction = -1) + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "") + 
      scale_y_continuous(expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x =element_text(size=12),axis.text.y =element_blank(),
            axis.title.y=element_blank(), axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),legend.key.width =  unit(1, "cm"))}
    
    #Generic Plotting for Attack Rate
    p2 <- ggplot(scendata, aes(x = tstart, y = t_dur, fill= cum)) + geom_tile() +
      scale_x_continuous(expand = c(0, 0)) + theme_bw()
    
    #Changing the axis for scenario 1
    if(parms["scen"] == 1){
      p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 2, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),legend.key.width =  unit(1, "cm")) 
    } else {
      p2 <- p2 + scale_fill_viridis_c(direction = -1, option = "magma")  + 
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x =element_text(size=12),axis.text.y =element_blank(),
              axis.title.y=element_blank(), axis.title.x = element_text(size=15),  plot.title = element_text(size = 2, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),legend.key.width =  unit(1, "cm"))}
    
    
    #Changing the axis for scenario 1
    if(parms["scen"] == 1){
      p2 <- p2 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "", title = "") +
        scale_y_continuous(expand = c(0,0), labels = formatter2)
    } else { p2 <- p2 + labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "", title = "") + 
      scale_y_continuous(expand = c(0,0))}
    
    #Have a combination of the peak prevalence and the attack rate for every scenario
    dump <- list(p1,p2,  scendata)
    return(dump)
  })
}

single <- ggarrange(scensens[[1]][[1]],scensens[[2]][[1]],scensens[[3]][[1]],scensens[[4]][[1]],scensens[[5]][[1]],
                      nrow = 1, ncol = 5,align = "h", widths = c(1.2,0.9,0.9,0.9,0.9))

singleattack <- ggarrange(scensens[[1]][[2]],scensens[[2]][[2]],scensens[[3]][[2]],scensens[[4]][[2]],scensens[[5]][[2]],
                    nrow = 1, ncol = 5,align = "h", widths = c(1.2,0.9,0.9,0.9,0.9))
# Multi Intervention - MultiParm Sensitivity Analysis ---------------------

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

#Set of ODEs for an SIR model with beta(t)
SIRmulti <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbetamult(scen, time, tstart1, t_dur1, tstart2, t_dur2, cmin1, cmin2)
    
    dS = - beta*S*(I)
    dE = beta*S*(I)- sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dE, dI, dR, dC)))
  })
} 


#Identify the range of values that I want to explore for the trigger point for intervention 1 and 2 - all possible combinations
optimdata <- expand.grid("tstart1" = seq(0,125, by = 5), "tstart2" = seq(0,125, by = 5))

#Define the intial conditions and the parameters
init <- c(S = 0.99999, E = 0, I = 0.00001, R = 0, C = 0)
times <- seq(0,1000,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 1,
          sigma = 1/3,
          tstart1 = 52,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          cmin1 = 0.4,
          cmin2 = 0.4)

#Create and empty dataframe to store results
outcomelist <- list()

for(j in 1:5) { #Run for the 5 scenarios 
  j = j
  
  parms["scen"] <- j #Set the explored NPI scenario for the model 
  
  if(parms["scen"] != 0 && parms["scen"] != 1) { #Double the explored dt for any scenario that is not 1 or 0 
    parms["t_dur1"] = 12*7
    parms["t_dur2"] = 12*7
  }
  
  outcomelist[[j]] <- local({
    
    #Create a empty dataframe to store results 
    optim <- data.frame(matrix(ncol = 6, nrow = nrow(optimdata)))
    
    for(i in 1:nrow(optimdata)) { #Run for each row of tp1 and tp2 combinations
      parms["tstart1"] <- optimdata[i,1] 
      parms["tstart2"] <- optimdata[i,2]
      
      #Run the model
      out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
      optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                     "tstart1" = parms[["tstart1"]], "tstart2" = parms[["tstart2"]],
                     "realstart2" = parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]])
    }
    
    colnames(optim) <- c("peak", "cum", "scen", "tstart1", "tstart2", "realstart2")
    
    #Basic plot framework for p1 (peak prevalence) and p2 (attack rate)
    p1 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= peak))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw()  + 
      scale_fill_viridis_c(direction = -1) +
      labs(x = bquote("Trigger Point 1 ("*italic(t[p1])*")"), y = bquote("Trigger Point 2 ("*italic(t[p2])*")"), fill = "", 
           title = paste("Scenario", j))
    
    p2 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= cum))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() + 
      scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.5, 0.8, by = (0.8-0.5)/3), limits = c(0.5, 0.8))  + 
      labs(x = bquote("Trigger Point 1 ("*italic(t[p1])*")"), y = bquote("Trigger Point 2 ("*italic(t[p2])*")"), fill = "", title = "")
    
    #Scenario specific plot alterations (removing axis titles etc)
    if(parms[["scen"]] == 1) {
      p1 <- p1 + theme(legend.position = "bottom", legend.title = element_text(size=12), legend.text=element_text(size=15), axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                       axis.title.y=element_text(size=15),axis.title.x = element_text(size=15), plot.title = element_blank(),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),
                       legend.key.width =  unit(1, "cm")) + guides(colour = guide_colourbar(title.vjust = 0.9))
    }
    
    else{
      p1 <- p1 + theme(legend.position = "bottom", legend.title = element_text(size=12), legend.text=element_text(size=15), axis.text.x=element_text(size=12),axis.text.y=element_blank(),
                       axis.title.y=element_blank(),axis.title.x = element_text(size=15),  plot.title = element_blank(),
                       legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),
                       legend.key.width =  unit(1, "cm")) 

    }
    
    if(parms[["scen"]] == 1) {
      p2 <- p2 + theme(legend.position = "bottom", legend.title = element_text(size=12), legend.text=element_text(size=15), axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                       axis.title.y=element_text(size=15),axis.title.x = element_text(size=15), plot.title = element_blank(),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),
                       legend.key.width =  unit(1, "cm")) + guides(colour = guide_colourbar(title.vjust = 0.9))
    }
    
    else{
      p2 <- p2 + theme(legend.position = "bottom", legend.title = element_text(size=12), legend.text=element_text(size=15), axis.text.x=element_text(size=12),axis.text.y=element_blank(),
                       axis.title.y=element_blank(),axis.title.x = element_text(size=15),  plot.title = element_blank(),
                       legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.3, "cm"),
                       legend.key.width =  unit(1, "cm")) 
      
    }
    print(paste0("Scenario: ", j, " Complete")) #Progress notification
    dump <- list(p1,p2, optim)
    return(dump)
  })
  
}

#Creating combinatory plots of all scenarios for peak prevalence and attack rate
multi <-ggarrange(outcomelist[[1]][[1]],outcomelist[[2]][[1]],outcomelist[[3]][[1]],outcomelist[[4]][[1]],outcomelist[[5]][[1]],
               ncol = 5, nrow = 1, align = "h", widths = c(1.2,0.9,0.9,0.9,0.9)) #peak

multiattack <-ggarrange(outcomelist[[1]][[2]],outcomelist[[2]][[2]],outcomelist[[3]][[2]],outcomelist[[4]][[2]],outcomelist[[5]][[2]],
                  ncol = 5, nrow = 1, align = "h", widths = c(1.2,0.9,0.9,0.9,0.9)) #

# Combining all 3 Sections Together ---------------------------------------

combinplot <- ggarrange(NULL, traj, single, singleattack, multi, multiattack, align = "v", ncol = 1, nrow = 6, labels = c("","A", "B", "C", "D", "E"),
                        font.label = c(size = 28), vjust = 0.1, hjust = -0.4, heights = c(0.1, 1,1,1,1,1))

ggsave(combinplot, filename = "SEIR_Supp.png", dpi = 300, type = "cairo", width = 14, height = 16, units = "in")
