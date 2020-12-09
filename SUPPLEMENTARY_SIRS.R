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

# SIRS MODEL - FUNCTIONS --------------------------------------------------------------
#Function defining beta(t) for the SIRS model 
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
  return(output)
}

plot(seq(0,400),combbeta(1, seq(0,400), 12, 24*7, 0.4))

#SIRS set of ODEs
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, cmin)
    
    dS = - beta*S*(I) + sigma*R
    dI = beta*S*(I) - gamma*I
    dR = gamma*I - sigma*R
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

# SIRS MODEL - Trajectory plot with baseline parameters  -----------------------------------

#Initial conditions and parameters
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,1000,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          sigma = 1/(30*3),
          t_dur = 12*7,
          cmin = 0.4)

#Initiate an empty dataframe
datalist <- list()

for(j in 1:length(seq(1,3))) { #For each of the 5 scenarios conduct a model run
  datalist[[j]] <- local({
    j=j
    data <- data.frame(matrix(nrow = 9, ncol = 0)) #Create empty dataframe 
    #Define the scenario based on the current model iteration
    
    parms["sigma"] <- c(1/(30*3),1/(6*30),1/(12*30))[j]
    
    for(i in 1:2) {
      
      parms["scen"] <- c(0,1)[i]
      
      #Run the model 
      out <- cbind(data.frame(ode(y = init, func = SIRS, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], 
                   "beta" = combbeta(c(0,1)[i], times, parms[["tstart"]], parms[["t_dur"]], parms[["cmin"]]))
      out$re <- out$beta/parms[["gamma"]]*out$S
      data <- rbind(data, out)
      
      print(max(out$I))
    }
    
    #Manipulate the data into a suitable format for plotting 
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    #Identify the peak prevalence and the attack rate for each of the scenarios + text for annotation on the plot
    peak <- round(max(out$I), 3)
    cum <- round(max(out$C), 3)
    datatext <- data.frame(x = c(500, 500), y = c(0.19, 0.16), label = c( paste0("italic(I)[italic(max)]", " ==", peak), 
                                                                           paste0("italic(I)[italic(c)](italic(t)[italic(max)])", " ==", cum)))
    
    #Generic plotting (common across scenarios) for the peak prevalence and the attack rate
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0, 0.2),expand = c(0,0)) + 
      scale_x_continuous( expand = c(0, 0)) + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkred", "darkred"))
    
    #Identifying the region on the plot where the intervention is occuring
    shade <- data.frame(xmin =  parms[["tstart"]], xmax = parms[["tstart"]]+(parms[["t_dur"]]), ymin = 0, ymax = Inf)
    
    #Scenario specific plotting (removing the axis labels and titles for certain panels)

      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
        labs(x ="Time (Days)", y = "Prevalence", col = "", title = paste0(c(3,6,12)[j], " Months")) +
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 4.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
   
      
      if(j == 1) { 
        p1 <- p1 + theme(legend.position = "none", 
                         plot.title = element_text(size = 24, vjust = 3, hjust = 0.5, face = "bold"), axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                         axis.title.y=element_text(size=15), axis.title.x = element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) 
      } else {
          p1 <- p1 + theme(legend.position = "none", 
                           plot.title = element_text(size = 24, vjust = 3, hjust = 0.5, face = "bold"), axis.text.x=element_text(size=12),axis.text.y=element_blank(),
                           axis.title.y=element_blank(), axis.title.x = element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) 
        }  
 
       dump <- list(p1, plotdata, max(out$C)) #Retaining the plots and dataframes needed outside of the local() for-loop
    return(dump)
  })
}

#Combine the plots together for each scenario
traj <- ggarrange(NULL,datalist[[1]][[1]], datalist[[2]][[1]], datalist[[3]][[1]], NULL,nrow = 1, ncol = 5, 
                  legend = "none", align = "h", widths = c(0.1,1,0.9,0.9,0.1))

# Single Intervention - MultiParm Sensitivity Analysis --------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,1000,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 1,
          tstart = 52,
          sigma = 1/(30*3),
          t_dur = 12*7,
          cmin = 0.4)

#Creating two parameter spaces to explore - With the bottom one used for scenario 1
#This is to create a relative scale for scenario 1 in terms of dt
parameterspacescen1 <- expand.grid("trigday" = seq(0,100, by =5), "length" = seq(1,126, by = 2.5))

#Create a empty list to store data
scensens <- list()

#Run for each of the 5 scenarios
for(j in 1:length(seq(1,3))) { #For each of the 5 scenarios conduct a model run
  scensens[[j]] <- local({
    j=j
    data <- data.frame(matrix(nrow = 9, ncol = 0)) #Create empty dataframe 
    #Define the scenario based on the current model iteration
    
    parms["sigma"] <- c(1/(30*3),1/(6*30),1/(12*30))[j]
    
    #If scenario 1 then create a vector "parameterspace" using the relative scale tp/dt parameter space
    parameterspace <- parameterspacescen1
    
    scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 5)) #Create an empty dataframe
    
    for(i in 1:nrow(parameterspace)) { # for each combination of dt and tp parameters identify the peak and attack rate
      
      print(paste0("Scenario ", j," - ", round(i/nrow(parameterspace), digits = 2))) #Progress notification 
      
      #Implement the explored parameter set in the model
      parms["tstart"] <- parameterspace[i,1]
      parms["t_dur"] <- parameterspace[i,2] 

      out <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
      

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
      scale_x_continuous(expand = c(0, 0)) + theme_bw() + scale_fill_viridis_c(direction = -1) + 
      labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "")
    
    p2 <- ggplot(scendata, aes(x = tstart, y = t_dur, fill= cum)) + geom_tile() +
      scale_x_continuous(expand = c(0, 0)) + theme_bw() + scale_fill_viridis_c(direction = -1, option = "magma") + 
      labs(x = bquote("Trigger Point ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "")
    
    
    if(j == 1) {
      p1 <- p1  + scale_y_continuous(expand = c(0,0), labels = formatter2) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_text(size=12), axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    } else {
      p1 <- p1  + scale_y_continuous(expand = c(0,0), labels = formatter2) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    }
    
    if(j == 1) {
      p2 <- p2  + scale_y_continuous(expand = c(0,0), labels = formatter2) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_text(size=12), axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    } else {
      p2 <- p2  + scale_y_continuous(expand = c(0,0), labels = formatter2) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    }
    #Have a combination of the peak prevalence and the attack rate for every scenario
    dump <- list(p1, p2,  scendata)
    return(dump)
  })
}

single <- ggarrange(NULL,scensens[[1]][[1]],scensens[[2]][[1]],scensens[[3]][[1]], NULL,
                    nrow = 1, ncol = 5,align = "h", widths = c(0.1, 1.1,0.9,0.9, 0.1))

singleattack <- ggarrange(NULL, scensens[[1]][[2]],scensens[[2]][[2]],scensens[[3]][[2]], NULL,
          ncol = 5, nrow = 1, align = "h", widths = c(0.1,1.1,0.9,0.9, 0.1)) #peak

# Multi Intervention - MultiParm Sensitivity Analysis ---------------------

combbetamult <- function(scen, time, tstart1, t_dur1, tstart2, t_dur2, cmin1, cmin2) {
  gamma <- 1/GenTime(3, 2.8)
  betascale <- (2.8*gamma)*0.7
  if(scen == 0) {
    output <- betascale
  }
  if(scen == 1) {
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse((time >= (tstart1 + t_dur1 + tstart2)),
                            betascale*cmin2,
                            betascale*cmin1),
                     betascale)
  }
  return(output)
}

plot(seq(0,365),combbetamult(1, seq(0,365), 52, 6*7, 20, 6*7, 0.4, 1))

#Set of ODEs for an SIR model with beta(t)
SIRSmulti <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbetamult(scen, time, tstart1, t_dur1, tstart2, t_dur2, cmin1, cmin2)
    
    dS = - beta*S*(I) + sigma*R
    dI = beta*S*(I) - gamma*I
    dR = gamma*I - sigma*R
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

#Identify the range of values that I want to explore for the trigger point for intervention 1 and 2 - all possible combinations
optimdata <- expand.grid("tstart1" = seq(0,100, by = 5), "tstart2" = seq(0,100, by = 5))

#Define the intial conditions and the parameters
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,1000,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 1,
          sigma = 1/(30*3),
          tstart1 = 52,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          cmin1 = 0.4,
          cmin2 = 0.4)

#Create and empty dataframe to store results
outcomelist <- list()

for(j in 1:3) { #Run for the 5 scenarios 
  j = j
  outcomelist[[j]] <- local({
    
    #Create a empty dataframe to store results 
    optim <- data.frame(matrix(ncol = 6, nrow = nrow(optimdata)))
    
    parms["sigma"] <- c(1/(30*3),1/(6*30),1/(12*30))[j]
    
    for(i in 1:nrow(optimdata)) { #Run for each row of tp1 and tp2 combinations
      parms["tstart1"] <- optimdata[i,1] 
      parms["tstart2"] <- optimdata[i,2]
      
      #Run the model
      out <- data.frame(ode(y = init, func = SIRSmulti, times = times, parms = parms))
      optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                     "tstart1" = parms[["tstart1"]], "tstart2" = parms[["tstart2"]],
                     "realstart2" = parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]])
    }
    
    colnames(optim) <- c("peak", "cum", "scen", "tstart1", "tstart2", "realstart2")
    
    #Basic plot framework for p1 (peak prevalence) and p2 (attack rate)
    p1 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= peak)) + geom_tile() +
      scale_x_continuous(expand = c(0, 0)) + theme_bw() + scale_fill_viridis_c(direction = -1) +
      labs(x = bquote("Trigger Point 1 ("*italic(t[p1])*")"), y = bquote("Trigger Point 2 ("*italic(t[p2])*")"), fill = "", 
           title = "")
    
    p2 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= cum)) + geom_tile() +
      scale_x_continuous(expand = c(0, 0)) + theme_bw() + scale_fill_viridis_c(direction = -1, option = "magma") +
      labs(x = bquote("Trigger Point 1 ("*italic(t[p1])*")"), y = bquote("Trigger Point 2 ("*italic(t[p2])*")"), fill = "", 
           title = "")
    
    #Scenario specific plot alterations (removing axis titles etc)
    
    if(j == 1) {
      p1 <- p1  + scale_y_continuous(expand = c(0,0)) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_text(size=12), axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 5, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    } else {
      p1 <- p1  + scale_y_continuous(expand = c(0,0)) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 5, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    }
    
    if(j == 1) {
      p2 <- p2  + scale_y_continuous(expand = c(0,0)) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_text(size=12), axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 5, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    } else {
      p2 <- p2  + scale_y_continuous(expand = c(0,0)) +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text.x=element_text(size=12), 
              axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_text(size=15),  
              plot.title = element_text(size = 5, vjust = 2, hjust = -0.2, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.45, "cm"),legend.key.width =  unit(1.2, "cm")) 
    }
    print(paste0("Scenario: ", j, " Complete")) #Progress notification
    dump <- list(p1,p2, optim)
    return(dump)
  })
  
}

#Creating combinatory plots of all scenarios for peak prevalence and attack rate
multi <-ggarrange(NULL, outcomelist[[1]][[1]],outcomelist[[2]][[1]],outcomelist[[3]][[1]], NULL,
                  ncol = 5, nrow = 1, align = "h", widths = c(0.1,1.1,0.9,0.9, 0.1)) #peak


multiattack <- ggarrange(NULL, outcomelist[[1]][[2]],outcomelist[[2]][[2]],outcomelist[[3]][[2]], NULL,
          ncol = 5, nrow = 1, align = "h", widths = c(0.1,1.1,0.9,0.9, 0.1)) #peak

# Combining all 3 Sections Together ---------------------------------------

combinplot <- ggarrange(NULL, traj, single, singleattack, multi, multiattack, align = "v", ncol = 1, nrow = 6, labels = c("","A", "B", "C", "D", "E"),
                        font.label = c(size = 28), vjust = 0.1, hjust = -0.4, heights = c(0.1, 1,1,1,1,1))

ggsave(combinplot, filename = "SIRS_Supp.png", dpi = 300, type = "cairo", width = 13, height = 16, units = "in")

