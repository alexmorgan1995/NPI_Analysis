library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")
start_time <- Sys.time()

# Model Functions ----------------------------------------------------------
#Function for the generation time - a function of R0 and the doubling time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function describing beta(t) for the 5 scenarios

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
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

# Epidemic Trajectory Plots for the Multiple Intervention Case Study ----------------------------------------

#Initial conditions and parameter values
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart1 = 52,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          cmin1 = 0.4,
          cmin2 = 0.4)

#Empty list for data storage
datalist <- list()

#For loop for the plotting of epidemic trajectory plots, beta(t) and Re(t) for the 5 scenarios
for(j in 1:length(seq(1,5))) {
  datalist[[j]] <- local({
    j=j
    
    data <- data.frame(matrix(nrow = 9, ncol = 0)) #Creating an empty dataframe
    explor_scen <- c(0, seq(1,5)[j]) #Select both unmitigated and the explored intervention scenario for each trajectory plot
    
    for(i in 1:2) {
      
      parms["scen"] <- explor_scen[i]
      
      if(parms["scen"] != 0 && parms["scen"] != 1) { #If the current is not 1 or 0, then double the dt value to compensate
        parms["t_dur1"] = 12*7
        parms["t_dur2"] = 12*7
      }
      
      #Run the model and obtain key data
      out <- cbind(data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], 
                   "beta" = combbetamult(explor_scen[i], times, parms[["tstart1"]], parms[["t_dur1"]], parms[["tstart2"]], parms[["t_dur2"]],
                                       parms[["cmin1"]], parms[["cmin2"]]))
      out$re <- (out$beta/parms["gamma"])*out$S
      data <- rbind(data, out)
    }
    
    #Melt the data into an appropriate format for plotting 
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    plotbeta <- melt(data, id.vars = c("time", "group"), measure.vars = ("beta"))
    plotre <- melt(data, id.vars = c("time", "group"), measure.vars = ("re"))
    
    #Obtaining the peak prevalence and the attack rate for plotting 
    peak <- round(max(out$I), 3); cum <- round(max(out$C), 3)
    datatext <- data.frame(x = c(200, 200), y = c(0.1875, 0.17), label = c( paste0("italic(I)[italic(max)]", " ==", peak), 
                                                                           paste0("italic(I)[italic(c)](italic(t)[italic(max)])", " ==", cum)))
    
    #Plotting the generic plotting format for plot 1 (traj plot), plot 2 (beta(t)) and plot 3 (re(t))
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.2),expand = c(0,0)) +
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
    
    #Create areas of the plot where the intervention is occuring
    shade <- data.frame(xmin =  c(parms[["tstart1"]], parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]]), 
                        xmax = c(parms[["tstart1"]]+parms[["t_dur1"]], parms[["tstart1"]]+parms[["t_dur1"]]+parms[["tstart2"]]+parms[["t_dur2"]]), 
                        ymin = 0, ymax = Inf)
    
    #Plot specific alterations (removing axis titles and labels etc)
    if(parms[["scen"]]  == 1) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", col = "", title = "Scenario 1")+ 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = bquote(beta* "(t)"), col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = expression(R[e](t)), col = "")
    }
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 2")+ 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 3")+ 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 4")+ 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    if(parms[["scen"]]  == 5) {
      
      #Scenario 5 requires a more complicated shading to denote the intervention period 
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
        labs(x ="", y = "", col = "",  title = "Scenario 5")+ 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "")
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "")
    }
    
    #Combining the 3 plots together for each scenario
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.45, 0.5))
    dump <- list(combplot, plotdata, max(out$C))
    return(dump)
  })
}

#Combining the 5 scenarios together 
combplotmulti <- ggarrange(datalist[[1]][[1]], datalist[[2]][[1]], datalist[[3]][[1]], datalist[[4]][[1]], datalist[[5]][[1]],nrow = 1, ncol = 5, 
                           legend = "none", align = "h")

ggsave(combplotmulti, filename = "5_scenarios_multi.png", dpi = 300, type = "cairo", width = 18, height = 10, units = "in")

#Obtaining the exact optimums for the sensitivity analysis
for(i in 1:5) {
  print(datalist[[i]][[2]][datalist[[i]][[2]]$group == "scenario",]
        [which.max(datalist[[i]][[2]]$value[datalist[[i]][[2]]$group == "scenario"]),])
  print(datalist[[i]][[3]])
}

for(i in 1:5) {
  print(datalist[[i]][[2]][datalist[[i]][[2]]$group == "scenario",]
        [which.max(datalist[[i]][[2]]$value[datalist[[i]][[2]]$group == "scenario"]),])
  print(datalist[[i]][[3]])
}

# Multiple Optimisations - Trigger Point --------------------------------------------------

#Identify the range of values that I want to explore for the trigger point for intervention 1 and 2 - all possible combinations
optimdata <- expand.grid("tstart1" = seq(0,100, by = 5), "tstart2" = seq(0,100, by = 5))

#Define the intial conditions and the parameters
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,1000,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 1,
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
      scale_fill_viridis_c(direction = -1, breaks=seq(0.04, 0.152, by = (0.152-0.04)/4), limits = c(0.04, 0.152)) +
      labs(x = bquote("Trigger Point 1 ("*italic(t[p1])*")"), y = bquote("Trigger Point 2 ("*italic(t[p2])*")"), fill = "I(t) Peak\n", 
           title = paste("Scenario", j))
    
    p2 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= cum))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() + 
      scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.5, 0.8, by = (0.8-0.5)/3), limits = c(0.5, 0.8))  + 
      labs(x = bquote("Trigger Point 1 ("*italic(t[p1])*")"), y = bquote("Trigger Point 2 ("*italic(t[p2])*")"), fill = "Attack\nRate", title = "")
    
    #Scenario specific plot alterations (removing axis titles etc)
    if(parms[["scen"]] != 5) {
      p1 <- p1 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_blank(),axis.text.y=element_text(size=18),
                       axis.title.y=element_text(size=18),axis.title.x = element_blank(),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(1, "cm")) + guides(colour = guide_colourbar(title.vjust = 0.9))
      p2 <- p2 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_blank(),axis.text.y=element_blank(),
                       axis.title.y=element_blank(),axis.title.x = element_blank(),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(0.8, "cm")) 
    }
    
    else{
      p1 <- p1 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_text(size=18),axis.text.y=element_text(size=18),
                       axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(1, "cm")) 
      
      p2 <- p2 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_text(size=18),axis.text.y=element_blank(),
                       axis.title.y=element_blank(),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(0.8, "cm")) 
    }
    print(paste0("Scenario: ", j, " Complete")) #Progress notification
    dump <- list(p1,p2, optim)
   return(dump)
  })

}

#Creating combinatory plots of all scenarios for peak prevalence and attack rate
p1 <-ggarrange(NULL, outcomelist[[1]][[1]],outcomelist[[2]][[1]],outcomelist[[3]][[1]],outcomelist[[4]][[1]],outcomelist[[5]][[1]],
               ncol = 1, nrow = 6, align = "v", common.legend = TRUE, legend = "bottom", heights = c(0.2,0.8,0.8,0.8,0.8,1)) #peak
p2 <-ggarrange(NULL,outcomelist[[1]][[2]],outcomelist[[2]][[2]],outcomelist[[3]][[2]],outcomelist[[4]][[2]],outcomelist[[5]][[2]],
               ncol = 1, nrow = 6, align = "v", common.legend = TRUE, legend = "bottom", heights = c(0.2,0.8,0.8,0.8,0.8,1)) #cum

#Combing the two sets of plots together  
multicombplot <- ggarrange(p1, NULL, p2, ncol = 3, nrow = 1, widths = c(1.05,0.03,0.85), align = "h")

# Multiple Optimisations - CMIN --------------------------------------------------

#Create sensitivity analysis range for c_min1 and c_min2 parameters - all possible combinations
cminoptim <- expand.grid("cmin1" = seq(0, 1, by = 0.05), "cmin2" = seq(0, 1, by = 0.05))

#Initial conditions and parameter values 
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,1000,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 1,
          tstart1 = 52,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          cmin1 = 0.4,
          cmin2 = 0.4)

outcomelistcmin <- list()

for(j in 1:5) { #Run for each of the 5 scenarios 
  j = j
  
  parms["scen"] <- j
  
  if(parms["scen"] != 0 && parms["scen"] != 1) { #If scenario is not 1 or 0 then double the dt 
    parms["t_dur1"] = 12*7
    parms["t_dur2"] = 12*7
  }
  
  outcomelistcmin[[j]] <- local({
    
    optim <- data.frame(matrix(ncol = 5, nrow = nrow(cminoptim))) #Creating empty dataframe
    
    for(i in 1:nrow(cminoptim)) { #Define cmin1 and cmin2 based on the iterated row of the sensitivity analysis dataframe
      parms["cmin1"] <- cminoptim[i,1]
      parms["cmin2"] <- cminoptim[i,2]
      
      #Run the model
      out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
      optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                     "cmin1" = parms[["cmin1"]], "cmin2" = parms[["cmin2"]])
    }
    
    colnames(optim) <- c("peak", "cum", "scen", "cmin1", "cmin2")
    
    #Generic plotting shared across all scenarios for the peak prevalence and the attack rate 
    p1 <- ggplot(optim, aes(x = cmin1, y = cmin2, fill= peak))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      labs(x = bquote(.(Intervention ~ 1 ~ italic(c[min]))), y = bquote(.(Intervention ~ 2 ~ italic(c[min]))), fill = "I(t) Peak\n", 
                                                        title = paste("Scenario", j)) + 
      scale_fill_viridis_c(direction = -1, breaks=seq(0.04, 0.152, by = (0.152-0.04)/4), limits = c(0.04, 0.152))
    
    p2 <- ggplot(optim, aes(x = cmin1, y = cmin2, fill= cum))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      labs(x = bquote(.(Intervention ~ 1 ~ italic(c[min]))), y = bquote(.(Intervention ~ 2 ~ italic(c[min]))), fill = "Attack\nRate", title = "") + 
      scale_fill_viridis_c(direction = -1, option = "magma", breaks=seq(0.48, 0.8, by = (0.8-0.48)/4), limits = c(0.48, 0.8))
    
    #Scenari-specific alterations to the plots (removal of plot labels and axis labels etc.)
    if(parms[["scen"]] != 5) {
      p1 <- p1 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_blank(),axis.text.y=element_text(size=18),
                       axis.title.y=element_text(size=18),axis.title.x = element_blank(),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(1, "cm"))
      
      p2 <- p2 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_blank(),axis.text.y=element_blank(),
                       axis.title.y=element_blank(),axis.title.x = element_blank(),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(0.8, "cm")) 
    } else{
      p1 <- p1 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_text(size=18),axis.text.y=element_text(size=18),
                       axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(1, "cm")) 
      
      p2 <- p2 + theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=13), axis.text.x=element_text(size=18),axis.text.y=element_blank(),
                       axis.title.y=element_blank(),axis.title.x = element_text(size=18),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
                       legend.spacing.x = unit(0.5, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.5, "cm"),
                       legend.key.width =  unit(0.8, "cm")) 
    }
    print(paste0("Scenario: ", j, " Complete"))
    dump <- list(p1,p2, optim)
    return(dump)
  })   
}

#Combine the plots for the peak prevalence and the attack rate 
pcmin1 <-ggarrange(NULL, outcomelistcmin[[1]][[1]],outcomelistcmin[[2]][[1]],outcomelistcmin[[3]][[1]],outcomelistcmin[[4]][[1]],outcomelistcmin[[5]][[1]],
               ncol = 1, nrow = 6, align = "v", common.legend = TRUE, legend = "bottom", heights = c(0.2,0.8,0.8,0.8,0.8,1)) #peak
pcmin2 <-ggarrange(NULL,outcomelistcmin[[1]][[2]],outcomelistcmin[[2]][[2]],outcomelistcmin[[3]][[2]],outcomelistcmin[[4]][[2]],outcomelistcmin[[5]][[2]],
               ncol = 1, nrow = 6, align = "v", common.legend = TRUE, legend = "bottom", heights = c(0.2,0.8,0.8,0.8,0.8,1)) #cum

#Combine the cmin plots together 
multicombplotcmin <- ggarrange(pcmin1, NULL, pcmin2, ncol = 3, nrow = 1, widths = c(1.05,0.03,0.85), align = "h")

end_time <- Sys.time()
end_time - start_time

# Combplot ----------------------------------------------------------------

#Combine the tp1/tp2 and cmin1/cmin2 plots together for figure 3
combplotcminpeak <- ggarrange(NULL, multicombplot, NULL, multicombplotcmin,NULL, ncol = 5, nrow = 1, widths = c(0.1,1,0.1,1.05,0.1), align = "h", labels = c("","A", "", "B",""), font.label = c(size = 35),
                              vjust = 1.2, hjust = 0.8, common.legend = TRUE)

#Save the final figure 3 in the working directory
ggsave(combplotcminpeak, filename = "Heat_5_comb.png", dpi = 300, type = "cairo", width = 15, height = 16, units = "in")
