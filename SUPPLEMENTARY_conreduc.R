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

# MODELLING SUSTAINABLE REDUCTIONS TO TRANSMISSION AFTER NPIS MODEL - FUNCTIONS ----------------------------------------------------------

#Function describing beta(t) for an initial optimisable NPI followed by a constant, indefinite reduction to transmission
#Note that the chracteristics of the constant reduction to transmission are defined by a start date and a magnitude - contstart and cmin2

combbetamult <- function(scen, time, tstart1, t_dur1, contstart, cmin1, cmin2) {
  gamma <- 1/GenTime(3, 2.8)
  betascale <- (2.8*gamma)*0.7
  if(scen == 0) {
    output <-  betascale
  }
  if(scen == 1) {
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (contstart) & time <= (contstart + Inf))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            betascale*cmin1,
                            betascale*cmin2),
                     betascale)
  }
  if(scen == 2) {
    betalin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(cmin1, 1), method="linear", rule =2)
    output <- ifelse((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (contstart) & time <= (contstart + Inf)),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            betalin1(time)*betascale,
                     cmin2*betascale),
                            betascale)
  }
  if(scen == 3) {
    betalin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(1, cmin1), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (contstart) & time <= (contstart + Inf))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            betalin1(time)*betascale,
                            cmin2*betascale),
                     betascale)
  }
  if(scen == 4) {
    betaincdec1 <-approxfun(x=c(tstart1, tstart1+(t_dur1/2), tstart1+t_dur1), y= c(1, cmin1, 1), method="linear", rule =2)
  output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (contstart) & time <= (contstart + Inf))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                           betaincdec1(time)*betascale,
                           cmin2*betascale),
                     betascale)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart1 & time <= tstart1+(t_dur1*0.16667)) | (time >= tstart1+(t_dur1*0.3333) & time <= tstart1+(t_dur1*0.5)) |
                       (time >= tstart1+(t_dur1*0.6667) & time <= tstart1+(t_dur1*0.83333)) |
                       (time >= (contstart) & time <= (contstart + Inf)),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            cmin1*betascale,
                            cmin2*betascale),
                     betascale)
  }
  return(output)
}

#Example plotting of beta(t)
plot(seq(0,200),combbetamult(3, seq(0,200), 20, 2*7, 100, 0.1, 0.5))

#SIR set of ODEs
SIRmulti <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbetamult(scen, time, tstart1, t_dur1, contstart, cmin1, cmin2)
    
    dS = - beta*S*I
    dI = beta*S*I- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

# MODELLING SUSTAINABLE REDUCTIONS TO TRANSMISSION AFTER NPIS MODEL - Trajectory plots and sensitivity analyses --------------------------------------------------

#Range of values explored in the sensitivity analysis
sens <- list("tstart1" = seq(0,100, by = 1),
             "cmin1" = seq(0,1, by = 0.01),
             "t_dur1" = seq(1,100, by = 1))

#Initial conditions and parameter values
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,300,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 1,
          tstart1 = 52,
          t_dur1 = 6*7,
          contstart = 94,
          cmin1 = 0.4,
          cmin2 = 0.5)

#Empty list for trajectory plot and sensitivity analysis data storage
outcomelistcmin <- list()

#This analysis loop both conducts a basic trasejectory plot for I(t) and beta(t) and
#Also explores the 


for(z in 1:5) { #Run for each of the 5 scenarios and store each trajectory plot separately
  outcomelistcmin[[z]] <- local({
    
    parms[["scen"]] <- z #Define the scenario model parameter for each analysis loop 
    
    ##TRAJECTORY PLOT ANALYSIS
    #Run the model and obtain I(t) and beta(t)
    outtraj <- cbind(data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms)),
                 "beta" =  combbetamult(parms[["scen"]], times, parms[["tstart1"]], parms[["t_dur1"]], parms[["contstart"]], parms[["cmin1"]], parms[["cmin2"]]))
    
    #Identify the peak prevalence and attack rate and define a text box for each scenario with these two outcome measures
    peak <- round(max(outtraj$I), 3)
    cum <- round(max(outtraj$C), 3)
    datatext <- data.frame(x = c(150, 150), y = c(0.07, 0.06), label = c(paste0("italic(I)[italic(max)]", " ==", peak), 
                                                                            paste0("italic(I)[italic(c)](italic(t)[italic(max)])", " ==", cum)))
    
    ##SENSITIVITY ANALYSIS
    #Create an empty dataframe to store sensitivity analysis results 
    datasens <- data.frame(matrix(nrow = 0, ncol = 5))
    
    for(i in 1:3) { #Run for each parameter explored 
      sensitivity <- sens[[i]] #Define what parameter you are exploring 
      data <- data.frame(matrix(nrow = 0, ncol = 5))
      
      for(j in 1:length(sensitivity)) { #For the identified parameter - run the model for the explored range
        parms1 <- parms #create a dummy vector for the model parameters 
        parms1[names(sens)[i]] <- sensitivity[j] #For the explored parameter and value - assign it to the actual parameter vector
        out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms1)) #Run the model
        data[j,1] <- max(out$I)
        data[j,2] <- max(out$C)
        data[j,3] <- sensitivity[j]
        data[j,4] <- names(sens)[i]
        data[j,5] <- as.character(z)
      }
      datasens <- rbind(datasens, data)
    }
    
    colnames(datasens) <- c("peak", "cum", "sensval", "sens", "scen")
    
    print(outtraj) # Progress bar
    
    #Define the period of the intervention for both the NPI (blue) and the constant reduction to transmission (red)
    shade <- data.frame(xmin =  c(parms[["tstart1"]], parms[["contstart"]]), 
                        xmax = c(parms[["tstart1"]]+parms[["t_dur1"]], Inf), ymin = c(0,0), ymax = c(Inf,Inf))
    
    #Generic Plotting for the trajectory plot and the beta(t) plot
    p1 <- ggplot(data = outtraj, aes(x = time, y = I)) + theme_bw() + scale_y_continuous(limits = c(0,0.075), expand = c(0,0)) + 
      scale_x_continuous(expand = c(0, 0)) + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                                                       fill = c("darkblue", "darkred")) + geom_line(size = 1.1, stat = "identity") 
    
    p2 <- ggplot(outtraj, aes(x = time, y = beta))  + theme_bw() + scale_y_continuous(limits = c(0 , 0.3), expand = c(0,0)) + 
    geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
              fill = c("darkblue", "darkred")) + geom_line(size = 1.1, stat = "identity") + 
      geom_line(size = 1.1, stat = "identity", col = "darkblue") + scale_x_continuous(expand = c(0, 0))
    
    
    #Scenario specific plotting for the peak prevalence and attack rate plots (removing axis labels for certain panels)
    if(parms[["scen"]] == 1) { 
           p1 <- p1 + theme(axis.text.x = element_blank(), axis.text.y = element_text(size=13), axis.title.y= element_text(size=15), 
                            axis.title.x = element_blank(), plot.margin=unit(c(0.2,0.4,0.2,0.4),"cm"),  
                            plot.title = element_text(size = 15, vjust = 3, hjust = 0.5, face = "bold")) + 
             labs(y = "Prevalence", title = paste("Scenario", z), x = "Time") + 
             geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
           p2 <- p2 + theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), axis.title.y=element_text(size=15),
                            axis.title.x = element_text(size=15), plot.margin=unit(c(0.2,0.4,0.2,0.4),"cm")) + 
             labs(y = bquote(italic(beta*"(t)")), x = "Time")
    } else { 
      p1 <- p1 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y= element_blank(), 
                       axis.title.x = element_blank(), plot.margin=unit(c(0.2,0.4,0.2,0.4),"cm"),  
                       plot.title = element_text(size = 15, vjust = 3, hjust = 0.5, face = "bold")) + 
        labs(y = bquote(italic(I*"(t)")), title = paste("Scenario", z))+ 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + theme(axis.text.x = element_text(size=13), axis.text.y = element_blank(), axis.title.y= element_blank(), 
                       axis.title.x = element_text(size=15), plot.margin=unit(c(0.2,0.4,0.2,0.4),"cm")) + 
        labs(x = "Time") 
    }
    
    combplot <- ggarrange(p1,p2, nrow = 2, ncol = 1, common.legend = TRUE, legend = "none", align = "v", heights = c(1, 0.6))
  
    dump <- list(combplot, datasens)
    #ggarange pase all scenarios together
    return(dump)
  })
}

#Combining the scenarios for the trajectory + beta(t) plots together
traj <- ggarrange(NULL,outcomelistcmin[[1]][[1]], outcomelistcmin[[2]][[1]], outcomelistcmin[[3]][[1]], outcomelistcmin[[4]][[1]], outcomelistcmin[[5]][[1]],NULL,
                  nrow = 1, ncol =7, widths = c(0.1,1,0.75,0.75,0.75,0.75,0.1))

#Formatting the sensitivity analysis for plotting
sensitivityanal <- data.frame(matrix(nrow = 0, ncol = 5))
for(i in 1:5) {
  sensitivityanal <- rbind(sensitivityanal, outcomelistcmin[[i]][[2]])
}

#Need to plot the sensitivity analysis for the explord parameter for the peak prevalence (p1) and the attack rate (p2) for each scenario

#MAGNITUDE (c_min)
cminpeak <- ggplot(sensitivityanal[sensitivityanal$sens == "cmin1",], aes(x = sensval, y = peak, col = scen)) + geom_line(size = 1.02) + theme_bw()  + 
  theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18), legend.spacing.x = unit(0.3, 'cm'),
        axis.text.x=element_blank(), axis.text.y=element_text(size=13), axis.title.y=element_text(size=15), axis.title.x = element_blank(), 
        plot.margin=unit(c(0.2,0.4,0.2,0.1),"cm")) + labs(y = expression(italic("I"["max"])), col = "Scenario") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.15),expand = c(0, 0)) + 
  scale_color_manual(values=c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))
  

cmincum <- ggplot(sensitivityanal[sensitivityanal$sens == "cmin1",], aes(x = sensval, y = cum, col = scen)) + geom_line(size = 1.02) + theme_bw()  + 
  theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18), legend.spacing.x = unit(0.3, 'cm'),
        axis.text=element_text(size=13),axis.title.x = element_text(size=15), 
        plot.margin=unit(c(0.2,0.4,0.2,0.1),"cm")) + 
  labs(x = expression("Magnitude "~"("~italic("c"["min"]~")")), y = expression(italic("I"["c"]~"("~tmax~")")), col = "Scenario")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.75),expand = c(0, 0))+ 
  scale_color_manual(values=c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))

#TSTART (t_p)
tstartpeak <- ggplot(sensitivityanal[sensitivityanal$sens == "tstart1",], aes(x = sensval, y = peak, col = scen)) + geom_line(size = 1.02) + theme_bw()  + 
  theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18), legend.spacing.x = unit(0.3, 'cm'),
        axis.text.x=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_blank(), 
        plot.margin=unit(c(0.2,0.4,0.2,0.1),"cm")) + 
  labs(y = expression(italic("I"["max"])), col = "Scenario") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.15),expand = c(0, 0))+ 
  scale_color_manual(values=c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))

tstartcum <- ggplot(sensitivityanal[sensitivityanal$sens == "tstart1",], aes(x = sensval, y = cum, col = scen)) + geom_line(size = 1.02) + theme_bw()  + 
  theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18), legend.spacing.x = unit(0.3, 'cm'),
        axis.text.x=element_text(size=13),axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_text(size=15), 
        plot.margin=unit(c(0.2,0.4,0.2,0.1),"cm")) + 
  labs(x = expression("Trigger "~"("~italic("t"["p"]~")")), y = expression(italic("I"["c"]~"("~infinity~")")), col = "Scenario")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.75),expand = c(0, 0))+ 
  scale_color_manual(values=c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))

#DURATION (d_t)
durpeak <- ggplot(sensitivityanal[sensitivityanal$sens == "t_dur1",], aes(x = sensval, y = peak, col = scen)) + geom_line(size = 1.02) + theme_bw()  + 
  theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18), legend.spacing.x = unit(0.3, 'cm'),
        axis.text.x=element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_blank(), 
        plot.margin=unit(c(0.2,0.4,0.2,0.1),"cm")) + 
  labs(y = expression(italic("I"["max"])), col = "Scenario")+
  scale_x_continuous(limits = c(0,100), expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.15),expand = c(0, 0))+ 
  scale_color_manual(values=c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))

durcum <- ggplot(sensitivityanal[sensitivityanal$sens == "t_dur1",], aes(x = sensval, y = cum, col = scen)) + geom_line(size = 1.02) + theme_bw()  + 
  theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18), legend.spacing.x = unit(0.3, 'cm'),
        axis.text.x=element_text(size=13),axis.text.y=element_blank(), axis.title.y=element_blank(),axis.title.x = element_text(size=15), 
        plot.margin=unit(c(0.2,0.4,0.2,0.1),"cm")) + 
  labs(x = expression("Duration "~"("~italic("d"["t"]~")")), y = expression(italic("I"["c"]~"("~infinity~")")), col = "Scenario")+
  scale_x_continuous(limits = c(0,100), expand = c(0, 0)) + scale_y_continuous(limits = c(0,0.75), expand = c(0, 0))+ 
  scale_color_manual(values=c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))

#Creating a unifying plot for the sensitivity analysis
peak <- ggarrange(NULL,cminpeak, tstartpeak, durpeak,NULL,
                  NULL,cmincum, tstartcum, durcum,NULL,
                  nrow = 2, ncol =5, common.legend = TRUE, legend = "bottom", align = "v", widths = c(0.1,1,1,1,0.1), heights = c(0.8,1))

#Combining the trajectory/beta(t) + sensitivity analysis plots together 
comb <- ggarrange(NULL, traj, NULL,peak, nrow = 4, ncol = 1, heights = c(0.08,1,0.1,0.8),
                  labels = c("","A","", "B"),font.label = c(size = 25), vjust = -0.3, hjust = -0.4)

#Save the plot to the working directory
ggsave(comb, filename = "5_scenarios_multiadapt_const.png", dpi = 300, type = "cairo", width = 11, height = 10, units = "in")


