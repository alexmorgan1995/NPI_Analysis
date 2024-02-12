rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur1, scaling) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1),
         (1.7*(1/(GenTime(3.3,2.8))))*scaling, #Intervention
         ifelse((time >= tstart1+tdur1 & time <= 281), 
                (1.7*(1/(GenTime(3.3,2.8)))), #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta1(seq(0,281), 71, (24*7), 0.5))

beta2 <- function(time, tstart1, tdur1, scaling) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         (1.7*(1/(GenTime(3.3,2.8)))), #Intervention
         ifelse((time >= tstart1+tdur1 & time <= 281), 
                (1.7*(1/(GenTime(3.3,2.8)))), #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta2(seq(0,281), 71, (24*7)))

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta1(time,tstart1,tdur1,scaling)*Iv*Sv - beta1(time,tstart1,tdur1,scaling)*Inv*Sv + zeta*Rv
    
    dSnv = - beta2*Inv*Snv - beta2*Iv*Snv + zeta*Rnv
    
    dIv = beta1(time,tstart1,tdur1,scaling)*Iv*Sv + beta1(time,tstart1,tdur1,scaling)*Inv*Sv - gamma*Iv
    
    dInv =  beta2*Inv*Snv + beta2*Iv*Snv - gamma*Inv
    
    dRv = gamma*Iv - zeta*Rv
    
    dRnv = gamma*Inv - zeta*Rnv
    
    dRv = gamma*Iv - zeta*Rv
    
    dCv = beta1(time,tstart1,tdur1,scaling)*Iv*Sv + beta1(time,tstart1,tdur1,scaling)*Inv*Sv

    return(list(c(dSv, dSnv, dIv, dInv, dRv, dRnv, dCv)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####

#Initial Conditions and Times

init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, Rv= 0, Rnv = 0, Cv = 0)
times <- seq(0,281,by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur1 = 24*7,
          scaling = 0, #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED
          beta2 = 1.7*(1/(GenTime(3.3,2.8))))

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
out1$Iv <- out1$Iv/0.20
out1$Inv <- out1$Inv/0.80
out1$Rv <- out1$Rv/0.20
out1$Rnv <- out1$Rnv/0.80
out1$Cv <- out1$Cv/0.20
out1$Beta1 <- beta1(seq(0,281), 71, (24*7), 0) #VARY THIS ASWELL
out1$Beta2 <- 1.7*(1/(GenTime(3.3,2.8)))
colnames(out1) <- c("Time", "Suscv", "Suscnv","Infected_Iv", "Infected_Inv", "Recovv", "Recovnv", "Cv", "Beta1", "Beta2")

#### Introducing the Intervention ####

scaling  <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
scalingdata <- data.frame(matrix(nrow = 0, ncol = 7))

for(i in 1:length(scaling)) {
  temp <- data.frame(matrix(nrow = length(seq(0,281)), ncol = 7))
  init <- c(Sv = 0.2, Snv = (0.8)-0.0001, Iv = 0, Inv = 0.0001, Rv= 0, Rnv = 0, Cv = 0)
  times <- seq(0,281,by = 1)
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            zeta = 1/365,
            tstart1 = 71, 
            tdur1 = 24*7,
            scaling = scaling[i],
            beta2 = 1.7*(1/(GenTime(3.3,2.8))))
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$Iv <- out1$Iv/0.2
  out1$Inv <- out1$Inv/0.8
  temp[,1] <- as.factor((1-scaling[i])*100)
  temp[,2] <- out1$time
  temp[,3] <- out1$Iv
  temp[,4] <- out1$Inv
  temp[,5] <- out1$Rv + out1$Rnv
  temp[,6] <- beta1(seq(0,281), 71, (24*7), scaling[i])
  temp[,7] <- 1.7*(1/(GenTime(3.3,2.8)))
  scalingdata <- rbind.data.frame(scalingdata, temp)
}

colnames(scalingdata) <- c("Intervention_Magnitude", "Time", "Infected_Iv", "Infected_Inv", "Recov", "Beta1", "Beta2")

statsinfecv <- melt(scalingdata, id.vars = c("Intervention_Magnitude", "Time"), measure.vars = c("Infected_Iv"))
statsinfecv$Intervention_Magnitude <- factor(statsinfecv$Intervention_Magnitude, levels=c("0", "20", "40", "80", "100", "60"))

statsinfecnv <- melt(scalingdata, id.vars = c("Intervention_Magnitude", "Time"), measure.vars = c("Infected_Inv"))
statsinfecnv$Intervention_Magnitude <- factor(statsinfecnv$Intervention_Magnitude, levels=c("0", "20", "40", "80", "100", "60"))

statsrecov <- melt(scalingdata, id.vars = c("Intervention_Magnitude", "Time"), measure.vars = c("Recov"))
statsrecov$Intervention_Magnitude <- factor(statsrecov$Intervention_Magnitude, levels=c("0", "20", "40", "80", "100", "60"))

statsbeta1 <- melt(scalingdata, id.vars =  c("Intervention_Magnitude", "Time"), measure.vars = c("Beta1"))
statsbeta1$Intervention_Magnitude <- factor(statsbeta1$Intervention_Magnitude, levels=c("0", "20", "40", "80", "100", "60"))

#### Aggregated Plots ####

pinfv <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.125) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=11),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + labs(color='Intervention Efficacy (%)') + 
  annotate("rect", xmin = as.numeric(parms[3]), xmax = as.numeric(parms[3])+24*7, ymin = 0, ymax = 0.125, alpha = .1, fill = "darkred") + 
  scale_color_manual(breaks = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
                    values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                    labels = c("0", "0.2", "0.4","0.6", "0.8","1"))

pinfnv <- ggplot(data = statsinfecnv, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.125) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=11),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  labs(color='Intervention Efficacy (%)')  + 
  annotate("rect", xmin = as.numeric(parms[3]), xmax = as.numeric(parms[3])+24*7, ymin = 0, ymax = 0.125, alpha = .1, fill = "darkred") + 
  scale_color_manual(breaks = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                     labels = c("0", "0.2", "0.4","0.6", "0.8","1"))

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=11), axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  annotate("rect", xmin = as.numeric(parms[3]), xmax = as.numeric(parms[3])+24*7, ymin = 0, ymax = 1, alpha = .1, fill = "darkred") + 
  scale_color_manual(breaks = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                     labels = c("0", "0.2", "0.4","0.6", "0.8","1"))

pbeta <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.35) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + labs(color='Intervention Efficacy') + 
  scale_color_manual(breaks = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                     labels = c("0", "0.2", "0.4","0.6", "0.8","1"))

ggarrange(pinfv,
          pinfnv,
          prec,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))


#### Seperated Plots ####

ggplot(data = statsinfecv, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.125) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + labs(color='Intervention Efficacy (%)') + 
  annotate("rect", xmin = as.numeric(parms[3]), xmax = as.numeric(parms[3])+24*7, ymin = 0, ymax = 0.125, alpha = .1, fill = "darkred") + labs(color='Intervention Efficacy') +  
  scale_color_manual(breaks = c("0", "20", "40", "60", "80", "100"),
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                     labels = c("0", "20", "40", "60", "80", "100"))

ggplot(data = statsinfecnv, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.125) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  labs(color='Intervention Efficacy (%)')  + 
  annotate("rect", xmin = as.numeric(parms[3]), xmax = as.numeric(parms[3])+24*7, ymin = 0, ymax = 0.125, alpha = .1, fill = "darkred") + labs(color='Intervention Efficacy (%)') + 
  scale_color_manual(breaks = c("0", "20", "40", "60", "80", "100"),
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                     labels = c("0", "20", "40", "60", "80", "100"))

ggplot(data = statsrecov, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14), axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + labs(color='Intervention Efficacy (%)') + 
  annotate("rect", xmin = as.numeric(parms[3]), xmax = as.numeric(parms[3])+24*7, ymin = 0, ymax = 1, alpha = .1, fill = "darkred") + 
  scale_color_manual(breaks = c("0", "20", "40", "60", "80", "100"),
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                     labels = c("0", "20", "40", "60", "80", "100"))

ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.35) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + labs(color='Intervention Efficacy (%)') + 
  scale_color_manual(breaks = c("0", "20", "40", "60", "80", "100"),
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F"), 
                     labels = c("0", "20", "40", "60", "80", "100"))

