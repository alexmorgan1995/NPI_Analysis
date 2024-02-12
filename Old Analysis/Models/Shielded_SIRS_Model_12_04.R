rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot");library("reshape2"); library("dplyr"); library("RColorBrewer"); library("ggpubr")

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur1, scaling) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         (1.7*(1/(GenTime(3.3,2.8))))*scaling,
         ifelse((time >= tstart1+tdur1 & time <= 281), 
                (1.7*(1/(GenTime(3.3,2.8)))), 
                (1.7*(1/(GenTime(3.3,2.8))))
         )
  )
}


plot(beta1(seq(0,281), 71, (24*7), 0.5))

beta2 <- function(time, tstart1, tdur1, scaling) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         (1.7*(1/(GenTime(3.3,2.8)))),
         ifelse((time >= tstart1+tdur1 & time <= 281), 
                (1.7*(1/(GenTime(3.3,2.8)))), 
                (1.7*(1/(GenTime(3.3,2.8))))
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
    
    dCv = 

    return(list(c(dSv, dSnv, dIv, dInv, dRv, dRnv)))
  })
}

#### Testing the Model Structure ####

#Initial Conditions and Times

init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, Rv= 0, Rnv = 0)
times <- seq(0,281,by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur1 = 6*7,
          scaling = 0.5,
          beta2 = 1.7*(1/(GenTime(3.3,2.8))))

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
out1$Iv <- out1$Iv/0.20
out1$Inv <- out1$Inv/0.80
out1$Rv <- out1$Rv/0.20
out1$Rnv <- out1$Rnv/0.80
out1$Beta1 <- beta1(seq(0,281), 71, (24*7), 0.5)
out1$Beta2 <- 1.7*(1/(GenTime(3.3,2.8)))
  
ggplot(out1, aes(x= time, y = Iv)) + geom_line()

colnames(out1) <- c("Time", "Suscv", "Suscnv","Infected_Iv", "Infected_Inv", "Recovv", "Recovnv", "Beta1", "Beta2")

#### Introducing the Intervention ####

scaling  <- seq(0,1, by= 0.2)
scalingdata <- data.frame(matrix(nrow = 0, ncol = 7))

for(i in 1:length(scaling)) {
  temp <- data.frame(matrix(nrow = length(seq(0,281)), ncol = 7))
  init <- c(Sv = 0.2, Snv = (0.8)-0.0001, Iv = 0, Inv = 0.0001, Rv= 0, Rnv = 0)
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
  print((i/length(scaling))*100)
}

colnames(scalingdata) <- c("Intervention_Magnitude", "Time", "Infected_Iv", "Infected_Inv", "Recov", "Beta1", "Beta2")

statsinfecv <- melt(scalingdata, id.vars = c("Intervention_Magnitude", "Time"), measure.vars = c("Infected_Iv"))
statsinfecnv <- melt(scalingdata, id.vars = c("Intervention_Magnitude", "Time"), measure.vars = c("Infected_Inv"))
statsrecov <- melt(scalingdata, id.vars = c("Intervention_Magnitude", "Time"), measure.vars = c("Recov"))

statsbeta1 <- melt(scalingdata, id.vars =  c("Intervention_Magnitude", "Time"), measure.vars = c("Beta1"))
statsbeta2 <- melt(scalingdata, id.vars =  c("Intervention_Magnitude", "Time"), measure.vars = c("Beta2"))

# Aggregated Plots

ggplot(data = statsinfecv, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Population Infected") + scale_y_continuous(limits = c(0,0.125) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette="Paired") + labs(color='Intervention Efficacy (%)') 

ggplot(data = statsinfecnv, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Population Infected") + scale_y_continuous(limits = c(0,0.125) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette="Paired") + labs(color='Intervention Efficacy (%)') 

ggplot(data = statsrecov, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette="Paired") + labs(color='Intervention Efficacy (%)') 

ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.35) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  +
  scale_color_brewer(palette="Paired") + labs(color='Intervention Efficacy (%)') 

ggplot(data = statsbeta2, aes(x = (Time), y = value, col = Intervention_Magnitude)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.35) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), axis.text=element_text(size=13), 
        axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  +
  scale_color_brewer(palette="Paired") + labs(color='Intervention Efficacy (%)') 

