rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur, scaling1) {
  ifelse((time >= tstart1 & time <= tstart1+tdur),
         (0.6*(1/(GenTime(3.3,2.8)))), #Intervention
         ifelse((time >= tstart1+tdur & time <= 730), 
                (0.6*(1/(GenTime(3.3,2.8))))*scaling1, #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta1(seq(0, 720), 71, (24*7), 0.4))

beta2 <- function(time, tstart1, tdur, scaling2) {
  ifelse((time >= tstart1 & time <= tstart1+tdur),
         (0.6*(1/(GenTime(3.3,2.8)))), #Intervention
         ifelse((time >= tstart1+tdur & time <= 730), 
                (1.7*(1/(GenTime(3.3,2.8))))*scaling2, #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta2(seq(0, 720), 71, (24*7), 0.4))

beta3 <- function(time, tstart1, tdur, scaling3) {
  ifelse((time >= tstart1 & time <= tstart1+tdur),
         (0.6*(1/(GenTime(3.3,2.8)))), #Intervention
         ifelse((time >= tstart1+tdur & time <= 730), 
                (2.8*(1/(GenTime(3.3,2.8))))*scaling3, #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta3(seq(0, 720), 71, (6*7), 0.8))


#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

    if(time < tstart1) {
      dSv = - beta1(time,tstart1,tdur,scaling1)*(Iv+Ih+Inv)*Sv + zeta*Rv
      
      dSh = - beta2(time,tstart1,tdur,scaling2)*(Iv+Ih+Inv)*Sh + zeta*Rh
      
      dSnv = - beta3(time,tstart1,tdur,scaling3)*(Iv+Ih+Inv)*Snv  + zeta*Rnv
      
      dIv = beta1(time,tstart1,tdur,scaling1)*(Iv+Ih+Inv)*Sv - gamma*Iv
      
      dIh =  beta2(time,tstart1,tdur,scaling2)*(Iv+Ih+Inv)*Sh - gamma*Ih
      
      dInv = beta3(time,tstart1,tdur,scaling3)*(Iv+Ih+Inv)*Snv - gamma*Inv
      
      dRv = gamma*Iv - zeta*Rv
      
      dRh = gamma*Ih - zeta*Rh
      
      dRnv = gamma*Inv - zeta*Rnv
    }
    else{
      dSv = - beta1(time,tstart1,tdur,scaling1)*Iv*Sv - beta1(time,tstart1,tdur,scaling1)*Ih*Sv + zeta*Rv
      
      dSh = - beta2(time,tstart1,tdur,scaling2)*Ih*Sh - beta1(time,tstart1,tdur,scaling2)*Iv*Sh - 
        beta2(time,tstart1,tdur,scaling2)*Inv*Sh + zeta*Rh
      
      dSnv = - beta3(time,tstart1,tdur,scaling3)*Inv*Snv - beta2(time,tstart1,tdur,scaling3)*Ih*Snv + zeta*Rnv
      
      
      
      dIv = beta1(time,tstart1,tdur,scaling1)*Iv*Sv + beta1(time,tstart1,tdur,scaling1)*Ih*Sv - gamma*Iv
      
      dIh =  beta2(time,tstart1,tdur,scaling2)*Ih*Sh + beta1(time,tstart1,tdur,scaling2)*Iv*Sh + 
        beta2(time,tstart1,tdur,scaling2)*Inv*Sh - gamma*Ih
      
      dInv = beta3(time,tstart1,tdur,scaling3)*Inv*Snv + beta2(time,tstart1,tdur,scaling3)*Ih*Snv - gamma*Inv
      
      
  
      dRv = gamma*Iv - zeta*Rv
      
      dRh = gamma*Ih - zeta*Rh
      
      dRnv = gamma*Inv - zeta*Rnv
    }
    return(list(c(dSv, dSh, dSnv, dIv, dIh, dInv, dRv, dRh, dRnv)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####

#Initial Conditions and Times

init <- c(Sv = 0.2, Sh = 0.2, Snv = 0.6-0.0001, Iv = 0, Ih = 0, Inv = 0.0001, Rv= 0, Rh = 0, Rnv = 0)
times <- seq(0,730,by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling1 = 0.4,
          scaling2 = 0.4,
          scaling3 = 0.4) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))

out1$Sv <- out1$Sv/0.20; out1$Sh <- out1$Sh/0.2; out1$Snv <- out1$Snv/0.60
out1$Iv <- out1$Iv/0.20; out1$Ih <- out1$Ih/0.2; out1$Inv <- out1$Inv/0.60
out1$Rv <- out1$Rv/0.20; out1$Rh <- out1$Rh/0.2; out1$Rnv <- out1$Rnv/0.60

out1$Beta1 <- beta1(seq(0,730), 71, (6*7), 0.4) #VARY THIS ASWELL
out1$Beta2 <- beta2(seq(0,730), 71, (6*7), 0.4)
out1$Beta3 <- beta3(seq(0,730), 71, (6*7), 0.4)

colnames(out1) <- c("Time", "Suscv", "Susch", "Suscnv", "Infected_Iv", "Infected_Ih", "Infected_Inv", 
                    "Recovv", "Recovh", "Recovnv", "Beta1", "Beta2", "Beta3")

ggplot(out1, aes(x = Time, y = Infected_Iv)) + geom_line()
#### Exploring the Parameters ####

scaling <- seq(0,1, by = 0.2)

scalingdata <- data.frame(matrix(nrow = 0, ncol = 7))

for(i in 1:length(scaling)) { 
  init <- c(Sv = 0.2, Sh = 0.2, Snv = 0.6-0.0001, Iv = 0, Ih = 0, Inv = 0.0001, Rv= 0, Rh = 0, Rnv = 0)
  times <- seq(0,730,by = 1)
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            zeta = 1/365,
            tstart1 = 71, 
            tdur = 6*7,
            scaling1 = scaling[i],
            scaling2 = (1.7*(1/(GenTime(3.3,2.8))) - ((1.7*(1/(GenTime(3.3,2.8))) - 0.6*(1/(GenTime(3.3,2.8))))*(1-scaling[i])))/
              (1.7*(1/(GenTime(3.3,2.8)))),
            scaling3 = (2.8*(1/(GenTime(3.3,2.8))) - ((2.8*(1/(GenTime(3.3,2.8))) - 1.7*(1/(GenTime(3.3,2.8))))*(1-scaling[i])))/
              (2.8*(1/(GenTime(3.3,2.8)))))
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$Sv <- out1$Sv/0.20; out1$Sh <- out1$Sh/0.2; out1$Snv <- out1$Snv/0.60
  out1$Iv <- out1$Iv/0.20; out1$Ih <- out1$Ih/0.2; out1$Inv <- out1$Inv/0.60
  out1$Rv <- out1$Rv/0.20; out1$Rh <- out1$Rh/0.2; out1$Rnv <- out1$Rnv/0.60
  out1$Beta1 <- beta1(seq(0,730), 71, (6*7), parms[5]) #VARY THIS ASWELL
  out1$Beta2 <- beta2(seq(0,730), 71, (6*7), parms[6])
  out1$Beta3 <- beta3(seq(0,730), 71, (6*7), parms[7])
  out1$efficacy <- as.factor(1-scaling[i])
  scalingdata <- rbind.data.frame(scalingdata, out1)
}


colnames(scalingdata) <- c("Time", "Suscv", "Susch", "Suscnv", "Infected_Iv", "Infected_Ih", "Infected_Inv", 
                    "Recovv", "Recovh", "Recovnv", "Beta1", "Beta2", "Beta3", "Efficacy")

statsinfecv <- melt(scalingdata, id.vars = c("Time", "Efficacy"), measure.vars = c("Infected_Iv"))
statsinfech <- melt(scalingdata, id.vars = c("Time","Efficacy"), measure.vars = c("Infected_Ih"))
statsinfecnv <- melt(scalingdata, id.vars = c("Time", "Efficacy"), measure.vars = c("Infected_Inv"))

statsbeta1 <- melt(scalingdata, id.vars =  c("Time", "Efficacy"), measure.vars = c("Beta1"))
statsbeta2 <- melt(scalingdata, id.vars =  c("Time", "Efficacy"), measure.vars = c("Beta2"))
statsbeta3 <- melt(scalingdata, id.vars =  c("Time", "Efficacy"), measure.vars = c("Beta3"))

#### Aggregated Plots ####

ggplot(data = statsinfecv, aes(x = (Time), y = value, col = Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired")
  
ggplot(data = statsinfech, aes(x = (Time), y = value, col = Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Shielders Infected") + scale_y_continuous(limits = c(0,0.125),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired")

ggplot(data = statsinfecnv, aes(x = (Time), y = value, col = Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.125),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired")

ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired")

ggplot(data = statsbeta2, aes(x = (Time), y = value, col = Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired")

ggplot(data = statsbeta3, aes(x = (Time), y = value, col = Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??3") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired")
