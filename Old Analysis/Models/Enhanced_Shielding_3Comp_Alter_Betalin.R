rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

(2.8*(1/(GenTime(3.3,2.8))))

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur, scaling) {
  beta1_2 <- (0.6*(1/(GenTime(3.3,2.8)))) * (1-scaling)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.6*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.6*(1/(GenTime(3.3,2.8))),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(1/(GenTime(3.3,2.8))))))}

plot(beta1(seq(0,730), 71, (6*7), 1))

beta2 <- function(time, tstart1, tdur, scaling) {
  beta1_2 <- (1.7*(1/(GenTime(3.3,2.8))) - ((1.7*(1/(GenTime(3.3,2.8))) - 0.6*(1/(GenTime(3.3,2.8))))*scaling))
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.6*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.6*(1/(GenTime(3.3,2.8))),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(1/(GenTime(3.3,2.8))))))}

plot(beta2(seq(0,730), 71, (6*7), 1))

beta3 <- function(time, tstart1, tdur, scaling) {
  beta1_2 <- (2.8*(1/(GenTime(3.3,2.8))) - ((2.8*(1/(GenTime(3.3,2.8))) - 1.7*(1/(GenTime(3.3,2.8))))*scaling))
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.6*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.6*(1/(GenTime(3.3,2.8))),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(1/(GenTime(3.3,2.8))))))}

plot(beta3(seq(0,730), 71, (6*7), 1))

beta4 <- function(time, tstart1,) {
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0, 
         1.7*(1/(GenTime(3.3,2.8)))
  )
}

plot(beta4(seq(0,730), 71))

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    if(time < tstart1) { #Same but introduce Ci
      dSv = - beta1(time,tstart1,tdur,scaling)*Iv*Sv - beta1(time,tstart1,tdur,scaling)*Ih*Sv -
        beta1(time,tstart1,tdur,scaling)*cb*In*Sv - beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Iv*Sv - 
        beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Ih*Sv + zeta*Rv
      
      dSh = - beta2(time,tstart1,tdur,scaling)*Ih*Sh - beta2(time,tstart1,tdur,scaling)*In*Sh -
        beta1(time,tstart1,tdur,scaling)*Iv*Sh + zeta*Rh
      
      dSn = - beta3(time,tstart1,tdur,scaling)*In*Sn - beta2(time,tstart1,tdur,scaling)*Ih*Sn -
        beta1(time,tstart1,tdur,scaling)*cb*Iv*Sn - beta3(time,tstart1,tdur,scaling)*bn*(1-cb)*In*Sn + zeta*Rn
      
      dIv = beta1(time,tstart1,tdur,scaling)*Iv*Sv + beta1(time,tstart1,tdur,scaling)*Ih*Sv +
        beta1(time,tstart1,tdur,scaling)*cb*In*Sv + beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Iv*Sv + 
        beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Ih*Sv - gamma*Iv
      
      dIh = beta2(time,tstart1,tdur,scaling)*Ih*Sh + beta2(time,tstart1,tdur,scaling)*In*Sh +
        beta1(time,tstart1,tdur,scaling)*Iv*Sh - gamma*Ih
      
      dIn = beta3(time,tstart1,tdur,scaling)*In*Sn + beta2(time,tstart1,tdur,scaling)*Ih*Sn +
        beta1(time,tstart1,tdur,scaling)*cb*Iv*Sn + beta3(time,tstart1,tdur,scaling)*bn*(1-cb)*In*Sn - gamma*In
      
      dRv = gamma*Iv - zeta*Rv
      
      dRh = gamma*Ih - zeta*Rh
      
      dRn = gamma*In - zeta*Rn
    }
    else{
      dSv = - beta1(time,tstart1,tdur,scaling)*Iv*Sv - beta1(time,tstart1,tdur,scaling)*Ih*Sv -
        beta1(time,tstart1,tdur,scaling)*ci*In*Sv - beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Iv*Sv - 
        beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Ih*Sv  + zeta*Rv
      
      dSh = - beta2(time,tstart1,tdur,scaling)*Ih*Sh - beta2(time,tstart1,tdur,scaling)*In*Sh -
        beta1(time,tstart1,tdur,scaling)*Iv*Sh + zeta*Rh
      
      dSn = - beta3(time,tstart1,tdur,scaling)*In*Sn - beta2(time,tstart1,tdur,scaling)*Ih*Sn -
        beta1(time,tstart1,tdur,scaling)*ci*Iv*Sn - beta3(time,tstart1,tdur,scaling)*bn*(1-cb)*In*Sn + zeta*Rn
      
      dIv = beta1(time,tstart1,tdur,scaling)*Iv*Sv + beta1(time,tstart1,tdur,scaling)*Ih*Sv +
        beta1(time,tstart1,tdur,scaling)*ci*In*Sv + beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Iv*Sv + 
        beta1(time,tstart1,tdur,scaling)*bv*0.5*(1-cb)*Ih*Sv - gamma*Iv
      
      dIh = beta2(time,tstart1,tdur,scaling)*Ih*Sh + beta2(time,tstart1,tdur,scaling)*In*Sh +
        beta1(time,tstart1,tdur,scaling)*Iv*Sh - gamma*Ih
      
      dIn = beta3(time,tstart1,tdur,scaling)*In*Sn + beta2(time,tstart1,tdur,scaling)*Ih*Sn +
        beta1(time,tstart1,tdur,scaling)*ci*Iv*Sn + beta3(time,tstart1,tdur,scaling)*bn*(1-cb)*In*Sn - gamma*In
      
      dRv = gamma*Iv - zeta*Rv
      
      dRh = gamma*Ih - zeta*Rh
      
      dRn = gamma*In - zeta*Rn
    }
    return(list(c(dSv, dSh, dSn, dIv, dIh, dIn, dRv, dRh, dRn)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####

#Initial Conditions and Times

init <- c(Sv = 0.2 - (0.0001*0.2), Sh = 0.2 - (0.0001*0.2), Sn = 0.6 - (0.0001*0.6), Iv = 0.0001*0.2, Ih = 0.0001*0.2, In = 0.0001*0.6, Rv= 0, Rh = 0, Rn = 0)
times <- seq(0,730,by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.4,
          cb = 0.5,
          ci = 0, 
          bv = 0.6/0.2,
          bn = 0.2/0.6) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))

out1$Sv <- out1$Sv/0.20; out1$Sh <- out1$Sh/0.2; out1$Sn <- out1$Sn/0.60
out1$Iv <- out1$Iv/0.20; out1$Ih <- out1$Ih/0.2; out1$In <- out1$In/0.60
out1$Rv <- out1$Rv/0.20; out1$Rh <- out1$Rh/0.2; out1$Rn <- out1$Rn/0.60

out1$Beta1 <- beta1(seq(0,730), 71, (6*7), 0.5) #VARY THIS ASWELL
out1$Beta2 <- beta2(seq(0,730), 71, (6*7), 0.5)
out1$Beta3 <- beta3(seq(0,730), 71, (6*7), 0.5)

colnames(out1) <- c("Time", "Suscv", "Susch", "Suscnv", "Infected_Iv", "Infected_Ih", "Infected_Inv", 
                    "Recovv", "Recovh", "Recovnv", "Beta1", "Beta2", "Beta3")


ggplot(out1, aes(x = Time, y = Infected_Iv)) + geom_line()


#### For Loop ####

scaling <- seq(0,1, by = 0.2)

scalingdata <- data.frame(matrix(nrow = 0, ncol = 13))

for(i in 1:length(scaling)) {
  init <- c(Sv = 0.2 - (0.0001*0.2), Sh = 0.2 - (0.0001*0.2), Sn = 0.6 - (0.0001*0.6), Iv = 0.0001*0.2, Ih = 0.0001*0.2, In = 0.0001*0.6, Rv= 0, Rh = 0, Rn = 0)
  times <- seq(0,730,by = 1)
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            zeta = 1/365,
            tstart1 = 71, 
            tdur = 6*7,
            scaling = scaling[i],
            cb = 0.5,
            ci = 0, 
            bv = 0.6/0.2,
            bn = 0.2/0.6) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$Sv <- out1$Sv/0.20; out1$Sh <- out1$Sh/0.2; out1$Sn <- out1$Sn/0.60
  out1$Iv <- out1$Iv/0.20; out1$Ih <- out1$Ih/0.2; out1$In <- out1$In/0.60
  out1$Rv <- out1$Rv/0.20; out1$Rh <- out1$Rh/0.2; out1$Rn <- out1$Rn/0.60
  out1$Beta1 <- beta1(seq(0,730), 71, (6*7), scaling[i]) #VARY THIS ASWELL
  out1$Beta2 <- beta2(seq(0,730), 71, (6*7), scaling[i])
  out1$Beta3 <- beta3(seq(0,730), 71, (6*7), scaling[i])
  out1$scaling <- as.factor(scaling[i])
  scalingdata <- rbind.data.frame(scalingdata, out1)
}

colnames(scalingdata) <- c("Time", "Suscv", "Susch", "Suscnv", "Infected_Iv", "Infected_Ih", "Infected_In", 
                    "Recovv", "Recovh", "Recovnv", "Beta1", "Beta2", "Beta3", "Scaling")

statsinfecv <- melt(scalingdata, id.vars = c("Time", "Scaling"), measure.vars = c("Infected_Iv"))
statsinfech <- melt(scalingdata, id.vars = c("Time","Scaling"), measure.vars = c("Infected_Ih"))
statsinfecn <- melt(scalingdata, id.vars = c("Time", "Scaling"), measure.vars = c("Infected_In"))

statsbeta1 <- melt(scalingdata, id.vars =  c("Time", "Scaling"), measure.vars = c("Beta1"))
statsbeta2 <- melt(scalingdata, id.vars =  c("Time", "Scaling"), measure.vars = c("Beta2"))
statsbeta3 <- melt(scalingdata, id.vars =  c("Time", "Scaling"), measure.vars = c("Beta3"))

#### Aggregated Plots ####

pinfv <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Vulnerable Infected") + scale_y_continuous(limits = c(0,0.035),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired") + labs(colour = "Efficacy")
  
pinfh <- ggplot(data = statsinfech, aes(x = (Time), y = value, col = Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Shielders Infected") + scale_y_continuous(limits = c(0,0.125),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired") + labs(colour = "Efficacy")

pinfn <- ggplot(data = statsinfecn, aes(x = (Time), y = value, col = Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.15),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired") + labs(colour = "Efficacy")

pbeta1 <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired") + labs(colour = "Efficacy")

pbeta2 <- ggplot(data = statsbeta2, aes(x = (Time), y = value, col = Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired") + labs(colour = "Efficacy")

pbeta3 <- ggplot(data = statsbeta3, aes(x = (Time), y = value, col = Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??3") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette ="Paired") + labs(colour = "Efficacy")


ggarrange(ggarrange(pinfv, pinfh, pinfn, ncol = 3, labels = c("A", "B", "C")), ggarrange(pbeta1, pbeta2, pbeta3, ncol = 3, labels = c("D", "E", "F")), nrow = 2, common.legend = TRUE, align = "v") 
