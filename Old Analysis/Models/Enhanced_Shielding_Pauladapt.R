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
  beta1_2 <- (0.5*(1/(GenTime(3.3,2.8)))) * (1-scaling)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.5*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.5*(1/(GenTime(3.3,2.8))),
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

beta4 <- function(time, tstart1) {
  ifelse((time >= tstart1 & time <= 730), #Phase 2
         0, 
         1.7*(1/(GenTime(3.3,2.8)))
  )
}

plot(beta4(seq(0,730), 71))

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - (beta1(time,tstart1,tdur,scaling) + (beta1(time,tstart1,tdur,scaling)*bv*(1-delta)))*Iv*Sv - 
      beta1(time,tstart1,tdur,scaling)*Ih*Sv -
      beta4(time,tstart1)*delta*Ir*Sv  + zeta*Rv
    
    dSh = - beta2(time,tstart1,tdur,scaling)*Ih*Sh - 
      beta2(time,tstart1,tdur,scaling)*Ir*Sh -
      beta1(time,tstart1,tdur,scaling)*Iv*Sh + zeta*Rh
    
    dSr = - (beta3(time,tstart1,tdur,scaling) + (beta3(time,tstart1,tdur,scaling)*bn*(1-delta)))*Ir*Sr - 
      beta2(time,tstart1,tdur,scaling)*Ih*Sr -
      beta4(time,tstart1)*delta*Iv*Sr + zeta*Rr
    
    
    dIv = (beta1(time,tstart1,tdur,scaling) + (beta1(time,tstart1,tdur,scaling)*bv*(1-delta)))*Iv*Sv + 
      beta1(time,tstart1,tdur,scaling)*Ih*Sv +
      beta4(time,tstart1)*delta*Ir*Sv - gamma*Iv
    
    dIh = beta2(time,tstart1,tdur,scaling)*Ih*Sh + 
      beta2(time,tstart1,tdur,scaling)*Ir*Sh +
      beta1(time,tstart1,tdur,scaling)*Iv*Sh - gamma*Ih
    
    dIr = (beta3(time,tstart1,tdur,scaling) + (beta3(time,tstart1,tdur,scaling)*bn*(1-delta)))*Ir*Sr  + 
      beta2(time,tstart1,tdur,scaling)*Ih*Sr +
      beta4(time,tstart1)*delta*Iv*Sr - gamma*Ir
    
    
    dRv = gamma*Iv - zeta*Rv
    
    dRh = gamma*Ih - zeta*Rh
    
    dRr = gamma*Ir - zeta*Rr
    return(list(c(dSv, dSh, dSr, dIv, dIh, dIr, dRv, dRh, dRr)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####

#Initial Conditions and Times

init <- c(Sv = 0.2 - (0.0001*0.2), Sh = 0.2 - (0.0001*0.2), Sr = 0.6 - (0.0001*0.6), 
          Iv = 0.0001*0.2, Ih = 0.0001*0.2, Ir = 0.0001*0.6, 
          Rv= 0, Rh = 0, Rr = 0)

times <- seq(0,730,by = 1)

parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5,
          delta = 0.5,
          bv = 0.6/0.2,
          bn = 0.2/0.6) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))

out1$Sv <- out1$Sv/0.20; out1$Sh <- out1$Sh/0.2; out1$Sr <- out1$Sr/0.60
out1$Iv <- out1$Iv/0.20; out1$Ih <- out1$Ih/0.2; out1$Ir <- out1$Ir/0.60
out1$Rv <- out1$Rv/0.20; out1$Rh <- out1$Rh/0.2; out1$Rr <- out1$Rr/0.60

out1$Beta_vv <- (beta1(seq(0,730),71,(6*7), 0.5) + beta1(seq(0,730),71,(6*7), 0.5)*as.numeric(parms[7])*(1-as.numeric(parms[6])))
out1$Beta_vh_hv <- beta1(seq(0,730), 71, (6*7), 0.5)
out1$Beta_hh_rh_hr <- beta2(seq(0,730), 71, (6*7), 0.5)
out1$Beta_rr <- (beta3(seq(0,730),71,(6*7), 0.5) + beta3(seq(0,730),71,(6*7), 0.5)*as.numeric(parms[8])*(1-as.numeric(parms[6])))
out1$Beta_rv_vr <- beta4(seq(0,730), 71)*as.numeric(parms[6])

colnames(out1) <- c("Time", "Suscv", "Susch", "Suscnv", "Infected_Iv", "Infected_Ih", "Infected_Ir", 
                    "Recovv", "Recovh", "Recovnv", "Beta_vv", "Beta_vh_hv", "Beta_hh_rh_hr", "Beta_rr", "Beta_rv_vr")

statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Ih", "Infected_Ir"))
statsbeta1 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_vv"))
statsbeta2 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_vh_hv"))
statsbeta3 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_hh_rh_hr"))
statsbeta4 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_rr"))
statsbeta5 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_rv_vr"))


#### Aggregated Plots ####
ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.1),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0))

pbeta1 <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta_vv") + scale_y_continuous(limits = c(0,0.5),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pbeta2 <- ggplot(data = statsbeta2, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta1_vh/hv") + scale_y_continuous(limits = c(0,0.5),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pbeta3 <- ggplot(data = statsbeta3, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta_hh/hr/rh") + scale_y_continuous(limits = c(0,0.5),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pbeta4 <- ggplot(data = statsbeta4, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta_rr") + scale_y_continuous(limits = c(0,0.5),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pbeta5 <- ggplot(data = statsbeta5, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta_rv/vr") + scale_y_continuous(limits = c(0,0.5),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

ggarrange(pbeta1, pbeta2, pbeta3, pbeta4, pbeta5, nrow = 5, ncol = 1, heights = c(0.2, 0.2, 0.2, 0.2, 0.2)) 
