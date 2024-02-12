rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur, scaling) {
  beta1_2 <- (0.8*(1/(GenTime(3.3,2.8))))*scaling
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(1/(GenTime(3.3,2.8))),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(1/(GenTime(3.3,2.8))))))}

plot(beta1(seq(0,730), 71, (6*7), 0.5), ylim = c(0,0.5))

beta2 <- function(time, tstart1, tdur, scaling) {
  beta1_2 <- (2.8*(1/(GenTime(3.3,2.8))) - ((2.8*(1/(GenTime(3.3,2.8))) - 0.9*(1/(GenTime(3.3,2.8))))*scaling))
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(1/(GenTime(3.3,2.8))),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(1/(GenTime(3.3,2.8))))))}

plot(beta2(seq(0,730), 71, (6*7), 0.5))

beta3 <- function(time, tstart1, tdur, scaling) {
  beta1_2 <- (2.8*(1/(GenTime(3.3,2.8))) - ((2.8*(1/(GenTime(3.3,2.8))) - 1.7*(1/(GenTime(3.3,2.8))))*scaling))
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(1/(GenTime(3.3,2.8))),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(1/(GenTime(3.3,2.8))))))}

plot(beta3(seq(0,730), 71, (6*7), 0.5))


beta4 <- function(time,tstart1,tdur,scaling) {
  beta1_2 <- (0.8*(1/(GenTime(3.3,2.8)))) *scaling
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(1/(GenTime(3.3,2.8))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(1/(GenTime(3.3,2.8))),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(1/(GenTime(3.3,2.8))))))}

plot(beta4(seq(0,730), 71, (6*7), 0.5))

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSv = - beta1(time,tstart1,tdur,scaling)*Iv*Sv - 
      beta1(time,tstart1,tdur,scaling)*Is*Sv - 
      beta4(time,tstart1,tdur,scaling)*Ir1*Sv - 
      beta4(time,tstart1,tdur,scaling)*Ir2*Sv - 
      beta4(time,tstart1,tdur,scaling)*Ir3*Sv + gamma*Iv
    
    dSs = - beta1(time,tstart1,tdur,scaling)*Iv*Ss - 
      beta1(time,tstart1,tdur,scaling)*Is*Ss -
      beta2(time,tstart1,tdur,scaling)*Ir1*Ss -
      beta2(time,tstart1,tdur,scaling)*Ir2*Ss -
      beta2(time,tstart1,tdur,scaling)*Ir3*Ss + gamma*Is

    dSr1 = - beta4(time,tstart1,tdur,scaling)*Iv*Sr1 - 
      beta2(time,tstart1,tdur,scaling)*Is*Sr1 - 
      beta3(time,tstart1,tdur,scaling)*Ir1*Sr1 - 
      beta3(time,tstart1,tdur,scaling)*Ir2*Sr1 - 
      beta3(time,tstart1,tdur,scaling)*Ir3*Sr1 + gamma*Ir1
      
    dSr2 = - beta4(time,tstart1,tdur,scaling)*Iv*Sr2 -
      beta2(time,tstart1,tdur,scaling)*Is*Sr2 -
      beta3(time,tstart1,tdur,scaling)*Ir1*Sr2 -
      beta3(time,tstart1,tdur,scaling)*Ir2*Sr2 -
      beta3(time,tstart1,tdur,scaling)*Ir3*Sr2 + gamma*Ir2
      
    dSr3 = - beta4(time,tstart1,tdur,scaling)*Iv*Sr3 -
      beta2(time,tstart1,tdur,scaling)*Is*Sr3 -
      beta3(time,tstart1,tdur,scaling)*Ir1*Sr3 -
      beta3(time,tstart1,tdur,scaling)*Ir2*Sr3 -
      beta3(time,tstart1,tdur,scaling)*Ir3*Sr3 + gamma*Ir3
    
    
    
    dIv = beta1(time,tstart1,tdur,scaling)*Iv*Sv + 
      beta1(time,tstart1,tdur,scaling)*Is*Sv + 
      beta4(time,tstart1,tdur,scaling)*Ir1*Sv + 
      beta4(time,tstart1,tdur,scaling)*Ir2*Sv + 
      beta4(time,tstart1,tdur,scaling)*Ir3*Sv - gamma*Iv
    
    dIs = beta1(time,tstart1,tdur,scaling)*Iv*Ss + 
      beta1(time,tstart1,tdur,scaling)*Is*Ss +
      beta2(time,tstart1,tdur,scaling)*Ir1*Ss +
      beta2(time,tstart1,tdur,scaling)*Ir2*Ss +
      beta2(time,tstart1,tdur,scaling)*Ir3*Ss - gamma*Is
    
    dIr1 = beta4(time,tstart1,tdur,scaling)*Iv*Sr1 +
      beta2(time,tstart1,tdur,scaling)*Is*Sr1 +
      beta3(time,tstart1,tdur,scaling)*Ir1*Sr1 + 
      beta3(time,tstart1,tdur,scaling)*Ir2*Sr1 + 
      beta3(time,tstart1,tdur,scaling)*Ir3*Sr1 - gamma*Ir1
    
    dIr2 = beta4(time,tstart1,tdur,scaling)*Iv*Sr2 +
      beta2(time,tstart1,tdur,scaling)*Is*Sr2 +
      beta3(time,tstart1,tdur,scaling)*Ir1*Sr2 +
      beta3(time,tstart1,tdur,scaling)*Ir2*Sr2 +
      beta3(time,tstart1,tdur,scaling)*Ir3*Sr2 - gamma*Ir2
    
    dIr3 = beta4(time,tstart1,tdur,scaling)*Iv*Sr3 +
      beta2(time,tstart1,tdur,scaling)*Is*Sr3 +
      beta3(time,tstart1,tdur,scaling)*Ir1*Sr3 +
      beta3(time,tstart1,tdur,scaling)*Ir2*Sr3 +
      beta3(time,tstart1,tdur,scaling)*Ir3*Sr3 - gamma*Ir3
    
    
    return(list(c(dSv, dSs, dSr1, dSr2, dSr3,
                  dIv, dIs, dIr1, dIr2, dIr3)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####

#Initial Conditions and Times

init <- c(Sv = 0.2 - 0.0001*0.2, Ss = 0.2 - 0.0001*0.2, 
          Sr1 = 0.2 - 0.0001*0.2, Sr2 = 0.2 - 0.0001*0.2, Sr3 = 0.2 - 0.0001*0.2,
          
          Iv = 0.0001*0.2, Is = 0.0001*0.2, Ir1 = 0.0001*0.2, Ir2 = 0.0001*0.2, Ir3 = 0.0001*0.2)

times <- seq(0, 478, by = 1)

parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))

out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3

out1$Sv <- out1$Sv/0.20; out1$Ss <- out1$Ss/0.20; out1$RemS <- out1$RemS/0.60
out1$Iv <- out1$Iv/0.20; out1$Is <- out1$Is/0.20; out1$RemI <- out1$RemI/0.60

out1$Beta1 <- beta1(times, 71, (6*7), as.numeric(parms[5])) 
out1$Beta2 <- beta2(times, 71, (6*7), as.numeric(parms[5])) 
out1$Beta3 <- beta3(times, 71, (6*7), as.numeric(parms[5]))
out1$Beta4 <- beta4(times, 71, (6*7), as.numeric(parms[5]))


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr1", "Suscr2","Suscr3", 
                    "Infected_Iv", "Infected_Is", "Infected_Ir1", "Infected_Ir2", "Infected_Ir3",
                    "RemSusc", "RemInf",
                    "Beta_1", "Beta_2", "Beta_3", "Beta_4")

statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsbeta1 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_1"))
statsbeta2 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_2"))
statsbeta3 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_3"))
statsbeta4 <- melt(out1, id.vars =  c("Time"), measure.vars = c("Beta_4"))

#### Aggregated Plots ####
ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.1), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('darkgreen','red','blue'), labels= c("Vulnerable", "Shielders", "Remainders"))

pbeta1 <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta_vv") + scale_y_continuous(limits = c(0,0.35),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pbeta2 <- ggplot(data = statsbeta2, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta1_vh/hv") + scale_y_continuous(limits = c(0,0.35),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pbeta3 <- ggplot(data = statsbeta3, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta_hh/hr/rh") + scale_y_continuous(limits = c(0,0.35),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pbeta4 <- ggplot(data = statsbeta4, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Beta_rr") + scale_y_continuous(limits = c(0,0.35),  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

ggarrange(pbeta1, pbeta2, pbeta3, pbeta4, nrow = 4, ncol = 1, heights = c(0.2, 0.2, 0.2, 0.2)) 

