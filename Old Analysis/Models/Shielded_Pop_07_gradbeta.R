rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot");library("reshape2"); library("dplyr"); library("RColorBrewer"); library("ggpubr")

#### Generic Functions ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur1, tdur2, beta1_2) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         0.064,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+tdur2), #Phase 2
                beta1_2,
                ifelse((time >= tstart1+tdur1+tdur2 & time <= 730), #Phase 4
                       (1.5*(1/(GenTime(4.6,2.4)))), #After Intervention
                       (1.5*(1/(GenTime(4.6,2.4))))#Before Intervention
                       )
                )
         )
  }

plot(beta1(seq(0,730), 100, (6*7), 24*7, 0.064))

beta12 <- function(time, tstart1, tdur1, tdur2, beta1_2) {
  betalin <- approxfun(x=c(tstart1+tdur1, tstart1+tdur1+(tdur2/2)),y = c(0.064, beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 2
         0.064,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+(tdur2/2)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur1+(tdur2/2) & time <= tstart1+tdur1+(tdur2)),
                       beta1_2,
                       ifelse((time >= tstart1+tdur1+tdur2 & time <= 730),
                              (1.5*(1/(GenTime(4.6,2.4)))), 
                              (1.5*(1/(GenTime(4.6,2.4))))
                       )
                )
         )
  )
}
                       
plot(beta12(seq(0,730), 100, (6*7), (24*7), 0))

beta13 <- function(time, tstart1, tdur1, tdur2, beta1_2) {
  betalin <- approxfun(x=c(tstart1+tdur1, tstart1+tdur1+tdur2),y = c(0.064, beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 2
         0.064,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+tdur2), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur1+tdur2 & time <= 730), #Phase 4
                       (1.5*(1/(GenTime(4.6,2.4)))), 
                       (1.5*(1/(GenTime(4.6,2.4)))) #Before Intervention
                )
         )
  )
}

plot(beta13(seq(0,730), 100, (6*7), (24*7), 0))


test <- data.frame("time" = seq(0, 730), "beta" = beta2(seq(0,730), 100, 6*7, 24*7))

ggplot(data = test, aes(x = (time), y = beta)) + geom_line(stat = "identity", col = "darkred", size = 1.02) +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) + scale_color_manual(values = c("darkred", "darkblue")) + 
  scale_linetype_manual(values = c(2, 1)) 

beta2 <- function(time, tstart1, tdur1, tdur2) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         0.064,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+tdur2), #Phase 2
                (1.5*(1/(GenTime(4.6,2.4)))),
                ifelse((time >= tstart1+tdur1+tdur2 & time <= 730), #Phase 3
                       (1.5*(1/(GenTime(4.6,2.4)))), #Phase 3
                       (1.5*(1/(GenTime(4.6,2.4)))) #Before Intervention
                )
         )
  )
}

plot(beta2(seq(0,730), 100, (6*7), (24*7)))

#### Baseline ####

#Function for Shielded/non-Shielded Pop
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta1(time,tstart1,tdur1,tdur2,beta1_2)*Iv*Sv - beta1(time,tstart1,tdur1,tdur2,beta1_2)*Inv*Sv 
    
    dSnv = - beta2(time,tstart1,tdur1,tdur2)*Inv*Snv - beta2(time,tstart1,tdur1,tdur2)*Iv*Snv 
    
    dIv = beta1(time,tstart1,tdur1,tdur2,beta1_2)*Iv*Sv + beta1(time,tstart1,tdur1,tdur2,beta1_2)*Inv*Sv - gamma*Iv
    
    dInv =  beta2(time,tstart1,tdur1,tdur2)*Inv*Snv + beta2(time,tstart1,tdur1,tdur2)*Iv*Snv - gamma*Inv
    
    dRv = gamma*Iv 
    
    dRnv = gamma*Inv
    return(list(c(dSv, dSnv, dIv, dInv, dRv, dRnv)))
  })
}

#Initial Conditions and Times
init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, Rv = 0, Rnv = 0)
times <- seq(0,730,by = 1)

#Beta Functions are specified in the ODE
parms = c(gamma = 1/(GenTime(4.6,2.4)), 
          tstart1 = 100, 
          tdur1 = 6*7,
          tdur2 = 24*7,
          beta1_2 = 0)

out1 <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

#Extra Data Columns in Dataframe
out1$Iv <- out1$Iv/0.20
out1$Inv <- out1$Inv/0.80
out1$Rv <- out1$Rv/0.20
out1$Rnv <- out1$Rnv/0.80
out1$beta1 <- beta1(seq(0,730), 100, (6*7), (24*7), 0)
out1$beta2 <- beta2(seq(0,730), 100, (6*7), (24*7))

#Cumulative Infections are scaled relative to the relative V and nV population sizes - 0.15 and 0.85
#You should not compare these numbers to the recovered fraction - they will not add up

colnames(out1) <- c("Time", "Suscv", "Suscnv","Infected_Iv", "Infected_Inv", "Recovv", "Recovnv", "Beta1", "Beta2")

statsinfec <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Inv"))
statsrecov <- melt(out1, id.vars = c("Time"), measure.vars = c("Recovv","Recovnv"))
statsbetav <- melt(out1, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbetav$inf <- "Vulnerable"
statsbetanv <- melt(out1, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbetanv$inf <- "non_Vulnerable" 

statsbetacomb <- rbind.data.frame(statsbetav, statsbetanv)

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Iv",], aes(x = (Time), y = value)) + geom_line(size = 1.02,  col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Inv",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.035) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,0.65) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_manual(values=c("darkred", "darkblue"))

pbeta <- ggplot(data = statsbetacomb, aes(x = (Time), y = value, col = inf, lty = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +  scale_linetype_manual(values=c(2, 1)) + scale_color_manual(values=c("darkblue", "darkred"))

ggarrange(pinfv, 
          pinfnv,
          prec,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))

#### Ramp up 12 Weeks ####

#Function for Shielded/non-Shielded Pop
SIR12 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta12(time,tstart1,tdur1,tdur2,beta1_2)*Iv*Sv - beta12(time,tstart1,tdur1,tdur2,beta1_2)*Inv*Sv 
    
    dSnv = - beta2(time,tstart1,tdur1,tdur2)*Inv*Snv - beta2(time,tstart1,tdur1,tdur2)*Iv*Snv 
    
    dIv = beta12(time,tstart1,tdur1,tdur2,beta1_2)*Iv*Sv + beta12(time,tstart1,tdur1,tdur2,beta1_2)*Inv*Sv - gamma*Iv
    
    dInv =  beta2(time,tstart1,tdur1,tdur2)*Inv*Snv + beta2(time,tstart1,tdur1,tdur2)*Iv*Snv - gamma*Inv
    
    dRv = gamma*Iv
    
    dRnv = gamma*Inv

    return(list(c(dSv, dSnv, dIv, dInv, dRv, dRnv)))
  })
}

#Initial Conditions and Times
init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, Rv = 0, Rnv = 0)
times <- seq(0,730,by = 1)

#Beta Functions are specified in the ODE
parms = c(gamma = 1/(GenTime(4.6,2.4)), 
          tstart1 = 100, 
          tdur1 = 6*7,
          tdur2 = 24*7,
          beta1_2 = 0)

out2 <- data.frame(ode(y = init, func = SIR12, times = times, parms = parms))

#Extra Data Columns in Dataframe
out2$Iv <- out2$Iv/0.20
out2$Inv <- out2$Inv/0.80
out2$Rv <- out2$Rv/0.20
out2$Rnv <- out2$Rnv/0.80
out2$beta1 <- beta12(seq(0,730), 100, (6*7), (24*7), 0)
out2$beta2 <- beta2(seq(0,730), 100, (6*7), (24*7))

#Cumulative Infections are scaled relative to the relative V and nV population sizes - 0.15 and 0.85
#You should not compare these numbers to the recovered fraction - they will not add up

colnames(out2) <- c("Time", "Suscv", "Suscnv","Infected_Iv", "Infected_Inv", "Recovv", "Recovnv", "Beta1", "Beta2")

statsinfec <- melt(out2, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Inv"))
statsrecov <- melt(out2, id.vars = c("Time"), measure.vars = c("Recovv","Recovnv"))
statsbetav <- melt(out2, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbetav$inf <- "Vulnerable"
statsbetanv <- melt(out2, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbetanv$inf <- "non_Vulnerable" 

statsbetacomb <- rbind.data.frame(statsbetav, statsbetanv)

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Iv",], aes(x = (Time), y = value)) + geom_line(size = 1.02,  col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Inv",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.035) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,0.65) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_manual(values=c("darkred", "darkblue"))

pbeta <- ggplot(data = statsbetacomb, aes(x = (Time), y = value, col = inf, lty = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +  scale_linetype_manual(values=c(2, 1)) + scale_color_manual(values=c("darkblue", "darkred"))

ggarrange(pinfv, 
          pinfnv,
          prec,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))

#### Ramp up 24 Weeks ####

#Function for Shielded/non-Shielded Pop
SIR13 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta13(time,tstart1,tdur1,tdur2,beta1_2)*Iv*Sv - beta13(time,tstart1,tdur1,tdur2,beta1_2)*Inv*Sv 
    
    dSnv = - beta2(time,tstart1,tdur1,tdur2)*Inv*Snv - beta2(time,tstart1,tdur1,tdur2)*Iv*Snv 
    
    dIv = beta13(time,tstart1,tdur1,tdur2,beta1_2)*Iv*Sv + beta13(time,tstart1,tdur1,tdur2,beta1_2)*Inv*Sv - gamma*Iv
    
    dInv =  beta2(time,tstart1,tdur1,tdur2)*Inv*Snv + beta2(time,tstart1,tdur1,tdur2)*Iv*Snv - gamma*Inv
    
    dRv = gamma*Iv
    
    dRnv = gamma*Inv
    
    return(list(c(dSv, dSnv, dIv, dInv, dRv, dRnv)))
  })
}

#Initial Conditions and Times
init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, Rv= 0, Rnv = 0)
times <- seq(0,730,by = 1)

#Beta Functions are specified in the ODE
parms = c(gamma = 1/(GenTime(4.6,2.4)), 
          tstart1 = 100, 
          tdur1 = 6*7,
          tdur2 = 24*7,
          beta1_2 = 0)

out3 <- data.frame(ode(y = init, func = SIR13, times = times, parms = parms))

#Extra Data Columns in Dataframe
out3$Iv <- out3$Iv/0.20
out3$Inv <- out3$Inv/0.80
out3$Rv <- out3$Rv/0.20
out3$Rnv <- out3$Rnv/0.80
out3$beta1 <- beta13(seq(0,730), 100, (6*7), (24*7), 0)
out3$beta2 <- beta2(seq(0,730), 100, (6*7), (24*7))

#Cumulative Infections are scaled relative to the relative V and nV population sizes - 0.15 and 0.85
#You should not compare these numbers to the recovered fraction - they will not add up

colnames(out3) <- c("Time", "Suscv", "Suscnv","Infected_Iv", "Infected_Inv", "Recovv", "Recovnv", "Beta1", "Beta2")

statsinfec <- melt(out3, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Inv"))
statsrecov <- melt(out3, id.vars = c("Time"), measure.vars = c("Recovv","Recovnv"))
statsbetav <- melt(out3, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbetav$inf <- "Vulnerable"
statsbetanv <- melt(out3, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbetanv$inf <- "non_Vulnerable" 

statsbetacomb <- rbind.data.frame(statsbetav, statsbetanv)

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Iv",], aes(x = (Time), y = value)) + geom_line(size = 1.02,  col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Inv",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.035) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,0.65) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_manual(values=c("darkred", "darkblue"))

pbeta <- ggplot(data = statsbetacomb, aes(x = (Time), y = value, col = inf, lty = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +  scale_linetype_manual(values=c(2, 1)) + scale_color_manual(values=c("darkblue", "darkred"))

ggarrange(pinfv, 
          pinfnv,
          prec,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))
