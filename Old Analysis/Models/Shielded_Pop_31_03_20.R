rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot");library("reshape2"); library("dplyr"); library("RColorBrewer"); library("ggpubr")

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
betastatdecrease <- function(time, int_timestart, int_timeend, scale) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*scale,
         (1.5*(1/(GenTime(4.6,2.4)))))
}

#### Different Beta ####
#Function for Shielded/non-Shielded Pop
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSs = - betastatdecrease(time, int_timestart, int_timeend,scale)*Is*Ss - 
      betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss 
    dSp = - betaPP*Ip*Sp - betaPS*Is*Sp 
    dIs = betastatdecrease(time, int_timestart, int_timeend,scale)*Is*Ss + 
      betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss - gamma*Is
    dIp = betaPP*Ip*Sp + betaPS*Is*Sp - gamma*Ip
    dR = gamma*Is + gamma*Ip
    dCv = betastatdecrease(time, int_timestart, int_timeend,scale)*Is*Ss + 
      betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss
    dCnv = betaPP*Ip*Sp + betaPS*Is*Sp 
    return(list(c(dSs, dSp, dIs, dIp, dR, dCv, dCnv)))
  })
}

#Initial Conditions and Times
init <- c(Ss = 0.15, Sp = 0.85-0.0001, Is = 0, Ip = 0.0001, R = 0, Cv= 0, Cnv = 0)
times <- seq(0,730,by = 1)

beta2scale <- c(1, 0.8, 0.6, 0.4, 0.2, 0)

#Run the ODE for 3 and 12 Week Interventions
stats1 <- data.frame(matrix(nrow = 0, ncol = 9))

for(i in 1:length(beta2scale)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 9))
  parms = c(betaPP = 1.5*(1/(GenTime(4.6,2.4))),
            betaPS = 1.5*(1/(GenTime(4.6,2.4))),
            gamma = 1/(GenTime(4.6,2.4)), 
            int_timestart = 100, 
            int_timeend = 100 + 24*7,
            scale = beta2scale[i])
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- paste0(as.character(beta2scale[i]))
  temp[,2] <- out1$time 
  temp[,3] <- (out1$Is)/0.15
  temp[,4] <- (out1$Ip)/0.85
  temp[,5] <- out1$R 
  temp[,6] <- as.numeric(parms[1])
  temp[,7] <- betastatdecrease(times, 
                               as.numeric(parms[4]), 
                               as.numeric(parms[5]),
                               beta2scale[i])
  temp[,8] <- out1$Cv /0.15
  temp[,9] <- out1$Cnv /0.85
  stats1 <- rbind(stats1, temp)
}

colnames(stats1) <- c("Beta2_Scaling", "Time", "Infected_Is", "Infected_Ip", "Recov", "Beta1", "Beta2", "CumV", "Cumnv")

#IV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Is[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#INV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Ip[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#CUM
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][731,]
  print(t)
}


statsinfec <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Infected_Is", "Infected_Ip"))
statsrecov <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Recov"))
statsbeta1 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta1"))
statsbeta2 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta2"))

statsinfec$Beta2_Scaling <- factor(statsinfec$Beta2_Scaling, levels = unique(statsinfec$Beta2_Scaling))
statsrecov$Beta2_Scaling <- factor(statsrecov$Beta2_Scaling, levels = unique(statsrecov$Beta2_Scaling))
statsbeta1$Beta2_Scaling <- factor(statsbeta1$Beta2_Scaling, levels = unique(statsbeta1$Beta2_Scaling))
statsbeta2$Beta2_Scaling <- factor(statsbeta2$Beta2_Scaling, levels = unique(statsbeta2$Beta2_Scaling))

#### Aggregated Plots

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Is",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.075) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Ip",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.075) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = Beta2_Scaling, factor = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,0.65) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))  + scale_color_brewer(palette="Paired")

pbeta1 <- ggplot(data = statsbeta2, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))   + scale_color_brewer(palette="Paired")

pbeta2 <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", axis.text=element_text(size=13), axis.title.y=element_text(size=14), 
        axis.title.x=element_text(size=14), legend.title = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + 
  scale_color_brewer(palette="Paired", name = "??1 Relative to ??2") 

ggarrange(ggarrange(pinfv, pinfnv, ncol = 2),
          prec,
          pbeta1,
          pbeta2,
          nrow = 4,
          heights = c(0.3,0.3, 0.2, 0.3))

#### No SS Modification ####

#Function for Shielded/non-Shielded Pop
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSs = - betaSS*Is*Ss - 
      betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss 
    dSp = - betaPP*Ip*Sp - betaPS*Is*Sp 
    dIs = betaSS*Is*Ss + betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss - gamma*Is
    dIp = betaPP*Ip*Sp + betaPS*Is*Sp - gamma*Ip
    dR = gamma*Is + gamma*Ip
    dCv = betaSS*Is*Ss + betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss
    dCnv = betaPP*Ip*Sp + betaPS*Is*Sp 
    return(list(c(dSs, dSp, dIs, dIp, dR, dCv, dCnv)))
  })
}

#Initial Conditions and Times
init <- c(Ss = 0.15, Sp = 0.85-0.0001, Is = 0, Ip = 0.0001, R = 0, Cv = 0, Cnv = 0)
times <- seq(0,730,by = 1)

beta2scale <- c(1, 0.8, 0.6, 0.4, 0.2, 0)

#Run the ODE for 3 and 12 Week Interventions
stats1 <- data.frame(matrix(nrow = 0, ncol = 9))

for(i in 1:length(beta2scale)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 9))
  parms = c(betaSS = 1.5*(1/(GenTime(4.6,2.4))),
            betaPP = 1.5*(1/(GenTime(4.6,2.4))),
            betaPS = 1.5*(1/(GenTime(4.6,2.4))),
            gamma = 1/(GenTime(4.6,2.4)), 
            int_timestart = 100, 
            int_timeend = 100 + 24*7,
            scale = beta2scale[i])
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- paste0(as.character(beta2scale[i]))
  temp[,2] <- out1$time 
  temp[,3] <- (out1$Is)/0.15
  temp[,4] <- (out1$Ip)/0.85
  temp[,5] <- out1$R 
  temp[,6] <- as.numeric(parms[1])
  temp[,7] <- betastatdecrease(times, 
                               as.numeric(parms[5]), 
                               as.numeric(parms[6]),
                               beta2scale[i])
  temp[,8] <- out1$Cv /0.15
  temp[,9] <- out1$Cnv /0.85
  stats1 <- rbind(stats1, temp)
}

colnames(stats1) <- c("Beta2_Scaling", "Time", "Infected_Is", "Infected_Ip", "Recov", "Beta1", "Beta2", "Cumv", "Cumnv")

#IV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Is[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#INV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Ip[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#CUM
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][731,]
  print(t)
}

statsinfec <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Infected_Is", "Infected_Ip"))
statsrecov <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Recov"))
statsbeta1 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta1"))
statsbeta2 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta2"))

statsinfec$Beta2_Scaling <- factor(statsinfec$Beta2_Scaling, levels = unique(statsinfec$Beta2_Scaling))
statsrecov$Beta2_Scaling <- factor(statsrecov$Beta2_Scaling, levels = unique(statsrecov$Beta2_Scaling))
statsbeta1$Beta2_Scaling <- factor(statsbeta1$Beta2_Scaling, levels = unique(statsbeta1$Beta2_Scaling))
statsbeta2$Beta2_Scaling <- factor(statsbeta2$Beta2_Scaling, levels = unique(statsbeta2$Beta2_Scaling))

#### Aggregated Plots
pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Is",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.075) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Ip",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.075) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = Beta2_Scaling, factor = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,0.65) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))  + scale_color_brewer(palette="Paired")

pbeta1 <- ggplot(data = statsbeta2, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))   + scale_color_brewer(palette="Paired")

pbeta2 <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", axis.text=element_text(size=13), axis.title.y=element_text(size=14), 
        axis.title.x=element_text(size=14), legend.title = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + 
  scale_color_brewer(palette="Paired", name = "??1 Relative to ??2") 

ggarrange(ggarrange(pinfv, pinfnv, ncol = 2),
          prec,
          pbeta1,
          pbeta2,
          nrow = 4,
          heights = c(0.3,0.3, 0.2, 0.3))


#### R0 kept at 2.4 ####

betastatdecrease1 <- function(time, int_timestart, int_timeend, scale) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (2.4*(1/(GenTime(4.6,2.4))))*scale,
         (2.4*(1/(GenTime(4.6,2.4)))))
}

#Function for Shielded/non-Shielded Pop
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSs = - betastatdecrease1(time, int_timestart, int_timeend,scale)*Is*Ss - 
      betastatdecrease1(time, int_timestart, int_timeend,scale)*Ip*Ss 
    dSp = - betaPP*Ip*Sp - betaPS*Is*Sp 
    dIs = betastatdecrease1(time, int_timestart, int_timeend,scale)*Is*Ss + 
      betastatdecrease1(time, int_timestart, int_timeend,scale)*Ip*Ss - gamma*Is
    dIp = betaPP*Ip*Sp + betaPS*Is*Sp - gamma*Ip
    dR = gamma*Is + gamma*Ip
    dCv = betastatdecrease1(time, int_timestart, int_timeend,scale)*Is*Ss + 
      betastatdecrease1(time, int_timestart, int_timeend,scale)*Ip*Ss
    dCnv = betaPP*Ip*Sp + betaPS*Is*Sp 
    return(list(c(dSs, dSp, dIs, dIp, dR, dCv, dCnv)))
  })
}

#Initial Conditions and Times
init <- c(Ss = 0.15, Sp = 0.85-0.0001, Is = 0, Ip = 0.0001, R = 0, Cv = 0, Cnv = 0)
times <- seq(0,730,by = 1)

beta2scale <- c(1, 0.8, 0.6, 0.4, 0.2, 0)

parms = c(betaPP = 2.4*(1/(GenTime(4.6,2.4))),
          betaPS = 2.4*(1/(GenTime(4.6,2.4))),
          gamma = 1/(GenTime(4.6,2.4)), 
          int_timestart = 100, 
          int_timeend = 100 + 24*7,
          scale = 1)
out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))

out1$Is <- out1$Is / 0.15
out1$Ip <- out1$Ip / 0.85

#Run the ODE for 3 and 12 Week Interventions
stats1 <- data.frame(matrix(nrow = 0, ncol = 9))

for(i in 1:length(beta2scale)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 9))
  parms = c(betaPP = 2.4*(1/(GenTime(4.6,2.4))),
            betaPS = 2.4*(1/(GenTime(4.6,2.4))),
            gamma = 1/(GenTime(4.6,2.4)), 
            int_timestart = 35, 
            int_timeend = 35 + 24*7,
            scale = beta2scale[i])
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- paste0(as.character(beta2scale[i]))
  temp[,2] <- out1$time 
  temp[,3] <- (out1$Is)/0.15
  temp[,4] <- (out1$Ip)/0.85
  temp[,5] <- out1$R 
  temp[,6] <- as.numeric(parms[1])
  temp[,7] <- betastatdecrease1(times, 
                               as.numeric(parms[4]), 
                               as.numeric(parms[5]),
                               beta2scale[i])
  temp[,8] <- out1$Cv /0.15
  temp[,9] <- out1$Cnv /0.85
  stats1 <- rbind(stats1, temp)
}

colnames(stats1) <- c("Beta2_Scaling", "Time", "Infected_Is", "Infected_Ip", "Recov", "Beta1", "Beta2", "Cumv", "Cumnv")

#IV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Is[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#INV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Ip[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#CUM
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][731,]
  print(t)
}

statsinfec <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Infected_Is", "Infected_Ip"))
statsrecov <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Recov"))
statsbeta1 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta1"))
statsbeta2 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta2"))

statsinfec$Beta2_Scaling <- factor(statsinfec$Beta2_Scaling, levels = unique(statsinfec$Beta2_Scaling))
statsrecov$Beta2_Scaling <- factor(statsrecov$Beta2_Scaling, levels = unique(statsrecov$Beta2_Scaling))
statsbeta1$Beta2_Scaling <- factor(statsbeta1$Beta2_Scaling, levels = unique(statsbeta1$Beta2_Scaling))
statsbeta2$Beta2_Scaling <- factor(statsbeta2$Beta2_Scaling, levels = unique(statsbeta2$Beta2_Scaling))

#### Aggregated Plots
pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Is",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.25) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Ip",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.25) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = Beta2_Scaling, factor = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))  + scale_color_brewer(palette="Paired")

pbeta1 <- ggplot(data = statsbeta2, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))   + scale_color_brewer(palette="Paired")

pbeta2 <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", axis.text=element_text(size=13), axis.title.y=element_text(size=14), 
        axis.title.x=element_text(size=14), legend.title = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + 
  scale_color_brewer(palette="Paired", name = "??1 Relative to ??2") 

ggarrange(ggarrange(pinfv, pinfnv, ncol = 2),
          prec,
          pbeta1,
          pbeta2,
          nrow = 4,
          heights = c(0.3,0.3, 0.2, 0.3))


#### Change the Trigger time ####

#Function for Shielded/non-Shielded Pop
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSs = - betastatdecrease(time, int_timestart, int_timeend,scale)*Is*Ss - 
      betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss 
    dSp = - betaPP*Ip*Sp - betaPS*Is*Sp 
    dIs = betastatdecrease(time, int_timestart, int_timeend,scale)*Is*Ss + 
      betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss - gamma*Is
    dIp = betaPP*Ip*Sp + betaPS*Is*Sp - gamma*Ip
    dR = gamma*Is + gamma*Ip
    dCv = betastatdecrease(time, int_timestart, int_timeend,scale)*Is*Ss + 
      betastatdecrease(time, int_timestart, int_timeend,scale)*Ip*Ss
    dCnv = betaPP*Ip*Sp + betaPS*Is*Sp 
    return(list(c(dSs, dSp, dIs, dIp, dR, dCv, dCnv)))
  })
}

#Initial Conditions and Times
init <- c(Ss = 0.15, Sp = 0.85-0.0001, Is = 0, Ip = 0.0001, R = 0, Cv= 0, Cnv = 0)
times <- seq(0,730,by = 1)

beta2scale <- c(1, 0.8, 0.6, 0.4, 0.2, 0)

#Run the ODE for 3 and 12 Week Interventions
stats1 <- data.frame(matrix(nrow = 0, ncol = 9))

for(i in 1:length(beta2scale)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 9))
  parms = c(betaPP = 1.5*(1/(GenTime(4.6,2.4))),
            betaPS = 1.5*(1/(GenTime(4.6,2.4))),
            gamma = 1/(GenTime(4.6,2.4)), 
            int_timestart = 125, 
            int_timeend = 125 + 24*7,
            scale = beta2scale[i])
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- paste0(as.character(beta2scale[i]))
  temp[,2] <- out1$time 
  temp[,3] <- (out1$Is)/0.15
  temp[,4] <- (out1$Ip)/0.85
  temp[,5] <- out1$R 
  temp[,6] <- as.numeric(parms[1])
  temp[,7] <- betastatdecrease(times, 
                               as.numeric(parms[4]), 
                               as.numeric(parms[5]),
                               beta2scale[i])
  temp[,8] <- out1$Cv /0.15
  temp[,9] <- out1$Cnv /0.85
  stats1 <- rbind(stats1, temp)
}

colnames(stats1) <- c("Beta2_Scaling", "Time", "Infected_Is", "Infected_Ip", "Recov", "Beta1", "Beta2", "CumV", "Cumnv")

#IV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Is[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#INV
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][which.max(stats1$Infected_Ip[stats1$Beta2_Scaling == as.character(beta2scale[i])]),]
  print(t)
}

#CUM
for(i in 1:length(beta2scale)){
  t <- stats1[stats1$Beta2_Scaling == as.character(beta2scale[i]),][731,]
  print(t)
}

statsinfec <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Infected_Is", "Infected_Ip"))
statsrecov <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Recov"))
statsbeta1 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta1"))
statsbeta2 <- melt(stats1, id.vars = c("Beta2_Scaling", "Time"), measure.vars = c("Beta2"))

statsinfec$Beta2_Scaling <- factor(statsinfec$Beta2_Scaling, levels = unique(statsinfec$Beta2_Scaling))
statsrecov$Beta2_Scaling <- factor(statsrecov$Beta2_Scaling, levels = unique(statsrecov$Beta2_Scaling))
statsbeta1$Beta2_Scaling <- factor(statsbeta1$Beta2_Scaling, levels = unique(statsbeta1$Beta2_Scaling))
statsbeta2$Beta2_Scaling <- factor(statsbeta2$Beta2_Scaling, levels = unique(statsbeta2$Beta2_Scaling))

#### Aggregated Plots

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Is",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.075) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Ip",], aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.075) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = Beta2_Scaling, factor = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,0.65) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))  + scale_color_brewer(palette="Paired")

pbeta1 <- ggplot(data = statsbeta2, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))   + scale_color_brewer(palette="Paired")

pbeta2 <- ggplot(data = statsbeta1, aes(x = (Time), y = value, col = Beta2_Scaling)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", axis.text=element_text(size=13), axis.title.y=element_text(size=14), 
        axis.title.x=element_text(size=14), legend.title = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + 
  scale_color_brewer(palette="Paired", name = "??1 Relative to ??2") 

ggarrange(ggarrange(pinfv, pinfnv, ncol = 2),
          prec,
          pbeta1,
          pbeta2,
          nrow = 4,
          heights = c(0.3,0.3, 0.2, 0.3))