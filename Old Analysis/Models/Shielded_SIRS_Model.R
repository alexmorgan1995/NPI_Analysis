rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot");library("reshape2"); library("dplyr"); library("RColorBrewer"); library("ggpubr")

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur1, tdur2, tdur3, beta1_2) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         (1.7*(1/(GenTime(3.3,2.8))))*0.353,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+tdur2), #Phase 2
                beta1_2,
                ifelse((time >= tstart1+tdur1+tdur2 & time <= tstart1+tdur1+tdur2+tdur3), #Phase 3
                       (1.7*(1/(GenTime(3.3,2.8)))),
                       ifelse((time >= tstart1+tdur1+tdur2+tdur3 & time <= 730), #Phase 3
                              (1.7*(1/(GenTime(3.3,2.8)))),
                              (1.7*(1/(GenTime(3.3,2.8))))
                              )
                       )
                )
         )
}

plot(beta1(seq(0,730), 100, (6*7), (24*7), (0*7), 0.07))

beta2 <- function(time, tstart1, tdur1, tdur2, tdur3) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         (1.7*(1/(GenTime(3.3,2.8))))*0.353,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+tdur2), #Phase 2
                (1.7*(1/(GenTime(3.3,2.8)))),
                ifelse((time >= tstart1+tdur1+tdur2 & time <= tstart1+tdur1+tdur2+tdur3), #Phase 3
                       (1.7*(1/(GenTime(3.3,2.8)))),
                       ifelse((time >= tstart1+tdur1+tdur2+tdur3 & time <= 730), #Phase 4 (after)
                              (1.7*(1/(GenTime(3.3,2.8)))),
                              (1.7*(1/(GenTime(3.3,2.8))))
                       )
                )
         )
  )
}

plot(beta2(seq(0,730), 100, (6*7), (24*7), (0*7)))

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Iv*Sv - beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Inv*Sv + zeta*Rv
    
    dSnv = - beta2(time,tstart1,tdur1,tdur2,tdur3)*Inv*Snv - beta2(time,tstart1,tdur1,tdur2,tdur3)*Iv*Snv + zeta*Rnv
    
    dIv = beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Iv*Sv + beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Inv*Sv - gamma*Iv
    
    dInv =  beta2(time,tstart1,tdur1,tdur2,tdur3)*Inv*Snv + beta2(time,tstart1,tdur1,tdur2,tdur3)*Iv*Snv - gamma*Inv
    
    dRv = gamma*Iv - zeta*Rv
    
    dRnv = gamma*Inv - zeta*Rnv

    return(list(c(dSv, dSnv, dIv, dInv, dRv, dRnv)))
  })
}

#### Optimising Phase 2 + 3 Duration ####

#Initial Conditions and Times

init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, Rv= 0, Rnv = 0)
times <- seq(0,730,by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur1 = 6*7,
          tdur2 = 24*7,
          tdur3 = 0,
          beta1_2 = 0.07)

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
out1$Iv <- out1$Iv/0.20
out1$Inv <- out1$Inv/0.80
out1$Rv <- out1$Rv/0.20
out1$Rnv <- out1$Rnv/0.80
out1$Beta1 <- beta1(seq(0,730), 71, (6*7), (24*7), (0*7), 0.07)
out1$Beta2 <- beta2(seq(0,730), 71, (6*7), (24*7), (0*7))
  
ggplot(out1, aes(x= time, y = Iv)) + geom_line()

colnames(out1) <- c("Time", "Suscv", "Suscnv","Infected_Iv", "Infected_Inv", "Recovv", "Recovnv", "Beta1", "Beta2")

statsinfec <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Inv"))
statsrecov <- melt(out1, id.vars = c("Time"), measure.vars = c("Recovv", "Recovnv"))
statsbeta1 <- melt(out1, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbeta1$comp <- "Vulnerable"
statsbeta2 <- melt(out1, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2")); statsbeta2$comp <- "non-Vulnerable"
statsbetacomb <- rbind.data.frame(statsbeta1, statsbeta2)

#### Aggregated Plots

pinfv <- ggplot(data = statsinfec, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.05) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))  + scale_color_manual(values=c("darkred","darkblue"))

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))  + scale_color_manual(values=c("darkred","darkblue"))

pbeta <- ggplot(data = statsbetacomb, aes(x = (Time), y = value, lty = variable, col = comp)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + scale_linetype_manual(values = c(4,1)) +  
  scale_color_manual(values=c("darkblue", "darkred"))

ggarrange(pinfv,
          prec,
          pbeta,
          nrow = 3,
          heights = c(0.3,0.3, 0.3))

#### Checking the Proportion Infected ####

vulfracdata <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
vulfrac <- data.frame(matrix(nrow = 0, ncol = 7))

for(i in 1:length(vulfracdata)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 7))
  init <- c(Sv = vulfracdata[i], Snv = (1-vulfracdata[i])-0.0001, Iv = 0, Inv = 0.0001, Rv= 0, Rnv = 0)
  times <- seq(0,730,by = 1)
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            zeta = 1/365,
            tstart1 = 71, 
            tdur1 = 6*7,
            tdur2 = 24*7,
            tdur3 = 0,
            beta1_2 = 0.07)
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$Iv <- out1$Iv/vulfracdata[i]
  out1$Inv <- out1$Inv/(1-vulfracdata[i])
  temp[,1] <- as.character(vulfracdata[i])
  temp[,2] <- out1$time
  temp[,3] <- out1$Iv
  temp[,4] <- out1$Inv
  temp[,5] <- out1$Rv + out1$Rnv
  temp[,6] <- beta1(seq(0,730), 71, (6*7), (24*7), (0*7), 0.07)
  temp[,7] <- beta2(seq(0,730), 71, (6*7), (24*7), (0*7))
  vulfrac <- rbind.data.frame(vulfrac, temp)
  print((i/length(vulfracdata))*100)
}

colnames(vulfrac) <- c("Vulnerable_Fraction", "Time", "Infected_Iv", "Infected_Inv", "Recov", "Beta1", "Beta2")

statsinfecv <- melt(vulfrac, id.vars = c("Vulnerable_Fraction", "Time"), measure.vars = c("Infected_Iv"))
statsinfecnv <- melt(vulfrac, id.vars = c("Vulnerable_Fraction", "Time"), measure.vars = c("Infected_Inv"))
statsrecov <- melt(vulfrac, id.vars = c("Vulnerable_Fraction", "Time"), measure.vars = c("Recov"))

statsbeta1 <- melt(vulfrac, id.vars =  c("Vulnerable_Fraction", "Time"), measure.vars = c("Beta1", "Beta2"))
statsbeta2 <- melt(vulfrac, id.vars =  c("Vulnerable_Fraction", "Time"), measure.vars = c("Beta1", "Beta2"))
statsbetacomb <- rbind.data.frame(statsbeta1, statsbeta2)

#### Seperated Plots ####

pinfv <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Vulnerable Infected") + scale_y_continuous(limits = c(0,0.05) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=11),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette="Paired")

pinfnv <- ggplot(data = statsinfecnv, aes(x = (Time), y = value, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=11),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette="Paired")

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Total Recovered") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

pbeta <- ggplot(data = statsbetacomb, aes(x = (Time), y = value, lty = variable, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + scale_linetype_manual(values = c(4,1)) +
  scale_color_brewer(palette="Paired")

ggarrange(pinfv,
          pinfnv,
          prec,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))

#### Aggregated Plots ####

pinfv <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Vulnerable Infected") + scale_y_continuous(limits = c(0,0.05) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=11),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette="Paired")

pinfnv <- ggplot(data = statsinfecnv, aes(x = (Time), y = value, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=11),axis.title.x= element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette="Paired")

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Total Recovered") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

pbeta <- ggplot(data = statsbetacomb, aes(x = (Time), y = value, lty = variable, col = Vulnerable_Fraction)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  + scale_linetype_manual(values = c(4,1)) +
  scale_color_brewer(palette="Paired")

ggarrange(pinfv,
          pinfnv,
          prec,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))
