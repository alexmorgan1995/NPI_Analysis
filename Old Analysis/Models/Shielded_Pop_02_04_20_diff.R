rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot");library("reshape2"); library("dplyr"); library("RColorBrewer"); library("ggpubr")

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         0.097,
         ifelse((time >= int_timeend & time <= int_timeend + ((30-6)*7)),
                0,
                ifelse((time >= int_timeend + (24*7) & time <= 730),
                       (1.5*(1/(GenTime(4.6,2.4)))), #after int
                       (1.5*(1/(GenTime(4.6,2.4)))) #before int
                )
         )
  )
}

plot(beta1(seq(0,730), 100, 100 + (6*7) ))

beta2 <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         0.097,
         ifelse((time >= int_timeend & time <= int_timeend + ((30-6)*7)),
                0.161,
                ifelse((time >= int_timeend + (24*7) & time <= 730),
                       (1.5*(1/(GenTime(4.6,2.4)))), #after int
                       (1.5*(1/(GenTime(4.6,2.4)))) #before int
                )
         )
  )
}

plot(beta2(seq(0,730), 100, 100 + (6*7) ))

#### Different Beta ####
#Function for Shielded/non-Shielded Pop
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta1(time, int_timestart, int_timeend)*Iv*Sv - beta1(time, int_timestart, int_timeend)*Inv*Sv 
    
    dSnv = -  beta2(time, int_timestart, int_timeend)*Inv*Snv -  beta2(time, int_timestart, int_timeend)*Iv*Snv 
    
    dIv = beta1(time, int_timestart, int_timeend)*Iv*Sv + beta1(time, int_timestart, int_timeend)*Inv*Sv - gamma*Iv
    
    dInv =  beta2(time, int_timestart, int_timeend)*Inv*Snv +  beta2(time, int_timestart, int_timeend)*Iv*Snv - gamma*Inv
    
    dR = gamma*Iv + gamma*Inv
    
    dCv = beta1(time, int_timestart, int_timeend)*Iv*Sv + beta1(time, int_timestart, int_timeend)*Inv*Sv
    
    dCnv =  beta2(time, int_timestart, int_timeend)*Inv*Snv +  beta2(time, int_timestart, int_timeend)*Iv*Snv 
    return(list(c(dSv, dSnv, dIv, dInv, dR, dCv, dCnv)))
  })
}

#Initial Conditions and Times
init <- c(Sv = 0.15, Snv = 0.85-0.0001, Iv = 0, Inv = 0.0001, R = 0, Cv= 0, Cnv = 0)
times <- seq(0,730,by = 1)

parms = c(gamma = 1/(GenTime(4.6,2.4)), 
          int_timestart = 100, 
          int_timeend = 100 + 6*7)

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))

out1$Iv <- out1$Iv/0.15
out1$Inv <- out1$Inv/0.85

out1$beta1 <- beta1(seq(0,730), 100, 100 + (6*7) )
out1$beta2 <- beta2(seq(0,730), 100, 100 + (6*7) )

colnames(out1) <- c("Time", "Suscv", "Suscnv","Infected_Is", "Infected_Ip", "Recov", "CumV", "Cumnv", "Beta1", "Beta2")

statsinfec <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Is", "Infected_Ip"))
statsrecov <- melt(out1, id.vars = c("Time"), measure.vars = c("Recov"))
statsbeta <- melt(out1, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2"))

#### Aggregated Plots

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Is",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Ip",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03) ,  expand = c(0,0)) +
  theme(legend.position = "none",legend.text=element_text(size=14), legend.title = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=11), axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

prec <- ggplot(data = statsrecov, aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkgreen") +
  labs(x ="Time (Days)", y = "Proportion Recovered") + scale_y_continuous(limits = c(0,0.65) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14), axis.text=element_text(size=14), axis.title.y=element_text(size=14),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))  

pbeta <- ggplot(data = statsbeta, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +  scale_color_manual(values=c("blue", "black"))

ggarrange(pinfv, 
          pinfnv,
          prec,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))

