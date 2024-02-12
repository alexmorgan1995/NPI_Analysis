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
                ifelse((time >= tstart1+tdur1+tdur2 & time <= 730), #Phase 3
                       (1.5*(1/(GenTime(4.6,2.4)))), #Phase 3
                       (1.5*(1/(GenTime(4.6,2.4)))) #Before Intervention
                )
         )
  )
}

plot(beta1(seq(0,730), 100, (6*7), (24*7), 0.064))


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

#### SS ####

#Function for Shielded/non-Shielded Pop
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta1(time, tstart1, tdur1, tdur2, beta1_2)*Iv*Sv - beta1(time, tstart1, tdur1, tdur2, beta1_2)*Inv*Sv 
    
    dSnv = - beta2(time,tstart1,tdur1,tdur2)*Inv*Snv - beta2(time,tstart1,tdur1,tdur2)*Iv*Snv 
    
    dIv = beta1(time, tstart1, tdur1, tdur2, beta1_2)*Iv*Sv + beta1(time, tstart1, tdur1, tdur2, beta1_2)*Inv*Sv - gamma*Iv
    
    dInv =  beta2(time,tstart1,tdur1,tdur2)*Inv*Snv + beta2(time,tstart1,tdur1,tdur2)*Iv*Snv - gamma*Inv
    
    dR = gamma*Iv + gamma*Inv
    
    dRv = gamma*Iv
    
    dRnv = gamma*Inv
    
    return(list(c(dSv, dSnv, dIv, dInv, dR, dRv, dRnv)))
  })
}

#Initial Conditions and Times
init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, R = 0, Rv= 0, Rnv = 0)
times <- seq(0,730,by = 1)

#Beta Functions are specified in the ODE
parms = c(gamma = 1/(GenTime(4.6,2.4)), 
          tstart1 = 100, 
          tdur1 = 6*7,
          tdur2 = 24*7,
          beta1_2 = 0.064)

out2 <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

#Extra Data Columns in Dataframe
out2$Iv <- out2$Iv/0.20
out2$Inv <- out2$Inv/0.80
out2$Rv <- out2$Rv/0.20
out2$Rnv <- out2$Rnv/0.80
out2$beta1 <- beta1(seq(0,730), 100, (6*7), (24*7), 0.064)
out2$beta2 <- beta2(seq(0,730), 100, (6*7), (24*7))

#Cumulative Infections are scaled relative to the relative V and nV population sizes - 0.15 and 0.85
#You should not compare these numbers to the recovered fraction - they will not add up

colnames(out2) <- c("Time", "Suscv", "Suscnv","Infected_Iv", "Infected_Inv", "Recov", "Rv", "Rnv", "Beta1", "Beta2")

statsinfec <- melt(out2, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Inv"))
statsrecov <- melt(out2, id.vars = c("Time"), measure.vars = c("Rv", "Rnv"))
statsbeta <- melt(out2, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2"))

#### Aggregated Plots

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Iv",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Inv",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.03) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0))

precov <- ggplot(data = statsrecov, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.8) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) + scale_color_manual(values = c("darkred", "darkblue")) 

pbeta <- ggplot(data = statsbeta, aes(x = (Time), y = value, lty = variable)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x = element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) + scale_linetype_manual(values = c(2, 1)) 


ggarrange(pinfv, 
          pinfnv,
          precov,
          pbeta,
          nrow = 4,
          heights = c(0.3,0.3, 0.3, 0.3))