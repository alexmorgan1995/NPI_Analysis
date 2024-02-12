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
         0.064,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+tdur2), #Phase 2
                beta1_2,
                ifelse((time >= tstart1+tdur1+tdur2 & time <= tstart1+tdur1+tdur2+tdur3), #Phase 3
                       (1.5*(1/(GenTime(4.6,2.4)))),
                       ifelse((time >= tstart1+tdur1+tdur2+tdur3 & time <= 730), #Phase 3
                              (2.4*(1/(GenTime(4.6,2.4)))),
                              (1.5*(1/(GenTime(4.6,2.4))))
                              )
                       )
                )
         )
}

plot(beta1(seq(0,730), 100, (12*7), (24*7), (24*7), 0.032))

beta2 <- function(time, tstart1, tdur1, tdur2, tdur3) {
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 1
         0.064,
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+tdur2), #Phase 2
                (1.5*(1/(GenTime(4.6,2.4)))),
                ifelse((time >= tstart1+tdur1+tdur2 & time <= tstart1+tdur1+tdur2+tdur3), #Phase 3
                       (1.5*(1/(GenTime(4.6,2.4)))),
                       ifelse((time >= tstart1+tdur1+tdur2+tdur3 & time <= 730), #Phase 3
                              (2.4*(1/(GenTime(4.6,2.4)))),
                              (1.5*(1/(GenTime(4.6,2.4))))
                       )
                )
         )
  )
}

plot(beta2(seq(0,730), 100, (12*7), (24*7), (24*7)))

#Function for Shielded/non-Shielded Pop
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSv = - beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Iv*Sv - beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Inv*Sv 
    
    dSnv = - beta2(time,tstart1,tdur1,tdur2,tdur3)*Inv*Snv - beta2(time,tstart1,tdur1,tdur2,tdur3)*Iv*Snv 
    
    dIv = beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Iv*Sv + beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Inv*Sv - gamma*Iv
    
    dInv =  beta2(time,tstart1,tdur1,tdur2,tdur3)*Inv*Snv + beta2(time,tstart1,tdur1,tdur2,tdur3)*Iv*Snv - gamma*Inv
    
    dR = gamma*Iv + gamma*Inv
    
    dCv = beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Iv*Sv + beta1(time,tstart1,tdur1,tdur2,tdur3,beta1_2)*Inv*Sv
    
    dCnv =  beta2(time,tstart1,tdur1,tdur2,tdur3)*Inv*Snv + beta2(time,tstart1,tdur1,tdur2,tdur3)*Iv*Snv 
    return(list(c(dSv, dSnv, dIv, dInv, dR, dCv, dCnv)))
  })
}

#### Optimising Phase 2 + 3 Duration ####

#Initial Conditions and Times
init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, R = 0, Cv= 0, Cnv = 0)
times <- seq(0,730,by = 1)

phase2dur <- seq(0, 300, by = 1)
phase3dur <- seq(0, 300, by = 1)
combdata <- expand.grid(phase2dur,phase3dur); colnames(combdata) <- c("phase2dur", "phase3dur")

start_time <- Sys.time()

phaseoptim <- data.frame(matrix(nrow = nrow(combdata), ncol = 6))

for(i in 1:nrow(combdata)) {
  temp <- numeric(6)
  init <- c(Sv = 0.20, Snv = 0.80-0.0001, Iv = 0, Inv = 0.0001, R = 0, Cv= 0, Cnv = 0)
  times <- seq(0,730,by = 1)
  parms = c(gamma = 1/(GenTime(4.6,2.4)), 
            tstart1 = 100, 
            tdur1 = 6*7,
            tdur2 = combdata$phase2dur[i],
            tdur3 = combdata$phase3dur[i],
            beta1_2 = 0)
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  out1$Iv <- out1$Iv/0.20
  out1$Inv <- out1$Inv/0.80
  out1$Cv <- out1$Cv/0.20
  out1$Cnv <- out1$Cnv/0.80
  temp[1] <- combdata$phase2dur[i]
  temp[2] <- combdata$phase3dur[i]
  temp[3] <- out1$time[which.max(out1$Iv)]
  temp[4] <- out1$Iv[which.max(out1$Iv)]
  temp[5] <- out1$Inv[which.max(out1$Iv)]
  temp[6] <- ifelse((out1$Iv[which.max(out1$Iv)] > out1$Iv[out1$time == 100]), "Yes", "No")
  phaseoptim[i,] <- temp
  print((i/nrow(combdata))*100)
}

colnames(phaseoptim) <- c("Phase2Duration", "Phase3Duration", "TimeofPeak", "InfVPeak", "InfnvPeak", "AboveLevel")

end_time <- Sys.time(); end_time - start_time

#HeatMap

phaseoptimtest <- phaseoptim
phaseoptimtest$Phase2Duration <- as.numeric(phaseoptimtest$Phase2Duration)
phaseoptimtest$Phase3Duration <- as.numeric(phaseoptimtest$Phase3Duration)
phaseoptimtest$InfVPeak <- as.numeric(phaseoptimtest$InfVPeak)

ggplot(phaseoptimtest, aes(x = Phase2Duration, y = Phase3Duration, z = InfVPeak)) +geom_contour_filled(breaks = c(0, 0.02, 0.04, 0.06,0.08,0.1,0.12,0.14,0.16,0.18)) + 
  scale_colour_distiller(palette = "Viridis", direction = 1)

#### Testing the Beta from the Optimiser ####

#Initial Conditions and Times
init <- c(Sv = 0.15, Snv = 0.85-0.0001, Iv = 0, Inv = 0.0001, R = 0, Cv= 0, Cnv = 0)
times <- seq(0,730,by = 1)

#Beta Functions are specified in the ODE
parms = c(gamma = 1/(GenTime(4.6,2.4)), 
          tstart1 = 100, 
          tdur1 = 6*7,
          tdur2 = 100,
          tdur3 = 100,
          beta1_2 = 0.064)

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))

#Extra Data Columns in Dataframe
out1$Iv <- out1$Iv/0.15
out1$Inv <- out1$Inv/0.85
out1$Cv <- out1$Cv/0.15
out1$Cnv <- out1$Cnv/0.85
out1$beta1 <- beta1(seq(0,730), 100, (6*7), 10, 10, 0.064)
out1$beta2 <- beta2(seq(0,730), 100, (6*7), 10, 10)

#Cumulative Infections are scaled relative to the relative V and nV population sizes - 0.15 and 0.85
#You should not compare these numbers to the recovered fraction - they will not add up

colnames(out1) <- c("Time", "Suscv", "Suscnv","Infected_Is", "Infected_Ip", "Recov", "CumV", "Cumnv", "Beta1", "Beta2")

statsinfec <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Is", "Infected_Ip"))
statsrecov <- melt(out1, id.vars = c("Time"), measure.vars = c("Recov"))
statsbeta <- melt(out1, id.vars = c("Time"), measure.vars = c("Beta1", "Beta2"))

#### Aggregated Plots

pinfv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Is",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of Vulnerable Infected") + scale_y_continuous(limits = c(0,0.1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14), axis.title.y=element_text(size=11),axis.title.x= element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

pinfnv <- ggplot(data = statsinfec[statsinfec$variable == "Infected_Ip",], aes(x = (Time), y = value)) + geom_line(size = 1.02, stat = "identity", col = "darkred") +
  labs(x ="Time (Days)", y = "Proportion of non-Vulnerable Infected") + scale_y_continuous(limits = c(0,0.1) ,  expand = c(0,0)) +
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
