rm(list=ls())

library("deSolve"); library("ggplot2"); library("cowplot")

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### SCENARIO 1 - 12 Weeks ####
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*0.4,
         (1.5*(1/(GenTime(4.6,2.4)))))
}

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - gamma*I - delta*I
    dC = delta*I - mu*C
    dR = gamma*I + mu*C
    dCumInf = betastatdecrease(time,int_timestart,int_timeend)*S*I
    dCumClin = delta*I
    return(list(c(dS,dI,dC,dR,dCumInf,dCumClin)))
  })
}

init <- c(S = 0.9999, I = 0.0001,  C = 0, R = 0, CumInf = 0, CumClin = 0)
times <- seq(0,730,by = 1)
trigday <- seq(1,200, by = 1)
stats1 <- data.frame(matrix(nrow = length(trigday), ncol = 7))
colnames(stats1) <- c("TrigDay", "TimePeakInf", "PeakInf_Inf", "TimePeakClin", "PeakInf_Clin", "CumInf","CumClin")

for(i in 1:length(trigday)) {
  parms = c(gamma = 1/(GenTime(4.6,2.4)),
            delta = (1/(GenTime(4.6,2.4)))/200,
            mu = 0.1,
            int_timestart = trigday[i],
            int_timeend = trigday[i] + (12*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  stats1[i,1] <- parms[4]
  stats1[i,2] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats1[i,3] <- out1[,3][which(out1[,3] == max(out1[,3]))] #No Inf
  stats1[i,4] <- out1[,1][which(out1[,4] == max(out1[,4]))]
  stats1[i,5] <- out1[,4][which(out1[,4] == max(out1[,4]))]
  stats1[i,6] <- max(out1[,6])
  stats1[i,7] <- max(out1[,7])
  print(i/length(trigday))
}

#### Peak Inf ####

stats1[which.min(stats1$PeakInf_Inf),]


parms = c(gamma = 1/(GenTime(4.6,2.4)),
          delta = (1/(GenTime(4.6,2.4)))/200,
          mu = 0.1,
          int_timestart = as.numeric(stats1[which.min(stats1$PeakInf_Inf),][1]),
          int_timeend = as.numeric(stats1[which.min(stats1$PeakInf_Inf),][1]) + (12*7))

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betastatdecrease(times, as.numeric(parms[4]), as.numeric(parms[5]))

#Plotting

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
                  data.frame("Compartment" = "ClinInf", "Times" = out1[,1], "Prev" = out1[,4]),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,8]))

p11 <- ggplot(out1, aes(x = time)) + geom_line(aes(y=I),size = 1.05, col = "darkblue") + geom_line(aes(y=(C*5000000)/6000), size=1.05, col="orange") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0, 0.15), expand = c(0,0), sec.axis = sec_axis(~.*6000, name = "Severe Cases")) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))

#### Peak Clin ####

stats1[which.min(stats1$PeakInf_Clin),]
parms = c(gamma = 1/(GenTime(4.6,2.4)),
          delta = (1/(GenTime(4.6,2.4)))/200,
          mu = 0.1,
          int_timestart = as.numeric(stats1[which.min(stats1$PeakInf_Clin),][1]),
          int_timeend = as.numeric(stats1[which.min(stats1$PeakInf_Clin),][1]) + (12*7))

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betastatdecrease(times, as.numeric(parms[4]), as.numeric(parms[5]))

#Plotting

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
                  data.frame("Compartment" = "ClinInf", "Times" = out1[,1], "Prev" = (out1[,4])*5000000),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,8]))

p11 <- ggplot(out1, aes(x = time)) + geom_line(aes(y=I),size = 1.05, col = "darkblue") + geom_line(aes(y=(C*5000000)/6000), size=1.05, col="orange") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.15), expand = c(0,0), sec.axis = sec_axis(~.*6000, name = "Severe Cases")) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11)) +
  scale_x_continuous(expand = c(0, 0))

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0))

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))