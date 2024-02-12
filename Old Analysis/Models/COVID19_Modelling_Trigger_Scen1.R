rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot")

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### SCENARIO 1.1 - 9 Weeks #### 
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (2*(1/(GenTime(6,2))))*0.5,
         (2*(1/(GenTime(6,2)))))
}

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - mu*I
    dR = mu*I 
    dC = betastatdecrease(time,int_timestart,int_timeend)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,365,by = 1)

trigday <- seq(1,100, by = 1)
stats <- data.frame(matrix(nrow = length(trigday), ncol = 4)); colnames(stats) <- c("TrigDay", "TimePeak", "PeakInf", "CumInf")

for(i in 1:length(trigday)) {
  parms = c(mu = 1/(GenTime(6,2)), int_timestart = trigday[i], int_timeend = trigday[i]+(9*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  stats[i,1] <- parms[2]
  stats[i,2] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats[i,3] <- out1[,3][which(out1[,3] == max(out1[,3]))] #No Inf
  stats[i,4] <- max(out1[,5]) 
  print(i/length(trigday))
}

stats[which.min(stats$PeakInf),]

parms = c(mu = 1/(GenTime(6,2)),
          int_timestart = as.numeric(stats[which.min(stats$PeakInf),][1]), 
          int_timeend = as.numeric(stats[which.min(stats$PeakInf),][1])+(9*7))

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betastatdecrease(times, as.numeric(parms[2]), as.numeric(parms[3]))

#Plotting 

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
                  data.frame("Compartment" = "Susc", "Times" = out1[,1], "Prev" = out1[,2]),
                  data.frame("Compartment" = "Rec", "Times" = out1[,1], "Prev" = out1[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = out1[,1], "Prev" = out1[,5]),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,6]))

p11 <- ggplot(data = outdata1[outdata1$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(expand = c(0, 0)) + annotate("rect", xmin = as.numeric(parms[2]), xmax = as.numeric(parms[3]), ymin = 0, ymax = 0.15, fill = "darkred", alpha = .2)

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))


#### SCENARIO 1.2 - 12 Weeks #### 
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (2*(1/(GenTime(6,2))))*0.625,
         (2*(1/(GenTime(6,2)))))
}

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - mu*I
    dR = mu*I 
    dC = betastatdecrease(time,int_timestart,int_timeend)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,365,by = 1)

trigday <- seq(1,100, by = 1)
stats <- data.frame(matrix(nrow = length(trigday), ncol = 4)); colnames(stats) <- c("TrigDay", "TimePeak", "PeakInf", "CumInf")

for(i in 1:length(trigday)) {
  parms = c(mu = 1/(GenTime(6,2)), int_timestart = trigday[i], int_timeend = trigday[i]+(12*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  stats[i,1] <- parms[2]
  stats[i,2] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats[i,3] <- out1[,3][which(out1[,3] == max(out1[,3]))] #No Inf
  stats[i,4] <- max(out1[,5]) 
  print(i/length(trigday))
}

stats[which.min(stats$PeakInf),]

parms = c(mu = 1/(GenTime(6,2)),
          int_timestart = as.numeric(stats[which.min(stats$PeakInf),][1]), 
          int_timeend = as.numeric(stats[which.min(stats$PeakInf),][1])+(12*7))

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betastatdecrease(times, as.numeric(parms[2]), as.numeric(parms[3]))

#Plotting 

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
                  data.frame("Compartment" = "Susc", "Times" = out1[,1], "Prev" = out1[,2]),
                  data.frame("Compartment" = "Rec", "Times" = out1[,1], "Prev" = out1[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = out1[,1], "Prev" = out1[,5]),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,6]))

p11 <- ggplot(data = outdata1[outdata1$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(expand = c(0, 0)) + annotate("rect", xmin = as.numeric(parms[2]), xmax = as.numeric(parms[3]), ymin = 0, ymax = 0.15, fill = "darkred", alpha = .2)

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))

#### SCENARIO 1.3 - 18 Weeks #### 
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (2*(1/(GenTime(6,2))))*0.75,
         (2*(1/(GenTime(6,2)))))
}

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - mu*I
    dR = mu*I 
    dC = betastatdecrease(time,int_timestart,int_timeend)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,365,by = 1)

trigday <- seq(1,100, by = 1)
stats <- data.frame(matrix(nrow = length(trigday), ncol = 4)); colnames(stats) <- c("TrigDay", "TimePeak", "PeakInf", "CumInf")

for(i in 1:length(trigday)) {
  parms = c(mu = 1/(GenTime(6,2)), int_timestart = trigday[i], int_timeend = trigday[i]+(18*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  stats[i,1] <- parms[2]
  stats[i,2] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats[i,3] <- out1[,3][which(out1[,3] == max(out1[,3]))] #No Inf
  stats[i,4] <- max(out1[,5]) 
  print(i/length(trigday))
}

stats[which.min(stats$PeakInf),]

parms = c(mu = 1/(GenTime(6,2)),
          int_timestart = as.numeric(stats[which.min(stats$PeakInf),][1]), 
          int_timeend = as.numeric(stats[which.min(stats$PeakInf),][1])+(18*7))

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betastatdecrease(times, as.numeric(parms[2]), as.numeric(parms[3]))

#Plotting 

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
                  data.frame("Compartment" = "Susc", "Times" = out1[,1], "Prev" = out1[,2]),
                  data.frame("Compartment" = "Rec", "Times" = out1[,1], "Prev" = out1[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = out1[,1], "Prev" = out1[,5]),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,6]))

p11 <- ggplot(data = outdata1[outdata1$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(expand = c(0, 0)) + annotate("rect", xmin = as.numeric(parms[2]), xmax = as.numeric(parms[3]), ymin = 0, ymax = 0.15, fill = "darkred", alpha = .2)

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))
#### SCENARIO 1.4 - 36 Weeks #### 
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (2*(1/(GenTime(6,2))))*0.875,
         (2*(1/(GenTime(6,2)))))
}

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - mu*I
    dR = mu*I 
    dC = betastatdecrease(time,int_timestart,int_timeend)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,365,by = 1)

trigday <- seq(1,100, by = 1)
stats <- data.frame(matrix(nrow = length(trigday), ncol = 4)); colnames(stats) <- c("TrigDay", "TimePeak", "PeakInf", "CumInf")

for(i in 1:length(trigday)) {
  parms = c(mu = 1/(GenTime(6,2)), int_timestart = trigday[i], int_timeend = trigday[i]+(36*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  stats[i,1] <- parms[2]
  stats[i,2] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats[i,3] <- out1[,3][which(out1[,3] == max(out1[,3]))] #No Inf
  stats[i,4] <- max(out1[,5]) 
  print(i/length(trigday))
}

stats[which.min(stats$PeakInf),]

parms = c(mu = 1/(GenTime(6,2)),
          int_timestart = as.numeric(stats[which.min(stats$PeakInf),][1]), 
          int_timeend = as.numeric(stats[which.min(stats$PeakInf),][1])+(36*7))

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betastatdecrease(times, as.numeric(parms[2]), as.numeric(parms[3]))

#Plotting 

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
                  data.frame("Compartment" = "Susc", "Times" = out1[,1], "Prev" = out1[,2]),
                  data.frame("Compartment" = "Rec", "Times" = out1[,1], "Prev" = out1[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = out1[,1], "Prev" = out1[,5]),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,6]))

p11 <- ggplot(data = outdata1[outdata1$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(expand = c(0, 0)) + annotate("rect", xmin = as.numeric(parms[2]), xmax = as.numeric(parms[3]), ymin = 0, ymax = 0.15, fill = "darkred", alpha = .2)

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))