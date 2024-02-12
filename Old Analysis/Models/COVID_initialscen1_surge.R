rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot")

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### Trigger Optim ####

betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*0.4,
         (1.5*(1/(GenTime(4.6,2.4)))))
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
times <- seq(0,730,by = 1)

trigday <- seq(1,300, by = 1)
stats <- data.frame(matrix(nrow = length(trigday), ncol = 4)); colnames(stats) <- c("TrigDay", "TimePeak", "PeakInf", "CumInf")

for(i in 1:length(trigday)) {
  parms = c(mu = 1/(GenTime(4.6,2.4)), 
            int_timestart = trigday[i], 
            int_timeend = trigday[i]+(12*7))
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



#### SCENARIO 1 - Beta tatic Decrease #### 

betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*0.4,
         (1.5*(1/(GenTime(4.6,2.4)))))
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
times <- seq(0,730,by = 1)
duration <- seq(1,12)

stats1 <- data.frame(matrix(nrow = 0, ncol = 4))


for(i in 1:length(duration)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 4))
  parms = c(mu = 1/(GenTime(4.6,2.4)), 
            int_timestart = 100, 
            int_timeend = 100 + (duration[i]*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- paste0(as.character(duration[i]), " Weeks")
  temp[,2] <- out1$time
  temp[,3] <- out1$I #Time
  temp[,4] <- betastatdecrease(times, 
                               as.numeric(parms[2]), 
                               as.numeric(parms[3]))
  stats1 <- rbind(stats1, temp)
}

colnames(stats1) <- c("WeekDuration", "Time", "Infec", "Beta")

stats1$WeekDuration <- factor(stats1$WeekDuration, levels = unique(stats1$WeekDuration))

#Plotting 

p11 <- ggplot(data = stats1, aes(x = (Time), y = Infec, col = WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

p21 <- ggplot(data = stats1, aes(x = (Time), y = Beta, col = WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "??") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + scale_color_brewer(palette="Paired")

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(1.5, -0.05, 1))

#Data Table 

surge <- data.frame(matrix(nrow = 0, ncol = 5))

for(i in 1:length(duration)) { 
  test <- stats1[stats1$WeekDuration == as.character(duration[i]*7),][which.max(stats1$Infec[stats1$WeekDuration == as.character(duration[i]*7)]),]
  surge <- rbind(surge, test)
  print(test)
}

surge$TimeBetPeaks <- surge$Time - 100
surge$RelInc <- surge$Infec / 0.0181504276
surge$dbltime <- c(35,53,71,89,107,125,143,160,177,193,209,224)
surge$WeekDuration <- c(1,2,3,4,5,6,7,8,9,10,11,12)

ggplot(surge, aes(x = as.numeric(WeekDuration), y=RelInc)) + geom_line(col = "darkblue", size = 1.05) +
  labs(y = "Proportional Increase in I(t) between 1st and 2nd Peak") + scale_y_continuous(limits = c(0,3.5) ,  expand = c(0,0)) +
  scale_x_continuous(name = "Duration of Intervention (Weeks)", breaks = seq(1,12, by = 1), 
                     limits = c(0,12))

ggplot(surge, aes(x = as.numeric(WeekDuration), y=dbltime)) + geom_line(col = "darkblue", size = 1.05) +
  labs(x ="Duration of Intervention (Weeks)", y = "Doubling Time (Days)") + scale_y_continuous(limits = c(0, 250), expand = c(0,0))  +
  scale_x_continuous(name = "Duration of Intervention (Weeks)", breaks = seq(1,12, by = 1), 
                     limits = c(0,12))


#### SIS ####

surgesis <- data.frame("Time" = seq(1,12),
                       "Doubling_Time" = c(27,40,52,65,78,91,103,116,129,141,154,167))
  
  
ggplot(surgesis, aes(x = as.numeric(Time), y=Doubling_Time)) + geom_line(col = "darkblue", size = 1.05) +
  labs(x ="Duration of Intervention (Weeks)", y = "Doubling Time (Days)") + scale_y_continuous(limits = c(0, 250), expand = c(0,0))  +
  scale_x_continuous(name = "Duration of Intervention (Weeks)", breaks = seq(1,12, by = 1), 
                     limits = c(0,12))
