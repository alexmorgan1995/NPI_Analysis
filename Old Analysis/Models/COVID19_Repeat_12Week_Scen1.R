rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot")

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### SCENARIO 1 - 12 Weeks #### 


betapulse <- function(time, int_timestart, int_timeend, timegap) {
  ifelse((time >= int_timestart & time <= int_timestart+(12*7)) | 
           (time >= (int_timestart + (12*7) + timegap) & time <= (int_timestart + (12*7) + timegap)+(12*7)), 
         1.5*(1/(GenTime(4.6,2.4)))*0.4, 
         1.5*(1/(GenTime(4.6,2.4))))
}

plot(times, betapulse(times, 41, 41+(12*7), 12))

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betapulse(time, int_timestart, int_timeend, timegap)*S*I
    dI = betapulse(time, int_timestart, int_timeend, timegap)*S*I - mu*I
    dR = mu*I 
    dC = betapulse(time, int_timestart, int_timeend, timegap)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,730,by = 1)
trigday <- seq(1,200, by = 1)
timegap <- seq(1,120, by = 1)
combparm <- expand.grid(trigday,timegap)
colnames(combparm) <- c("trigday", "timegap")

stats1 <- data.frame(matrix(nrow = length(combparm), ncol = 5)); colnames(stats1) <- c("TrigDay", "TimeGap", "TimePeak", "PeakInf", "CumInf")

for(i in 1:nrow(combparm)) {
  parms = c(mu = 1/(GenTime(4.6,2.4)), int_timestart = combparm[i,1], int_timeend = combparm[i,1]+(12*7), timegap = combparm[i,2])
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  stats1[i,1] <- combparm[i,1]
  stats1[i,2] <- combparm[i,2]
  stats1[i,3] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats1[i,4] <- out1[,3][which(out1[,3] == max(out1[,3]))] #No Inf
  stats1[i,5] <- max(out1[,5]) 
  print(i/nrow(combparm))
}

stats1[which.min(stats1$PeakInf),]

parms = c(mu = 1/(GenTime(4.6,2.4)),
          int_timestart = as.numeric(stats1[which.min(stats1$PeakInf),][1]), 
          int_timeend = as.numeric(stats1[which.min(stats1$PeakInf),][1])+(12*7),
          timegap = as.numeric(stats1[which.min(stats1$PeakInf),][2]))


parms = c(mu = 1/(GenTime(4.6,2.4)),
          int_timestart = 100, 
          int_timeend = 100+(12*7),
          timegap = 104)

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betapulse(times, as.numeric(parms[2]), as.numeric(parms[3]), as.numeric(parms[4]))

#Plotting 

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,6]))

p11 <- ggplot(data = outdata1[outdata1$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(expand = c(0, 0)) 

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))