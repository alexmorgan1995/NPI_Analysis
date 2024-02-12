rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot"); library("shiny")

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### SCENARIO 1 - 12 Weeks #### 

betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*0.625,
         (1.5*(1/(GenTime(4.6,2.4)))))
}

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - gamma*I - delta*I
    dC = delta*I - mu*C
    dR = gamma*I + mu*C
    return(list(c(dS,dI,dC,dR)))
  })
}

init <- c(S = 0.9999, I = 0.0001,  C = 0, R = 0)
times <- seq(0,365,by = 1)
trigday <- seq(1,100, by = 1)
stats1 <- data.frame(matrix(nrow = length(trigday), ncol = 4))
colnames(stats1) <- c("TrigDay", "TimePeak", "PeakInf", "CumInf")

for(i in 1:length(trigday)) {
  parms = c(gamma = 1/(GenTime(4.6,2.4)),
            delta = (1/(GenTime(4.6,2.4)))/200,
            mu = 0.1,
            int_timestart = trigday[i], 
            int_timeend = trigday[i]+(12*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  stats1[i,1] <- parms[2]
  stats1[i,2] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats1[i,3] <- out1[,3][which(out1[,3] == max(out1[,3]))] #No Inf
  stats1[i,4] <- out1[,1][which(out1[,3] == max(out1[,3]))] #Time
  stats1[i,5] <- out1[,3][which(out1[,3] == max(out1[,3]))] 
  stats1[i,6] <- max(out1[,5]) 
  print(i/length(trigday))
}

stats1[which.min(stats1$PeakInf),]

parms = c(gamma = 1/(GenTime(4.6,2.4)),
          int_timestart = as.numeric(stats1[which.min(stats1$PeakInf),][1]), 
          int_timeend = as.numeric(stats1[which.min(stats1$PeakInf),][1])+(12*7))

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betastatdecrease(times, as.numeric(parms[2]), as.numeric(parms[3]))

#Plotting 

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3]),
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
