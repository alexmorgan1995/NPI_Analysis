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

plot(seq(0,365), betapulse(seq(0,365), 41, 41+(12*7), 12))

SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betapulse(time, int_timestart, int_timeend, timegap)*S*I + mu*I 
    dI = betapulse(time, int_timestart, int_timeend, timegap)*S*I - mu*I
    dC = betapulse(time, int_timestart, int_timeend, timegap)*S*I
    return(list(c(dS,dI,dC)))
  })
}

SIR2 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betapulse(time, int_timestart, int_timeend, timegap)*S*I
    dI = betapulse(time, int_timestart, int_timeend, timegap)*S*I - mu*I
    dR = mu*I 
    dC = betapulse(time, int_timestart, int_timeend, timegap)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, C = 0)
init1 <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)

times <- seq(0,730,by = 1)

parms = c(mu = 1/(GenTime(4.6,2.4)),
          int_timestart = 100, 
          int_timeend = 100+(12*7),
          timegap = 100)

out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out1$beta <- betapulse(times, as.numeric(parms[2]), as.numeric(parms[3]), as.numeric(parms[4]))

out2 <- data.frame(ode(y = init1, func = SIR2, times = times, parms = parms))
out2$beta <- betapulse(times, as.numeric(parms[2]), as.numeric(parms[3]), as.numeric(parms[4]))

#Plotting 

outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = out1[,1], "Prev" = out1[,3], "State" = "SIS"),
                  data.frame("Compartment" = "Beta", "Times" = out1[,1], "Prev" = out1[,5], "State" = "SIS"),
                  data.frame("Compartment" = "Infec", "Times" = out2[,1], "Prev" = out2[,3], "State" = "SIR"),
                  data.frame("Compartment" = "Beta", "Times" = out2[,1], "Prev" = out2[,6], "State" = "SIR"))

p11 <- ggplot(data = outdata1[outdata1$Compartment == "Infec",], aes(x = Times, y = Prev, col = State)) + geom_line(size = 1.05) +
  labs(y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.5) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(expand = c(0, 0)) 

p21 <- ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

plot_grid(p11, NULL, p21, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))