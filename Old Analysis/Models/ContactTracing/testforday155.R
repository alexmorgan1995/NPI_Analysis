setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/ContactTracing") 
rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("Cairo")

#### Model Functions - Gen Time/ Betas/ODEs ####

#Generation Time

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Model Betas

beta <- function(time, tstart1, R0) {
  gamma <- 1/(GenTime(3.3,2.8))
  ifelse((time >= tstart1 & time <= Inf), 
         R0*gamma, #Phase 2 - Before this was R0 = 0.9
         1.7*gamma #Phase 1
  )
} 

#Baseline SIRS To identify time at inc = 1000

SIRS1 <- function(time, init, parms) {
  
  beta <- beta(time, parms["tstart1"], parms["R0"])
  
  dS = - beta*init["I"]*init["S"]
  dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
  dR = parms["gamma"]*init["I"]
  
  inc_tim_inc <- c(time, beta*init["I"]*init["S"]*parms["N"])
  inc_cont <<- rbind(inc_cont, inc_tim_inc) # Tracking the daily incidence
  inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
  
  return(list(c(dS, dI, dR)))
}


N <- 5.5*10^6

init <- c(S = (N-(N*0.00001))/N, I = (N*0.00001)/N, R= 0)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 64, 
          N = 5.5*10^6, 
          R0 = 0.85)

#Baseline 
times <- seq(0, 200, by = 0.5)
inc_cont <<- data.frame(0, 0)
out1 <- data.frame(rk(y = init, func = SIRS1, times = times, parms = parms, method = "rk4"))

ggplot(data = inc_cont, aes(x = X0, y = X0.1))+ theme_bw() +
  labs(x ="Time (Days)", y = "Daily Incidence", color = "Population") + scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_line(size = 1.02, stat = "identity")

timecomp <<- inc_cont[,1][which( inc_cont[,1] > parms["tstart1"] & abs((inc_cont[,2] - 1000) - 0) == min(abs((inc_cont[,2] - 1000) - 0)))]
t <- data.frame(abs((inc_cont[,2] - 1000) - 0) == min(abs((inc_cont[,2] - 1000) - 0)))
