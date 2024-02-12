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

SIRS1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    beta <- beta(time,tstart1, R0)
    
    dS = - beta*I*S + zeta*R
    dI = beta*I*S - gamma*I
    dR = gamma*I - zeta*R 
    
    inc_tim_inc <- c(time, beta*I*S)
    inc_cont <<- rbind(inc_cont, inc_tim_inc) # Tracking the daily incidence
    
    return(list(c(dS, dI, dR)))
  })
}



#Parameters
N <- 5.5*10^6
init <- c(S = (N-(N*0.0001))/N, I = (N*0.0001)/N, R= 0)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          N = 5.5*10^6,
          R0 = 0.85)


#Baseline 
times <- seq(0, 478, by = 1)
inc_cont <- data.frame(0, 0)
out1 <- data.frame(rk(y = init, func = SIRS1, times = times, parms = parms, method = "rk4"))
inc_cont <- inc_cont[!duplicated(trunc(inc_cont$X0)), ]
