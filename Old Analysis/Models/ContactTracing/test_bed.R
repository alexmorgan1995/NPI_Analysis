setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/ContactTracing") 
rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("Cairo")

#### Model Functions - Gen Time/ Betas/ODEs ####
#Function to calculate the Generation time (1/gamma) parameter

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

beta <- function(time, tstart1) {
  gamma <- 1/(GenTime(3.3,2.8))
  ifelse((time >= tstart1 & time <= Inf), 
         1.2*gamma, #Phase 2
         1.7*gamma #Phase 1
  )
} 

plot(beta(seq(0,730), 71))


lin_cont <- function(x,m,comp) {
    y = m*(x-comp) + 0
  return(y)
}


inc <- data.frame(0,0) # initiate at 0
cont_trac <- data.frame(0,0) # initiate at 0

#inc <- data.frame("time" = seq(0, 15, by = 1), "inc" = seq(0, 1500, by =  100)/(5.5*10^6))
#inc$time[which(inc$inc >= 1000/(5.5*10^6))[1]]
#lin_cont(seq(0,20),  50, inc$time[which(inc$inc >= 1000/(5.5*10^6))[1]])

SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    ifelse(is.na(inc[,1][which(inc[,2] >= 1000/(5.5*10^6))[1]]) == TRUE,
           timecomp <<- Inf,
           timecomp <<- inc[,1][which(inc[,2] >= 1000/(5.5*10^6))[1]]
           )

    if(time >= timecomp) {
      
      beta <- beta(time,tstart1)
      inctime <- c(time, beta*I*S)
      inc <<- rbind(inc, inctime)
      print(time)
      cont_trac <<- rbind(cont_trac, lin_cont(time, 50, timecomp))
      
      dS = - beta*I*S + zeta*R
      dI = beta*I*S - gamma*I
      dR = gamma*I - zeta*R 
      
    } else {
      
      timecomp <<- inc[,2][which(inc[,1] >= 1000/(5.5*10^6))[1]]
      
      beta <- beta(time,tstart1)
      print(c(time, "PE"))
      cont_trac <<- rbind(cont_trac, 0)
      inctime <- c(time, beta*I*S)
      inc <<- rbind(inc, inctime)
      
      dS = - beta*I*S + zeta*R
      dI = beta*I*S - gamma*I
      dR = gamma*I - zeta*R 
      
    }
    return(list(c(dS, dI, dR)))
  })
}


N <- 5.5*10^6
init <- c(S = (N-1)/N,
          I = 1/N,   
          R= 0)

times <- seq(0, 478, by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71)

out1 <- data.frame(euler(y = init, func = SIRS, times = times, parms = parms))
