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

#ODEs for the Enhanced Shielding Model - Remainder is split into 3 sub-groups

#First if function if incidence is over 1000 initiate the if statement - initiate contact tracing function here at 0 but do not change anything
#Nested in this if function - have the contact tracing increase and the effect on beta with the efficacy modifier on the beta
#nested in this if function - if beta*I*S*conacttacing*fficacy > contact tracing -then i need to specify contact tracing to be = beta*I*S*contact tracing*effiacy
#have a function that increases bh 50 every day and affects the beta
#then also have an if function that states that if beta*I*S reaches 
#tokepp track of contact tracing and inidnence just have a <<-rbind at each timestep.

inc <- data.frame(0) # initiate at 0
cont_trac <- data.frame(0) # initiate at 0

lin_cont <- function(x,m) {
  y = m*x + 0
  return(y)
}


SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    if(beta*I*S > 1000/N) { # incidence over 1000 # find time at which cont = inc then from this point to infinity run this model 
      
      if(beta*I*S < cont_trac[length(cont_trac)]) { #If Incidence = Contacts
        
        beta <- beta(time,tstart1) # find time at which cont = inc then from this point to infinity run this model 
        
        cont_trac <<- rbind(cont_trac, beta*I*S)
        inc <<- rbind(inc, beta*I*S)
        
        dS = - beta*I*S + zeta*R
        dI = beta*I*S - gamma*I
        dR = gamma*I - zeta*R
        
      } else{ #If not at contacts = 
        
        cont_trac <<- rbind(cont_trac, lin_cont(time,m))# do soemthing like time - time at which incidence = 1000 - and take the 1st vector 
        
        beta <- beta(time,tstart1)*cont_trac[length(cont_trac)]*efficacy
        inc <<- rbind(inc, beta*I*S)
        
        dS = - beta*I*S + zeta*R
        dI = beta*I*S - gamma*I
        dR = gamma*I - zeta*R
        
      }
      
    } else{
      
      cont_trac <<- rbind(cont_trac, 0)
      
      beta <- beta(time,tstart1)
      inc <<- rbind(inc, beta*I*S)
      
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
          tstart1 = 71,
          m = 50,
          efficacy = 0.6)

out1 <- data.frame(euler(y = init, func = SIRS, times = times, parms = parms))



out1$inc <- c(0, diff(out1$S-out1$R*1/365)*-1)


####

linear <- function(x,m) {
  y = m*x + 0
  return(y)
}

linear(20, 50)
