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
         0.9*gamma, #Phase 2
         1.7*gamma #Phase 1
  )
} 

plot(beta(seq(0,730), 71))


lin_cont <- function(x,m,comp) {
    y = m*(x-comp) + 0
  return(y)
}


inc_cont <- data.frame(0,0, 0) # initiate at 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
#inc <- data.frame("time" = seq(0, 15, by = 1), "inc" = seq(0, 1500, by =  100)/(5.5*10^6))
#inc$time[which(inc$inc >= 1000/(5.5*10^6))[1]]
#lin_cont(seq(0,20),  50, inc$time[which(inc$inc >= 1000/(5.5*10^6))[1]])

SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    timecomp = 176
    
    

    ifelse(is.na(inc_cont[,1][which(inc_cont[,2] < (inc_cont[nrow(inc_cont),2]/N))[1]]) == TRUE,
           equals1 <<- Inf,
           equals1 <<- inc_cont[,1][which(inc_cont[,2] < (inc_cont[nrow(inc_cont),2]/N))[1]]
    )
    
    if(time >= timecomp) {
      
      beta <- beta(time,tstart1)*(1-(inc_cont[nrow(inc_cont),3]/N)*0.6)

      inctime <- c(time, beta*I*S, lin_cont(time, 50, timecomp))
      inc_cont <<- rbind(inc_cont, inctime)
      print(time)
    
      dS = - beta*I*S + zeta*R
      dI = beta*I*S - gamma*I
      dR = gamma*I - zeta*R 
      
    } else {
      
      beta <- beta(time,tstart1)

      print(c(time, "PE"))
      
      inctime <- c(time, beta*I*S , 0)
      inc_cont <<- rbind(inc_cont, inctime)

      dS = - beta*I*S + zeta*R
      dI = beta*I*S - gamma*I
      dR = gamma*I - zeta*R 
      
    }
    return(list(c(dS, dI, dR)))
  })
}



N <- 5.5*10^6
init <- c(S = (N-(N*0.0001))/N,
          I = (N*0.0001)/N,   
          R= 0)

times <- seq(0, 478, by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          effic = 0.6,
          N = 5.5*10^6)

out1 <- data.frame(euler(y = init, func = SIRS, times = times, parms = parms))

plot(out1$time, out1$I)

inc_cont$X0.1 <- inc_cont$X0.1*N

plot(inc_cont$X0, inc_cont$X0.1)



#### Equal Model ####


lin_cont <- function(x,m,comp) {
  y = m*(x-comp) + 0
  return(y)
}


inc_cont <- data.frame(0,0, 0, 0,0) # initiate at 0

#inc <- data.frame("time" = seq(0, 15, by = 1), "inc" = seq(0, 1500, by =  100)/(5.5*10^6))
#inc$time[which(inc$inc >= 1000/(5.5*10^6))[1]]
#lin_cont(seq(0,20),  50, inc$time[which(inc$inc >= 1000/(5.5*10^6))[1]])

SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    timecomp = 176
    
    timeequal = 187
    
    
    if(time >= timecomp) {
      
      if(time >= timeequal) { # if time is equal
        
        beta <- beta(time,tstart1)*(1-(1*efficacy))
        
        inctime <- c(time, beta*I*S*N, beta*I*S*N , beta*I*S*N/ beta*I*S*N)
        inc_cont <<- rbind(inc_cont, inctime)
        
        print(c(time, "Equal"))
        
        dS = - beta*I*S + zeta*R
        dI = beta*I*S - gamma*I
        dR = gamma*I - zeta*R 
        
      } else{ #if inc is over 1000
        
        beta <- beta(time,tstart1)*(1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]*N))*efficacy)
        
        inctime <- c(time, beta*I*S*N, lin_cont(time, 50, timecomp), beta, lin_cont(time, 50, timecomp)/ beta*I*S*N)
        inc_cont <<- rbind(inc_cont, inctime)
        print(time)
        
        dS = - beta*I*S + zeta*R
        dI = beta*I*S - gamma*I
        dR = gamma*I - zeta*R 
      }

      } else { #if inc < 1000
      
      beta <- beta(time,tstart1)
      
      print(c(time, "PE"))
      
      inctime <- c(time, beta*I*S*N , 0, beta, 0/ beta*I*S*N)
      inc_cont <<- rbind(inc_cont, inctime)
      
      dS = - beta*I*S + zeta*R
      dI = beta*I*S - gamma*I
      dR = gamma*I - zeta*R 
      
    }
    
    return(list(c(dS, dI, dR)))
  })
}

colnames(inc_cont) <- c("time", "incidence", "contacttraccap", "beta", "con/inc rat")

N <- 5.5*10^6
init <- c(S = (N-(N*0.0001))/N,
          I = (N*0.0001)/N,   
          R= 0)

times <- seq(0, 478, by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          effic = 0.5,
          N = 5.5*10^6,
          efficacy = .8)

out1 <- data.frame(euler(y = init, func = SIRS, times = times, parms = parms))

plot(out1$time, out1$I)

inc_cont$X0.1 <- inc_cont$X0.1*N

plot(inc_cont$X0, inc_cont$X0.1)

