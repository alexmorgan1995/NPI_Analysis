setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/ContactTracing") 
rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("Cairo")

#### Model Functions - Gen Time/ Betas/ODEs ####
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

beta <- function(time, tstart1) {
  gamma <- 1/(GenTime(3.3,2.8))
  ifelse((time >= tstart1 & time <= Inf), 
         0.85*gamma, #Phase 2 - Before this was R0 = 0.9
         1.7*gamma #Phase 1
  )
} 

plot(beta(seq(0,730), 71))

#### Baseline w/ Lockdown #####

inc <- data.frame(0, 0) 

SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    beta <- beta(time,tstart1)
    
    dS = - beta*I*S + zeta*R
    dI = beta*I*S - gamma*I
    dR = gamma*I - zeta*R 
      
    inc_tim_inc <- c(time, beta*I*S)
    inc <<- rbind(inc, inc_tim_inc) # Tracking the daily incidence
     
    return(list(c(dS, dI, dR)))
  })
}

#Parameters
N <- 5.5*10^6
init <- c(S = (N-(N*0.0001))/N,
          I = (N*0.0001)/N,   
          R= 0)

times <- seq(0, 478, by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71,
          N = 5.5*10^6)

#Model Run
out1 <- data.frame(euler(y = init, func = SIRS, times = times, parms = parms))

plot(out1$time, out1$I*N)

plot(inc[,1], inc[,2]*N)

#Find the Downwards Point at which inc = 1000
test <- inc[,1][which(inc[,2]*N > 1000)]
test[length(test)] # This is the trigger day for contact tracing capacity increase


#### Model to Find Point of Intersect after day 176 ####

lin_cont <- function(x,m,comp) {
  y = m*(x-comp) + 0
  return(y)
}

timecomp <- test[length(test)]

inc_cont <- data.frame(0, 0, 0, 0, 0) # initiate at 0

SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    timecomp = timecomp
    timeequal = 165
    
    if(time >= timecomp) {
      
      if(time >= timeequal) { # if time is equal
        
        beta <- beta(time,tstart1)*(1-(1*efficacy))
        
        
        print(c(time, "Equal"))
        
        dS = - beta*I*S + zeta*R
        dI = beta*I*S - gamma*I
        dR = gamma*I - zeta*R 
        
        
        inctime <- c(time+1, beta*I*S*N, beta*I*S*N ,beta,  (beta*I*S*N)/ (beta*I*S*N))
        inc_cont <<- rbind(inc_cont, inctime)
        
      } else{ #if inc is over 1000
        
        beta <- beta(time,tstart1)*(1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*efficacy)
       
        print(c(time, "During"))
        
        dS = - beta*I*S + zeta*R
        dI = beta*I*S - gamma*I
        dR = gamma*I - zeta*R 
        
        
        inctime <- c(time+1, beta*I*S*N, lin_cont(time, 50, timecomp), beta, lin_cont(time, 50, timecomp)/ (beta*I*S*N))
        inc_cont <<- rbind(inc_cont, inctime)
      }

      } else { #if inc < 1000
      
      beta <- beta(time,tstart1)
      
      print(c(time, "Before"))
      
      dS = - beta*I*S + zeta*R
      dI = beta*I*S - gamma*I
      dR = gamma*I - zeta*R 
      
      inctime <- c(time+1, beta*I*S*N , 0, beta, 0/ (beta*I*S*N))
      print(inctime)
      inc_cont <<- rbind(inc_cont, inctime)
      
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
          N = 5.5*10^6,
          efficacy = 0.8)

out1 <- data.frame(euler(y = init, func = SIRS, times = times, parms = parms))

inc_cont <- inc_cont[-1,]
inc_cont[,1][which(inc_cont[,2] < 10)][2]
inc_cont[,1][which(inc_cont[,2] < 10)][2] - timecomp
inc_cont[1,1] <- 0; inc_cont[1,2] <- 0

plot(out1$time, out1$I, lty = 2)

plot(inc_cont$time,inc_cont$incidence, lty = 2)


statsinc <- melt(inc_cont, id.vars =  c("time"), measure.vars = c("incidence", "contacttraccap"))

ggplot(data = statsinc, aes(x = (time), y = value, col = variable))+ theme_bw() +
  labs(x ="Time (Days)", y = "Daily Incidence", color = "Population") + scale_y_continuous(limits = c(0,1000), expand = c(0,0)) + 
  scale_x_continuous(limits= c(177, 210), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_line(size = 1.02, stat = "identity")

ggplot(data = out1, aes(x = (time), y = I*N))+ theme_bw() +
  labs(x ="Time (Days)", y = "Total Infected", color = "Population") + scale_y_continuous(limits = c(0, 20000),expand = c(0,0)) + 
  scale_x_continuous(limits = c(176, 210), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_line(size = 1.02, stat = "identity")



ggplot(data = inc_cont, aes(x = (time), y = incidence))+ theme_bw() +
  labs(x ="Time (Days)", y = "Total Infected", color = "Population") + scale_y_continuous(limits = c(0, 25000),expand = c(0,0)) + 
  scale_x_continuous( expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_line(size = 1.02, stat = "identity")
