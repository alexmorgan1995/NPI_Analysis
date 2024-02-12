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

beta <- function(R0) {
  gamma <- 1/(GenTime(3.3,2.8))
  return(R0*gamma)
} 

#Linear Increase in contact tracing

lin_cont <- function(x,m,comp) {
  y = m*(x) + 0
  return(y)
}

#SIRS Model To Identify where contact tracing = incidence

SIRS2 <- function(time, init, parms) {
  
  beta <- beta(parms["R0"])
  linear <<- lin_cont(time, parms["ramp_rate"])
  
  incidence <- beta*init["I"]*init["S"]*parms["N"]
  inctime <- c(time, incidence, linear, beta, linear/ (incidence))
  inc_cont <<- rbind(inc_cont, inctime)
  inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
  
  beta <- beta(parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
  
  dS = - beta*init["I"]*init["S"] 
  dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
  dR = parms["gamma"]*init["I"]
  
  return(list(c(dS, dI, dR)))
}


#SIRS Model used for Final Model Output

SIRS3 <- function(time, init, parms) {
  
  timeequal = timeequal
  
  if(time >= timeequal) { # if time is equal
    
    beta <- beta(parms["R0"])*(1-(1*parms["efficacy"]))
    incidence <- beta*init["I"]*init["S"]*parms["N"]
    
    inctime <- c(time, incidence, incidence, beta, 1)
    inc_cont <<- rbind(inc_cont, inctime)
    inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
    
    dS = - beta*init["I"]*init["S"]
    dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
    dR = parms["gamma"]*init["I"] 
    
  } else{ #if inc is over 1000
    
    beta <- beta(parms["R0"])
    linear <<- lin_cont(time, parms["ramp_rate"])
    
    incidence <- beta*init["I"]*init["S"]*parms["N"]
    inctime <- c(time, incidence, linear, beta, 
                 linear/ (incidence))
    inc_cont <<- rbind(inc_cont, inctime)
    inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
    
    beta <- beta(parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
    
    dS = - beta*init["I"]*init["S"] 
    dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
    dR = parms["gamma"]*init["I"]
    
  }
  
  return(list(c(dS, dI, dR)))
}

#### Big Function ####

N <- 5.5*10^6
init <- c(S = (N-(N*0.0001))/N, I = (N*0.0001)/N, R= 0)

contact_trac_run <- function(efficacy, ramp_rate, R0, target, iinput, rinput) {
  
  #Parameters
  N <- 5.5*10^6
  init <- c(S = (N-(iinput + rinput))/N, I = iinput/N, R = rinput/N)
  
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            N = 5.5*10^6,
            efficacy = efficacy,
            ramp_rate = ramp_rate, 
            R0 = R0,
            target = target)
  
  #Intersect
  inc_cont <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"] , 0, 0, 0)
  times <- seq(0, 200, by = 1)
  out1 <- data.frame(rk(y = init, func = SIRS2, times = times, parms = parms, method = "rk4"))
  
  t <- inc_cont[,1][which(inc_cont[,5] < 1 & inc_cont[,5] > 0.5)]
  t1 <- split(t, cumsum(c(1, diff(t) != 0.5)))[[1]]
  timeequal <<- t1[length(t1)]
  
  
  #Actually Run the Model 
  inc_cont <<- data.frame(0, beta(parms["R0"])*init["I"]*init["S"]*parms["N"] , 0, 0, 0)
  times <- seq(0, 200, by = 1)
  out1 <<- data.frame(rk(y = init, func = SIRS3, times = times, parms = parms, method = "rk4"))
  
  #inc_cont[1,1] <<- 0; inc_cont[1,2] <<- 0
  colnames(inc_cont) <<- c("time", "incidence", "contacttraccap", "beta", "con/inc rat")
  
  timesince <- inc_cont[,1][which(inc_cont[,2] < target)][1]
  return(c("Time1000" = timeequal, "RelTimetoTarget" = timesince))
}

#Specify the efficacy (as a fraction), ramp_rate, lockdown R0 and target incidence
#inc_cont will always output during the run - view for a more detailed overview of the model run.

test <- contact_trac_run(0.8, 50, 0.9, 10, 10880, 459384.8)

statsinc <- melt(inc_cont, id.vars =  c("time"), measure.vars = c("incidence", "contacttraccap"))

ggplot(data = statsinc, aes(x = (time), y = value, col = variable))+ theme_bw() +
  labs(x ="Time (Days)", y = "Daily Incidence", color = "Population") + scale_y_continuous(limits = c(0,2000), expand = c(0,0)) + 
  scale_x_continuous(limits= c(test["TimeContactTracEqualInc"], test["TimeTargetReach"]+10), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_line(size = 1.02, stat = "identity")

ggplot(data = out1, aes(x = (time), y = I*N))+ theme_bw() +
  labs(x ="Time (Days)", y = "Total Infected", color = "Population") + scale_y_continuous(limits = c(0, 20000),expand = c(0,0)) + 
  scale_x_continuous(limits = c(test["TimeContactTracEqualInc"], test["TimeTargetReach"]+10), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_line(size = 1.02, stat = "identity")

