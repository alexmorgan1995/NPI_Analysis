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
    
    dS = - beta*init["I"]*init["S"] + parms["zeta"]*init["R"]
    dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
    dR = parms["gamma"]*init["I"] - parms["zeta"]*init["R"] 
    
    inc_tim_inc <- c(time, beta*init["I"]*init["S"])
    inc_cont <<- rbind(inc_cont, inc_tim_inc) # Tracking the daily incidence
    inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
    
    return(list(c(dS, dI, dR)))
}

#Linear Increase in contact tracing

lin_cont <- function(x,m,comp) {
  y = m*(x-comp+1) + 0
  return(y)
}

#SIRS Model To Identify where contact tracing = incidence

SIRS2 <- function(time, init, parms) {
    timecomp = timecomp
    
    if(time >= timecomp) {
      
      beta <- beta(time, parms["tstart1"], parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
      linear <<- lin_cont(time, parms["ramp_rate"], timecomp)
      
      dS = - beta*init["I"]*init["S"] + parms["zeta"]*init["R"]
      dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
      dR = parms["gamma"]*init["I"] - parms["zeta"]*init["R"] 
      
      incidence <- beta*init["I"]*init["S"]*parms["N"]
      
      inctime <- c(time, incidence, linear, beta, linear/ (incidence))
      inc_cont <<- rbind(inc_cont, inctime)
      inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
      
    } else { #if inc < 1000
      
      beta <- beta(time, parms["tstart1"], parms["R0"])
      
      dS = - beta*init["I"]*init["S"] + parms["zeta"]*init["R"]
      dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
      dR = parms["gamma"]*init["I"] - parms["zeta"]*init["R"]
      
      incidence <- beta*init["I"]*init["S"]*parms["N"]
      
      inctime <- c(time, incidence, 0, beta, 0/ (incidence))
      inc_cont <<- rbind(inc_cont, inctime)
      inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
    }
    
    return(list(c(dS, dI, dR)))
}

#SIRS Model used for Final Model Output

SIRS3 <- function(time, init, parms) {

    timecomp = timecomp
    timeequal = timeequal
    
    if(time >= timecomp) {
      
      if(time >= timeequal) { # if time is equal
        
        beta <- beta(time,parms["tstart1"], parms["R0"])*(1-(1*parms["efficacy"]))
        
        dS = - beta*init["I"]*init["S"] + parms["zeta"]*init["R"]
        dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
        dR = parms["gamma"]*init["I"] - parms["zeta"]*init["R"] 
        
        incidence <- beta*init["I"]*init["S"]*parms["N"]
        
        inctime <- c(time, incidence, incidence, beta, 1)
        inc_cont <<- rbind(inc_cont, inctime)
        inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
        
      } else{ #if inc is over 1000
        
        beta <- beta(time,parms["tstart1"], parms["R0"]) * (1-(inc_cont[nrow(inc_cont),3]/(inc_cont[nrow(inc_cont),2]))*parms["efficacy"])
        linear <- lin_cont(time, parms["ramp_rate"], timecomp)
        
        dS = - beta*init["I"]*init["S"] + parms["zeta"]*init["R"]
        dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
        dR = parms["gamma"]*init["I"] - parms["zeta"]*init["R"] 
        
        incidence <- beta*init["I"]*init["S"]*parms["N"]
        
        inctime <- c(time, incidence, linear, beta, 
                     linear/ (incidence))
        
        inc_cont <<- rbind(inc_cont, inctime)
        inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
      }
      
    } else { #if inc < 1000
      
      beta <- beta(time, parms["tstart1"], parms["R0"])
      
      dS = - beta*init["I"]*init["S"] + parms["zeta"]*init["R"]
      dI = beta*init["I"]*init["S"] - parms["gamma"]*init["I"]
      dR = parms["gamma"]*init["I"] - parms["zeta"]*init["R"] 
      
      incidence <- beta*init["I"]*init["S"]*parms["N"]
      
      inctime <<- c(time, incidence , 0, beta, 0)
      inc_cont <<- rbind(inc_cont, inctime)
      inc_cont <<- inc_cont[!duplicated(inc_cont$X0), ]
      
    }
    return(list(c(dS, dI, dR)))
}

#### Big Function ####

N <- 5.5*10^6
init <- c(S = (N-(N*0.0001))/N, I = (N*0.0001)/N, R= 0)

contact_trac_run <- function(efficacy, ramp_rate, R0, target, initial) {
  
  #Parameters
  N <- 5.5*10^6
  init <- c(S = (N-(N*0.0001))/N, I = (N*0.0001)/N, R= 0)
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            zeta = 1/365,
            tstart1 = 71, 
            N = 5.5*10^6,
            efficacy = efficacy,
            ramp_rate = ramp_rate, 
            R0 = R0,
            target = target)
  
  #Baseline 
  times <- seq(0, 400, by = 1)
  inc_cont <<- data.frame(0, 0)
  out1 <- data.frame(rk(y = init, func = SIRS1, times = times, parms = parms, method = "rk4"))
  timecomp <<- inc_cont[,1][which(abs((inc_cont[,2]*N - initial) - 0) == min(abs((inc_cont[,2]*N - initial) - 0)) & inc_cont[,1] > parms["tstart1"])]
  
  #Intersect
  inc_cont <<- data.frame(0, 0, 0, 0, 0)
  times <- seq(0, 400, by = 1)
  out1 <- data.frame(rk(y = init, func = SIRS2, times = times, parms = parms, method = "rk4"))
  
  t <- inc_cont[,1][which(inc_cont[,5] < 1 & inc_cont[,5] > 0.5)]
  t1 <- split(t, cumsum(c(1, diff(t) != 0.5)))[[1]]
  timeequal <<- t1[length(t1)]
  
  #Actually Run the Model 
  inc_cont <<- data.frame(0, 0, 0, 0, 0)
  times <- seq(0, 400, by = 1)
  out1 <<- data.frame(rk(y = init, func = SIRS3, times = times, parms = parms, method = "rk4"))
  
  inc_cont[1,1] <<- 0; inc_cont[1,2] <<- 0
  colnames(inc_cont) <<- c("time", "incidence", "contacttraccap", "beta", "con/inc rat")
  
  time10 <- inc_cont[,1][which(inc_cont[,2] < target)][2]
  timesince <- inc_cont[,1][which(inc_cont[,2] < target)][2] - timecomp
  return(c("Time1000" = timeequal, "TimeContactTracEqualInc" = timecomp,"TimeTargetReach" = time10, "RelTimetoTarget" = timesince))
}

#Specify the efficacy (as a fraction), ramp_rate, lockdown R0 and target incidence
#inc_cont will always output during the run - view for a more detailed overview of the model run.

test <- contact_trac_run(0.8, 50, 0.9, 10, 1000)

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


#### Heat Maps #### 

heatmap <- expand.grid(seq(1,50),seq(0.7,1, by = 0.05))

heatdata <- data.frame(matrix(nrow = nrow(heatmap), ncol = 6))

for (i in 1:nrow(heatmap)) {
  heatdata[i,] <- c(heatmap[i,1],heatmap[i,2],contact_trac_run(0.8, heatmap[i,1], heatmap[i,2], 10))
  print((i/nrow(heatmap)*100))
}

colnames(heatdata) <- c("ContRampUp", "Lock_R0", "Time1000", "TimeContactTracEqualInc", "TimeTargetReach", "RelTimetoTarget")

ggplot(heatdata, aes(ContRampUp, Lock_R0, z = RelTimetoTarget)) + geom_contour_filled()
