rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSv = - beta1*Iv*Sv - beta1*Ih*Sv + zeta*Rv
    
    dSh = - beta2*Ih*Sh - beta2*Iv*Sh - beta2*Inv*Sh + zeta*Rh
    
    dSnv = - beta3*Inv*Snv - beta3*Ih*Snv + zeta*Rnv
    
    
    
    dIv = beta1*Iv*Sv + beta1*Ih*Sv - gamma*Iv
    
    dIh =  beta2*Ih*Sh + beta2*Iv*Sh +  beta2*Inv*Sh - gamma*Ih
      
    dInv = beta3*Inv*Snv + beta3*Ih*Snv - gamma*Inv
    
    
    
    dRv = gamma*Iv - zeta*Rv
    
    dRh = gamma*Ih - zeta*Rh
    
    dRnv = gamma*Inv - zeta*Rnv
    
    return(list(c(dSv, dSh, dSnv, dIv, dIh, dInv, dRv, dRh, dRnv)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####

init <- c(Sv = 0.1-0.0001, Sh = 0.7-0.0001, Snv = 0.1-0.0001, Iv = 0.0001, Ih =  0.0001, Inv = 0.0001, Rv= 0, Rh = 0, Rnv = 0)
times <- seq(0,730,by = 1)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          beta1 = 0.15,
          beta2 = 0.15,
          beta3 = 0.15)


out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))

out1$Sv <- out1$Sv/0.10; out1$Sh <- out1$Sh/0.7; out1$Snv <- out1$Snv/0.1
out1$Iv <- out1$Iv/0.10; out1$Ih <- out1$Ih/0.7; out1$Inv <- out1$Inv/0.10
out1$Rv <- out1$Rv/0.10; out1$Rh <- out1$Rh/0.7; out1$Rnv <- out1$Rnv/0.10


colnames(out1) <- c("Time", "Suscv", "Susch", "Suscnv", "Infected_Iv", "Infected_Ih", "Infected_Inv", 
                    "Recovv", "Recovh", "Recovnv")


statsSUSC <- melt(out1, id.vars = c("Time"), measure.vars = c("Suscv", "Susch", "Suscnv"))

ggplot(data = statsSUSC, aes(x = (Time), y = value, col = variable)) + geom_line(size = 1.02, stat = "identity")

