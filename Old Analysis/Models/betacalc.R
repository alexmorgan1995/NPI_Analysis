rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

(1.7*(1/(GenTime(3.3,2.8))) - ((1.7*(1/(GenTime(3.3,2.8))) - 0.6*(1/(GenTime(3.3,2.8))))*0.5))
(2.8*(1/(GenTime(3.3,2.8))) - ((2.8*(1/(GenTime(3.3,2.8))) - 1.7*(1/(GenTime(3.3,2.8))))*0.5))

beta1 <- (0.5*(1/(GenTime(3.3,2.8))))
beta2 <- (0.25*(1/(GenTime(3.3,2.8))))
beta3 <- (1.7*(1/(GenTime(3.3,2.8))))

bv <- 0.6/0.2
br <- 0.2/0.6
cb <- 0.5
#vv
beta1 + beta1*bv*0.5*(1-cb)

#vr + rv
beta4*cb

#rr
beta3 + beta3*br*(1-cb)
(0.25*(1/(GenTime(3.3,2.8))))
