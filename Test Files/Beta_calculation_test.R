library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

combbetamult <- function(scen, time, tstart1, t_dur1, tstart2, t_dur2, R0Dec1, R0Dec2) {
  gamma <- 1/GenTime(3.3, 2.8)
  if(scen == 0) {
    output <-  1.7*gamma
  }
  if(scen == 1) {
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse((time >= (tstart1 + t_dur1 + tstart2)),
                            R0Dec2*gamma,
                            R0Dec1*gamma),
                     1.7*gamma)
  }
  if(scen == 2) {
    R0lin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(R0Dec1, 1.7), method="linear", rule =2)
    R0lin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(R0Dec2, 1.7), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0lin1(time)*gamma,
                            R0lin2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 3) {
    R0lin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(1.7, R0Dec1), method="linear", rule =2)
    R0lin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(1.7, R0Dec2), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0lin1(time)*gamma,
                            R0lin2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 4) {
    R0incdec1 <-approxfun(x=c(tstart1, tstart1+(t_dur1/2), tstart1+t_dur1), y= c(1.7, R0Dec1, 1.7), method="linear", rule =2)
    R0incdec2 <-approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 +(t_dur2/2)), (tstart1 + t_dur1 + tstart2 +t_dur2)), y= c(1.7, R0Dec2, 1.7), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0incdec1(time)*gamma,
                            R0incdec2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart1 & time <= tstart1+(t_dur1*0.16667)) | (time >= tstart1+(t_dur1*0.3333) & time <= tstart1+(t_dur1*0.5)) |
                       (time >= tstart1+(t_dur1*0.6667) & time <= tstart1+(t_dur1*0.83333)) |
                       (time >= tstart1 + t_dur1 + tstart2 & time <= tstart1 + t_dur1 + tstart2+(t_dur2*0.16667)) | 
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.3333) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.5)) |
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.6667) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.83333)),
                     ifelse((time >= tstart1 + t_dur1 + tstart2),
                            R0Dec2*gamma,
                            R0Dec1*gamma),
                     1.7*gamma)
  }
  return(output)
}

plot(seq(0,365),combbetamult(5, seq(0,365), 71, 12*7, 20, 12*7, 0.8, 0.6))
