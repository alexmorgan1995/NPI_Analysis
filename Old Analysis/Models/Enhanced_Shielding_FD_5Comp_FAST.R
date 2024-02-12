rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("RColorBrewer"); library("sensitivity");library("fast")

#### Model Functions ####

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

gamma <- 1/(GenTime(3.3,2.8))

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur, scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (0.8*(gamma))*scaling
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta1(seq(0,730), 71, (6*7), 0.5), ylim = c(0,0.5))

beta2 <- function(time, tstart1, tdur, scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (2.8*(gamma) - ((2.8*(gamma) - 0.9*(gamma))*scaling))
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta2(seq(0,730), 71, (6*7), 0.5))

beta3 <- function(time, tstart1, tdur, scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (2.8*(gamma) - (2.8*(gamma) - 1.7*(gamma))*scaling)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta3(seq(0,730), 71, (6*7), 0.5))


beta4 <- function(time,tstart1,tdur,scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.8*(gamma) *scaling
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta4(seq(0,730), 71, (6*7), 0.5))

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- beta1(time,tstart1,tdur,scaling)
    beta2 <- beta2(time,tstart1,tdur,scaling)
    beta3 <- beta3(time,tstart1,tdur,scaling)
    beta4 <- beta4(time,tstart1,tdur,scaling)
    
    dSv = - beta1*Iv*Sv*scalingAA - 
      beta1*Is*Sv*scalingBB - 
      beta4*Ir1*Sv*scalingCC - 
      beta4*Ir2*Sv*scalingCC - 
      beta4*Ir3*Sv*scalingCC + zeta*Rv
    
    dSs = - beta1*Iv*Ss*scalingDD - 
      beta1*Is*Ss*scalingEE -
      beta2*Ir1*Ss*scalingFF -
      beta2*Ir2*Ss*scalingFF -
      beta2*Ir3*Ss*scalingFF + zeta*Rs

    dSr1 = - beta4*Iv*Sr1*scalingGG - 
      beta2*Is*Sr1*scalingHH - 
      beta3*Ir1*Sr1*scalingII - 
      beta3*Ir2*Sr1*scalingII - 
      beta3*Ir3*Sr1*scalingII + zeta*Rr1
      
    dSr2 = - beta4*Iv*Sr2*scalingGG -
      beta2*Is*Sr2*scalingHH -
      beta3*Ir1*Sr2*scalingII -
      beta3*Ir2*Sr2*scalingII -
      beta3*Ir3*Sr2*scalingII + zeta*Rr2
      
    dSr3 = - beta4*Iv*Sr3*scalingGG -
      beta2*Is*Sr3*scalingHH -
      beta3*Ir1*Sr3*scalingII -
      beta3*Ir2*Sr3*scalingII -
      beta3*Ir3*Sr3*scalingII + zeta*Rr3
    
    
    dIv = beta1*Iv*Sv*scalingAA + 
      beta1*Is*Sv*scalingBB + 
      beta4*Ir1*Sv*scalingCC + 
      beta4*Ir2*Sv*scalingCC + 
      beta4*Ir3*Sv*scalingCC - gamma*Iv
    
    dIs = beta1*Iv*Ss*scalingDD + 
      beta1*Is*Ss*scalingEE +
      beta2*Ir1*Ss*scalingFF +
      beta2*Ir2*Ss*scalingFF +
      beta2*Ir3*Ss*scalingFF - gamma*Is
    
    dIr1 = beta4*Iv*Sr1*scalingGG + 
      beta2*Is*Sr1*scalingHH + 
      beta3*Ir1*Sr1*scalingII + 
      beta3*Ir2*Sr1*scalingII + 
      beta3*Ir3*Sr1*scalingII - gamma*Ir1
    
    dIr2 = beta4*Iv*Sr2*scalingGG +
      beta2*Is*Sr2*scalingHH +
      beta3*Ir1*Sr2*scalingII +
      beta3*Ir2*Sr2*scalingII +
      beta3*Ir3*Sr2*scalingII - gamma*Ir2
    
    dIr3 = beta4*Iv*Sr3*scalingGG +
      beta2*Is*Sr3*scalingHH +
      beta3*Ir1*Sr3*scalingII +
      beta3*Ir2*Sr3*scalingII +
      beta3*Ir3*Sr3*scalingII - gamma*Ir3
    
    
    dRv = gamma*Iv - zeta*Rv
    dRs = gamma*Is - zeta*Rs
    dRr1 = gamma*Ir1 - zeta*Rr1
    dRr2 = gamma*Ir2 - zeta*Rr2
    dRr3 = gamma*Ir3 - zeta*Rr3
    
    dCIv = beta1*Iv*Sv*scalingAA + 
      beta1*Is*Sv*scalingBB + 
      beta4*Ir1*Sv*scalingCC + 
      beta4*Ir2*Sv*scalingCC + 
      beta4*Ir3*Sv*scalingCC
    
    return(list(c(dSv, dSs, dSr1, dSr2, dSr3,
                  dIv, dIs, dIr1, dIr2, dIr3,
                  dRv, dRs, dRr1, dRr2, dRr3, dCIv)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####

start_time <- Sys.time()

parms = fast_parameters(minimum = c(0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75), 
                        
                        maximum = c(1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25), 
                        factor=9, names = c("beta_AA", "beta_BB" ,"beta_CC", "beta_DD", 
                                            "beta_EE", "beta_FF", "beta_GG", "beta_HH",
                                            "beta_II"))

init <- c(Sv = 0.2 - 0.0001*0.2, Ss = 0.2 - 0.0001*0.2, 
          Sr1 = 0.2 - 0.0001*0.2, Sr2 = 0.2 - 0.0001*0.2, Sr3 = 0.2 - 0.0001*0.2,
          Iv = 0.0001*0.2, Is = 0.0001*0.2, Ir1 = 0.0001*0.2, Ir2 = 0.0001*0.2, Ir3 = 0.0001*0.2,   
          Rv= 0, Rs = 0, Rr1 = 0, Rr2 = 0, Rr3 = 0,
          CIv = 0)
times <- seq(0, 478, by = 1)

output <- data.frame(matrix(ncol = 6, nrow = nrow(parms)))
colnames(output) <- c("CumIvYr","TimeSecPeak","HeightSecPeak","TimeFirPeak","HeightFirPeak","HigherFirPeak")

for (i in 1:nrow(parms)) {
  temp <- numeric(6)
  
  parms1 = c(gamma = 1/(GenTime(3.3,2.8)), 
             zeta = 1/365,
             tstart1 = 71, 
             tdur = 6*7,
             scaling = 0.5,
             scalingAA = parms$beta_AA[i],
             scalingBB = parms$beta_BB[i],
             scalingCC = parms$beta_CC[i],
             scalingDD = parms$beta_DD[i],
             scalingEE = parms$beta_EE[i],
             scalingFF = parms$beta_FF[i],
             scalingGG = parms$beta_GG[i],
             scalingHH = parms$beta_HH[i],
             scalingII = parms$beta_II[i])
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms1))
  out1$Iv <- out1$Iv/0.20
  out1$CIv <- out1$CIv/0.20
  temp[1] <- out1$CIv[max(out1$time)]
  temp[2] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][2]
  temp[3] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][2]
  temp[4] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][1]
  temp[5] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][1]
  temp[6] <- ifelse((out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)][2] > out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)][1]), 1, 0)
  output[i,] <- temp
  print(i/nrow(parms))
}

end_time <- Sys.time(); end_time - start_time

#### Plotting Of Sensitivity Analysis ##### 

#ICombH
sensit1 <- NULL; df.equilibrium <- NULL
sensit1 <- output$HeightSecPeak
sens1 <- sensitivity(x=sensit1, numberf=9, make.plot=T, names = c("beta_AA", "beta_BB" ,"beta_CC", "beta_DD", 
                                                                   "beta_EE", "beta_FF", "beta_GG", "beta_HH",
                                                                   "beta_II"))
df.equilibrium <- data.frame(parameter=rbind("beta_AA", "beta_BB" ,"beta_CC", "beta_DD", 
                                             "beta_EE", "beta_FF", "beta_GG", "beta_HH",
                                             "beta_II"), value=sens1)
ggplot(df.equilibrium, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
