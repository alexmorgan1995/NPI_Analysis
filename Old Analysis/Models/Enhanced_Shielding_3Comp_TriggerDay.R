rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur, scaling1) {
  ifelse((time >= tstart1 & time <= tstart1+tdur),
         (0.6*(1/(GenTime(3.3,2.8)))), #Intervention
         ifelse((time >= tstart1+tdur & time <= 730), 
                (0.6*(1/(GenTime(3.3,2.8))))*scaling1, #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta1(seq(0, 720), 71, (24*7), 0.4))

beta2 <- function(time, tstart1, tdur, scaling2) {
  ifelse((time >= tstart1 & time <= tstart1+tdur),
         (0.6*(1/(GenTime(3.3,2.8)))), #Intervention
         ifelse((time >= tstart1+tdur & time <= 730), 
                (1.7*(1/(GenTime(3.3,2.8))))*scaling2, #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta2(seq(0, 720), 71, (24*7), 0.4))

beta3 <- function(time, tstart1, tdur, scaling3) {
  ifelse((time >= tstart1 & time <= tstart1+tdur),
         (0.6*(1/(GenTime(3.3,2.8)))), #Intervention
         ifelse((time >= tstart1+tdur & time <= 730), 
                (2.8*(1/(GenTime(3.3,2.8))))*scaling3, #After Intervention
                (1.7*(1/(GenTime(3.3,2.8)))) #Before Intervention
         )
  )
}

plot(beta3(seq(0, 720), 71, (6*7), 0.8))


#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

    if(time < tstart1) {
      dSv = - beta1(time,tstart1,tdur,scaling1)*(Iv+Ih+Inv)*Sv + zeta*Rv
      
      dSh = - beta2(time,tstart1,tdur,scaling2)*(Iv+Ih+Inv)*Sh + zeta*Rh
      
      dSnv = - beta3(time,tstart1,tdur,scaling3)*(Iv+Ih+Inv)*Snv  + zeta*Rnv
      
      dIv = beta1(time,tstart1,tdur,scaling1)*(Iv+Ih+Inv)*Sv - gamma*Iv
      
      dIh =  beta2(time,tstart1,tdur,scaling2)*(Iv+Ih+Inv)*Sh - gamma*Ih
      
      dInv = beta3(time,tstart1,tdur,scaling3)*(Iv+Ih+Inv)*Snv - gamma*Inv
      
      dRv = gamma*Iv - zeta*Rv
      
      dRh = gamma*Ih - zeta*Rh
      
      dRnv = gamma*Inv - zeta*Rnv
      
      dCv = beta1(time,tstart1,tdur,scaling1)*(Iv+Ih+Inv)*Sv
    }
    else{
      dSv = - beta1(time,tstart1,tdur,scaling1)*Iv*Sv - beta1(time,tstart1,tdur,scaling1)*Ih*Sv + zeta*Rv
      
      dSh = - beta2(time,tstart1,tdur,scaling2)*Ih*Sh - beta1(time,tstart1,tdur,scaling2)*Iv*Sh - 
        beta2(time,tstart1,tdur,scaling2)*Inv*Sh + zeta*Rh
      
      dSnv = - beta3(time,tstart1,tdur,scaling3)*Inv*Snv - beta2(time,tstart1,tdur,scaling3)*Ih*Snv + zeta*Rnv
      
      
      
      dIv = beta1(time,tstart1,tdur,scaling1)*Iv*Sv + beta1(time,tstart1,tdur,scaling1)*Ih*Sv - gamma*Iv
      
      dIh =  beta2(time,tstart1,tdur,scaling2)*Ih*Sh + beta1(time,tstart1,tdur,scaling2)*Iv*Sh + 
        beta2(time,tstart1,tdur,scaling2)*Inv*Sh - gamma*Ih
      
      dInv = beta3(time,tstart1,tdur,scaling3)*Inv*Snv + beta2(time,tstart1,tdur,scaling3)*Ih*Snv - gamma*Inv
      
      
  
      dRv = gamma*Iv - zeta*Rv
      
      dRh = gamma*Ih - zeta*Rh
      
      dRnv = gamma*Inv - zeta*Rnv
      
      dCv = beta1(time,tstart1,tdur,scaling1)*Iv*Sv + beta1(time,tstart1,tdur,scaling1)*Ih*Sv 
    }
    return(list(c(dSv, dSh, dSnv, dIv, dIh, dInv, dRv, dRh, dRnv, dCv)))
  })
}

#### Testing the Model Structure + Obtaining Specific Information from Model Runs ####
TriggerDay <- seq(0,100, by = 1)

Trigdaydata <- data.frame(matrix(nrow = length(TriggerDay), ncol = 4))

for(i in 1:length(TriggerDay)) {
  temp <- numeric(4)
  init <- c(Sv = 0.2, Sh = 0.2, Snv = 0.6-0.0001, Iv = 0, Ih = 0, Inv = 0.0001, Rv= 0, Rh = 0, Rnv = 0, Cv = 0)
  times <- seq(0,730,by = 1)
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            zeta = 1/365,
            tstart1 = TriggerDay[i], 
            tdur = 6*7,
            scaling1 = 0,
            scaling2 = 0,
            scaling3 = 0)
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$Iv <- out1$Iv/0.20; out1$Ih <- out1$Ih/0.2; out1$Inv <- out1$Inv/0.60
  out1$R <- out1$Rv + out1$Rh + out1$Rnv
  out1$Beta1 <- beta1(seq(0,730), TriggerDay[i], (6*7), parms[5]) #VARY THIS ASWELL
  out1$Beta2 <- beta2(seq(0,730), TriggerDay[i], (6*7), parms[6])
  out1$Beta3 <- beta3(seq(0,730), TriggerDay[i], (6*7), parms[7])
  temp[1] <- TriggerDay[i]
  temp[2] <- TriggerDay[i] + 6
  temp[3] <- out1$R[out1$time == TriggerDay[i] + 7]
  temp[4] <- out1$Iv[out1$time == TriggerDay[i] + 7]
  Trigdaydata[i,] <- temp
}

colnames(Trigdaydata) <- c("TriggerDay", "1Week", "Rt", "Iv")

ggplot(Trigdaydata, aes(x = TriggerDay, y = Rt)) + geom_line() + geom_hline(yintercept = 0.06, lty = 2) +
  labs(x ="Trigger Day", y = "Proportion Total Recovered - 1 Week After") + scale_y_continuous(expand = c(0,0)) +
  theme(legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 


#### Optimisation ####

start_time <- Sys.time()

beta1t <- seq(0,1, by = 0.01)
beta3t <- seq(0,1, by = 0.01)
betacomb <- expand.grid(beta1t, beta3t)
beta2t <- seq(0,1, by = 0.2)

betadata <- data.frame(matrix(nrow = 0, ncol = 6))

for(j in 1:length(beta2t)) { 
  tempdata <- data.frame(matrix(nrow = nrow(betacomb), ncol = 6))
  for(i in 1:nrow(betacomb)) { 
    temp <- numeric(6)
    init <- c(Sv = 0.2, Sh = 0.2, Snv = 0.6-0.0001, Iv = 0, Ih = 0, Inv = 0.0001, Rv= 0, Rh = 0, Rnv = 0, Cv = 0)
    times <- seq(0,436,by = 1)
    parms = c(gamma = 1/(GenTime(3.3,2.8)), 
              zeta = 1/365,
              tstart1 = 71, 
              tdur = 6*7,
              scaling1 = betacomb[i,1],
              scaling2 = (1.7*(1/(GenTime(3.3,2.8))) - ((1.7*(1/(GenTime(3.3,2.8))) - 0.6*(1/(GenTime(3.3,2.8))))*(1-beta2t[j])))/
                (1.7*(1/(GenTime(3.3,2.8)))),
              scaling3 = (2.8*(1/(GenTime(3.3,2.8))) - ((2.8*(1/(GenTime(3.3,2.8))) - 1.7*(1/(GenTime(3.3,2.8))))*(1-betacomb[i,2])))/
                (2.8*(1/(GenTime(3.3,2.8)))))
    out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
    out1$Sv <- out1$Sv/0.20; out1$Sh <- out1$Sh/0.2; out1$Snv <- out1$Snv/0.60
    out1$Iv <- out1$Iv/0.20; out1$Ih <- out1$Ih/0.2; out1$Inv <- out1$Inv/0.60
    out1$Rv <- out1$Rv/0.20; out1$Rh <- out1$Rh/0.2; out1$Rnv <- out1$Rnv/0.60
    out1$Cv <- out1$Cv/0.20
    out1$Beta1 <- beta1(seq(0,436), 71, (6*7), parms[5]) #VARY THIS ASWELL
    out1$Beta2 <- beta2(seq(0,436), 71, (6*7), parms[6])
    out1$Beta3 <- beta3(seq(0,436), 71, (6*7), parms[7])
    temp[1] <- out1$Beta1[out1$time == max(out1$time)]
    temp[2] <- out1$Beta2[out1$time == max(out1$time)]
    temp[3] <- out1$Beta3[out1$time == max(out1$time)]
    temp[4] <- out1$Cv[out1$time == max(out1$time)]
    temp[5] <- out1$Iv[which.max(out1$Iv)]
    temp[6] <- out1$Ih[which.max(out1$Ih)]
    tempdata[i,] <- temp
    print(paste0(signif(i/nrow(betacomb), digits = 2), signif(out1$Beta2[out1$time == max(out1$time)], digits = 3)))
  }
  betadata <- rbind.data.frame(betadata, tempdata)
}

colnames(betadata) <- c("Beta1", "Beta2","Beta3", "Cv", "Iv","Ih")

end_time <- Sys.time(); end_time - start_time


#HeatMap
#Beta2 = 0.07
betadattest <- betadata

ggplot(betadattest[1:10201,], aes(x = Beta1, y = Beta3, z = Cv)) + geom_contour_filled(breaks = seq(0,0.2,by = 0.01)) 

#Beta2 = 0.0957
ggplot(betadattest[10202:20402,], aes(x = Beta1, y = Beta3, z = Cv))+ geom_contour_filled(breaks = seq(0,0.2,by = 0.01)) 



#Beta2 = 0.1213591
ggplot(betadattest[20403:30603,], aes(x = Beta1, y = Beta3, z = Cv)) + geom_contour_filled(breaks = seq(0,0.2,by = 0.01)) 


#Beta2 = 0.1470312
ggplot(betadattest[30604:40804,], aes(x = Beta1, y = Beta3, z = Cv)) + geom_contour_filled(breaks = seq(0,0.2,by = 0.01)) 


#Beta2 = 0.1727033
ggplot(betadattest[40805:51005,], aes(x = Beta1, y = Beta3, z = Cv))  + geom_contour_filled(breaks = seq(0,0.2,by = 0.01)) 

#Beta2 = 0.1983755
ggplot(betadattest[51006:61206,], aes(x = Beta1, y = Beta3, z = Cv)) + geom_contour_filled(breaks = seq(0,0.2,by = 0.01)) 

#HeatMap
#Beta2 = 0.07
betadattest <- betadata

ggplot(betadattest[1:10201,], aes(x = Beta1, y = Beta3, z = Ih)) + geom_contour_filled(breaks = seq(0,0.03,by = 0.001)) 

#Beta2 = 0.0957
ggplot(betadattest[10202:20402,], aes(x = Beta1, y = Beta3, z = Ih))+ geom_contour_filled(breaks = seq(0,0.03,by = 0.001)) 



#Beta2 = 0.1213591
ggplot(betadattest[20403:30603,], aes(x = Beta1, y = Beta3, z = Ih)) + geom_contour_filled(breaks = seq(0,0.03,by = 0.001)) 


#Beta2 = 0.1470312
ggplot(betadattest[30604:40804,], aes(x = Beta1, y = Beta3, z = Ih)) + geom_contour_filled(breaks = seq(0,0.03,by = 0.001)) 


#Beta2 = 0.1727033
ggplot(betadattest[40805:51005,], aes(x = Beta1, y = Beta3, z = Ih))  + geom_contour_filled(breaks = seq(0,0.03,by = 0.001)) 

#Beta2 = 0.1983755
ggplot(betadattest[51006:61206,], aes(x = Beta1, y = Beta3, z = Ih)) + geom_contour_filled(breaks = seq(0,0.03,by = 0.001)) 
