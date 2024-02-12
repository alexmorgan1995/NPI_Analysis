rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot"); library("reshape2"); library("dplyr"); library("RColorBrewer")

#### Model Functions ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model intervention - currently set at baseline - added additional functionality to it
beta1 <- function(time, tstart1, tdur1,  beta1_2) {
  betalin <- approxfun(x=c(tstart1+tdur1, tstart1+tdur1+(12*7)),y = c((0.5*(1/(GenTime(3.3,2.8)))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 2
         (0.5*(1/(GenTime(3.3,2.8)))),
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur1+(12*7) & time <= tstart1+tdur1+(24*7)),
                       beta1_2,
                       ifelse((time >= tstart1+tdur1+(24*7) & time <= 730),
                              beta1_2, 
                              (1.7*(1/(GenTime(3.3,2.8))))
                       )
                )
         )
  )
}

beta2 <- function(time, tstart1, tdur1,  beta1_2) {
  betalin <- approxfun(x=c(tstart1+tdur1, tstart1+tdur1+(12*7)),y = c((0.6*(1/(GenTime(3.3,2.8)))), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur1), #Phase 2
         (0.6*(1/(GenTime(3.3,2.8)))),
         ifelse((time >= tstart1+tdur1 & time <= tstart1+tdur1+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur1+(12*7) & time <= tstart1+tdur1+(24*7)),
                       beta1_2,
                       ifelse((time >= tstart1+tdur1+(24*7) & time <= 730),
                              beta1_2, 
                              (1.7*(1/(GenTime(3.3,2.8))))
                       )
                )
         )
  )
}

plot(beta2(seq(0,730), 71, (6*7), 0.1))

beta4 <- function(time, tstart1) {
  ifelse((time >= tstart1 & time <= 730), #Phase 2
         0, 
         1.7*(1/(GenTime(3.3,2.8)))
  )
}

plot(beta4(seq(0,730), 71))

temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 5))

temp[,1] <- beta1(seq(0,730), 71, (6*7), (0.25*(1/(GenTime(3.3,2.8)))))
temp[,2] <- beta2(seq(0,730), 71, (6*7), (1.15*(1/(GenTime(3.3,2.8)))))
temp[,3] <- beta2(seq(0,730), 71, (6*7), (2.25*(1/(GenTime(3.3,2.8)))))
temp[,4] <- beta4(seq(0,730), 71)
temp[,5] <- seq(0,730)

ggplot(data = temp, aes(x = (X5), y = X1)) + geom_line(size = 1.02, stat = "identity", col = "darkred") +
  labs(x ="Time (Days)", y = "??1") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

ggplot(data = temp, aes(x = (X5), y = X2)) + geom_line(size = 1.02, stat = "identity", col = "darkblue") +
  labs(x ="Time (Days)", y = "??2") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0))

ggplot(data = temp, aes(x = (X5), y = X3)) + geom_line(size = 1.02, stat = "identity", col = "darkgreen") +
  labs(x ="Time (Days)", y = "??3") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 

ggplot(data = temp, aes(x = (X5), y = X4)) + geom_line(size = 1.02, stat = "identity", col = "darkorange") +
  labs(x ="Time (Days)", y = "??4") + scale_y_continuous(limits = c(0,0.4),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) 