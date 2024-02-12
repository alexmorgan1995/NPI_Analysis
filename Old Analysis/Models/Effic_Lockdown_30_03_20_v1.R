rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot");library("reshape2"); library("dplyr"); library("RColorBrewer")

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### SIR Model ####
#### Exploring Week 3 and 12 Interventions - SIR Model - Will need to tweak code for each model output#### 

#Function to model intervention - currently set at baseline - added additional functionality to it
betastatdecrease <- function(time, int_timestart, int_timeend, effic_lock) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*effic_lock,
         (1.5*(1/(GenTime(4.6,2.4)))))
}

#Function for SIR ODEs with cumulative infection compartment
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend, effic_lock)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend, effic_lock)*S*I - mu*I
    dR = mu*I 
    dC = betastatdecrease(time,int_timestart,int_timeend, effic_lock)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

#Initial Conditions and Times
init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,162,by = 1)
duration <- 3 # Duration of the intervention
effic <- c((1.5-(1.2*1))/1.5,
           (1.5-(1.2*0.75))/1.5,
           (1.5-(1.2*0.5))/1.5,
           (1.5-(1.2*0.25))/1.5,
           (1.5-(1.2*0))/1.5)

#Run the ODE for 3 and 12 Week Interventions
stats1 <- data.frame(matrix(nrow = 0, ncol = 6))
EfficPerc <- c("100%","75%","50%","25%","0%")

for(i in 1:length(effic)) {
  temp <- data.frame(matrix(nrow = length(seq(0,162)), ncol = 6))
  parms = c(mu = 1/(GenTime(4.6,2.4)), 
            int_timestart = 100, 
            int_timeend = 100 + (duration*7),
            effic_lock = effic[i])
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- EfficPerc[i]
  temp[,2] <- c(seq(-100,0), seq(1,162-100))
  temp[,3] <- out1$I 
  temp[,4] <- out1$R 
  temp[,5] <- out1$C 
  temp[,6] <- betastatdecrease(times, 
                               as.numeric(parms[2]), 
                               as.numeric(parms[3]),
                               as.numeric(parms[4]))
  stats1 <- rbind(stats1, temp)
}



colnames(stats1) <- c("Lockdown_Efficacy", "Time", "Infec", "Recov", "Cum", "Beta")

stats1 <- melt(stats1, id.vars = c("Lockdown_Efficacy", "Time"), measure.vars = c("Infec", "Recov", "Cum", "Beta"))
stats1$Lockdown_Efficacy <- factor(stats1$Lockdown_Efficacy, levels = unique(stats1$Lockdown_Efficacy))

#Plotting 
p11 <- ggplot(data = stats1[which(stats1$variable == "Infec"),], aes(x = (Time), y = value, col = Lockdown_Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.075) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0, col = "black", size = 1, lty = 2)  +
  scale_color_manual(values=c("black", "#1F78B4", "#B2DF8A", "#33A02C", "darkred"))

p12 <- ggplot(data = stats1[which(stats1$variable == "Recov"),], aes(x = (Time), y = value, col = Lockdown_Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recoved") + scale_y_continuous(limits = c(0,0.6) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_color_manual(values=wes_palette(n=5, name="Darjeeling1")) +
  scale_x_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0, col = "black", size = 1, lty = 2) +
  scale_color_manual(values=c("black", "#1F78B4", "#B2DF8A", "#33A02C",  "darkred"))

p21 <- ggplot(data = stats1[which(stats1$variable == "Beta"),], aes(x = (Time), y = value, col = Lockdown_Efficacy)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time Relative to Lockdown Date (Days)", y = expression(beta)) + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom" , axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.title =element_text(size=14), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + 
  guides(size = guide_legend(title.position="top", title.hjust = 0.5, fill=guide_legend(nrow=2,byrow=TRUE)))  +
  scale_x_continuous(expand = c(0, 0)) + labs(color='Level of Compliance') +
  scale_color_manual(values=c("black", "#1F78B4", "#B2DF8A", "#33A02C",  "darkred"))

plot_grid(p11, p12, p21, align = "v", nrow = 3, rel_heights = c(0.3, 0.3, 0.25))
