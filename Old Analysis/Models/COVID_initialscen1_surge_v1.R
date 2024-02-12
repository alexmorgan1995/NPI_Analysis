rm(list=ls())
library("deSolve"); library("ggplot2"); library("cowplot");library("reshape2"); library("dplyr")

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#R code sections with "Explore" requires the user to manually select parameters to correspond with the desired scenario
#To replicate the results in the technical document - the lockdown date must change
#THis requires code alteration to the for loop and the formatting data code sections
#and the R0 must also change in the beta alteration function

#### SIR Model ####
#### Exploring Week 3 and 12 Interventions - SIR Model - Will need to tweak code for each model output#### 

#Function to model intervention - currently set at baseline
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*0.4,
         (1.5*(1/(GenTime(4.6,2.4)))))
}

#Function for SIR ODEs with cumulative infection compartment
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - mu*I
    dR = mu*I 
    dC = betastatdecrease(time,int_timestart,int_timeend)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

#Initial Conditions and Times
init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,730,by = 1)
duration <- c(3,12) #3 and 12 Week Intervention duration

#Run the ODE for 3 and 12 Week Interventions
stats1 <- data.frame(matrix(nrow = 0, ncol = 6))
for(i in 1:length(duration)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 6))
  parms = c(mu = 1/(GenTime(4.6,2.4)), 
            int_timestart = 100, 
            int_timeend = 100 + (duration[i]*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- paste0(as.character(duration[i]), " Weeks")
  temp[,2] <- out1$time 
  temp[,3] <- out1$I 
  temp[,4] <- out1$R 
  temp[,5] <- out1$C 
  temp[,6] <- betastatdecrease(times, 
                               as.numeric(parms[2]), 
                               as.numeric(parms[3]))
  stats1 <- rbind(stats1, temp)
}

colnames(stats1) <- c("WeekDuration", "Time", "Infec", "Recov","Cum", "Beta")

#Formatting the Data - Have shifted the x axis so that 0 is the lockdown initiation date
combdatainf <- rbind(data.frame("WeekDuration" = "3 Weeks", "Time" = c(seq(-100,0), seq(1,730-100)), "Value" = stats1$Infec[stats1$WeekDuration == "3 Weeks"]),
                  data.frame("WeekDuration" = "12 Weeks", "Time" = c(seq(-100,0), seq(1,730-100)), "Value" = stats1$Infec[stats1$WeekDuration == "12 Weeks"]))
combdatainf$WeekDuration <- factor(combdatainf$WeekDuration, levels = unique(combdatainf$WeekDuration))

combdatarec <- rbind(data.frame("WeekDuration" = "3 Weeks", "Time" = c(seq(-100,0), seq(1,730-100)), "Value" = stats1$Recov[stats1$WeekDuration == "3 Weeks"]),
                   data.frame("WeekDuration" = "12 Weeks", "Time" = c(seq(-100,0), seq(1,730-100)), "Value" = stats1$Recov[stats1$WeekDuration == "12 Weeks"]))
combdatarec$WeekDuration <- factor(combdatarec$WeekDuration, levels = unique(combdatarec$WeekDuration))

combdatabeta <- rbind(data.frame("WeekDuration" = "3 Week Intervention", "Time" = c(seq(-100,0), seq(1,730-100)), "Value" = stats1$Beta[stats1$WeekDuration == "3 Weeks"]),
                      data.frame("WeekDuration" = "12 Week Intervention", "Time" = c(seq(-100,0), seq(1,730-100)), "Value" = stats1$Beta[stats1$WeekDuration == "12 Weeks"]))
combdatabeta$WeekDuration <- factor(combdatabeta$WeekDuration, levels = unique(combdatabeta$WeekDuration))

#Quick Analysis of Highest Peak
surge <- data.frame(matrix(nrow = 0, ncol = 5))

for(i in 1:length(duration)) { 
  test <- stats1[stats1$WeekDuration == as.character(paste0(duration[i], " Weeks")),][which.max(stats1$Infec[stats1$WeekDuration == as.character(paste0(duration[i], " Weeks"))]),]
  surge <- rbind(surge, test)
  print(test)
}

#Plotting 
p11 <- ggplot(data = combdatainf, aes(x = (Time), y = Value, col = WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0, col = "black", size = 1, lty = 2) +
  scale_color_manual(values=c("darkred","royalblue"))

p12 <- ggplot(data = combdatarec, aes(x = (Time), y = Value, col = WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recoved") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0, col = "black", size = 1, lty = 2)+
  scale_color_manual(values=c("darkred", "royalblue"))

p21 <- ggplot(data = combdatabeta, aes(x = (Time), y = Value, col = WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time Relative to Lockdown Date (Days)", y = "??") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  scale_linetype_manual(values=c("solid", "dotdash", "solid", "dotdash")) + 
  scale_color_manual(values=c("darkred", "royalblue"))

plot_grid(p11, p12, p21, align = "v", nrow = 3, rel_heights = c(0.3, 0.3, 0.25))

#### Exploring Weeks 1-12 Intervention Duration ####

#Function to model intervention - currently set at baseline
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*0.4,
         (1.5*(1/(GenTime(4.6,2.4)))))
}

#Initial Conditions
init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,730,by = 1)
duration <- seq(1,12) #Now do the analysis for all 1-12 week scenarios

#Run the ODEs for ALL 1-12 week intervention scenarios
stats1 <- data.frame(matrix(nrow = 0, ncol = 6))

for(i in 1:length(duration)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 6))
  parms = c(mu = 1/(GenTime(4.6,2.4)), 
            int_timestart = 100, 
            int_timeend = 100 + (duration[i]*7))
  out1 <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  temp[,1] <- paste0(as.character(duration[i]), " Weeks")
  temp[,2] <- c(seq(-100,0), seq(1,730-100)) #Have shifted the x axis so that 0 is the lockdown initiation date
  temp[,3] <- out1$I 
  temp[,4] <- out1$R 
  temp[,5] <- out1$C 
  temp[,6] <- betastatdecrease(times, 
                               as.numeric(parms[2]), 
                               as.numeric(parms[3]))
  stats1 <- rbind(stats1, temp)
}

colnames(stats1) <- c("WeekDuration", "Time", "Infec", "Recov","Cum", "Beta")

#Format the data for ggplot
stats1 <- melt(stats1, id.vars = c("WeekDuration", "Time"), measure.vars = c("Infec", "Recov", "Cum", "Beta"))
stats1$WeekDuration <- factor(stats1$WeekDuration, levels = unique(stats1$WeekDuration))

#Plotting - Have shifted the x axis so that 0 is the lockdown initiation date
p11 <- ggplot(data = stats1[which(stats1$variable == "Infec"),], aes(x = Time, y = value, col = WeekDuration, alpha=WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0, col = "black", size = 1, lty = 2) +
  scale_alpha_manual(values = c(0.4, 0.4, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1)) +
  scale_color_manual(values=c("#A6CEE3", "#1F78B4","darkred", "#B2DF8A", "#33A02C", "#FB9A99",
                              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A","royalblue"))

p12 <- ggplot(data = stats1[which(stats1$variable == "Recov"),], aes(x = (Time), y = value, col = WeekDuration, alpha=WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Recoved") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0, col = "black", size = 1, lty = 2)  +
  scale_alpha_manual(values = c(0.4, 0.4, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1)) +
  scale_color_manual(values=c("#A6CEE3", "#1F78B4","darkred", "#B2DF8A", "#33A02C", "#FB9A99",
                              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A","royalblue"))

p21 <- ggplot(data = stats1[which(stats1$variable == "Beta"),], aes(x = (Time), y = value, col = WeekDuration, alpha=WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time Relative to Lockdown Date (Days)", y = "??") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  +
  scale_alpha_manual(values = c(0.4, 0.4, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1)) +
  scale_color_manual(values=c("#A6CEE3", "#1F78B4","darkred", "#B2DF8A", "#33A02C", "#FB9A99",
                              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A","royalblue"))

plot_grid(p11, p12, p21, align = "v", nrow = 3, rel_heights = c(0.3, 0.3, 0.25))

#### SIS Model ####
#### Exploring Weeks 1-12 Intervention Duration - SIS Model ####

#Function for the SDM intervention
betastatdecrease <- function(time, int_timestart, int_timeend) {
  ifelse((time >= int_timestart & time <= int_timeend),
         (1.5*(1/(GenTime(4.6,2.4))))*0.4,
         (1.5*(1/(GenTime(4.6,2.4)))))
}

#SIS Model ODEs
SIS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betastatdecrease(time,int_timestart,int_timeend)*S*I + mu*I 
    dI = betastatdecrease(time,int_timestart,int_timeend)*S*I - mu*I
    dC = betastatdecrease(time,int_timestart,int_timeend)*S*I
    return(list(c(dS,dI,dC)))
  })
}

#Initial Conditions
init <- c(S = 0.9999, I = 0.0001, C = 0)
times <- seq(0,730,by = 1)
duration <- seq(1,12)

#Run the ODEs for all 12 intervention durations
#I have set lockdown initiation date to 98 to correspond with I(t) = ~0.0182
#This means that I(t) at lockdown initiation is the same for the SIS and SIR models 

stats2 <- data.frame(matrix(nrow = 0, ncol = 5))

for(i in 1:length(duration)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 5))
  parms = c(mu = 1/(GenTime(4.6,2.4)), 
            int_timestart = 98, 
            int_timeend = 98 + (duration[i]*7))
  out1 <- data.frame(ode(y = init, func = SIS, times = times, parms = parms))
  temp[,1] <- paste0(as.character(duration[i]), " Weeks")
  temp[,2] <- c(seq(-98,0), seq(1,730-98))
  temp[,3] <- out1$I 
  temp[,4] <- out1$C 
  temp[,5] <- betastatdecrease(times, 
                               as.numeric(parms[2]), 
                               as.numeric(parms[3]))
  stats2 <- rbind(stats2, temp)
}

colnames(stats2) <- c("WeekDuration", "Time", "Infec", "Cum", "Beta")

#Format the data for ggplot
stats2 <- melt(stats2, id.vars = c("WeekDuration", "Time"), measure.vars = c("Infec", "Cum", "Beta"))
stats2$WeekDuration <- factor(stats2$WeekDuration, levels = unique(stats2$WeekDuration))

#Plotting 
p11 <- ggplot(data = stats2[which(stats2$variable == "Infec"),], aes(x = Time, y = value, col = WeekDuration, alpha=WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.4) ,  expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), axis.title.x = element_blank(), axis.text=element_text(size=14), axis.title.y=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0, col = "black", size = 1, lty = 2) +
  scale_alpha_manual(values = c(0.4, 0.4, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1)) +
  scale_color_manual(values=c("#A6CEE3", "#1F78B4","darkred", "#B2DF8A", "#33A02C", "#FB9A99",
                              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A","royalblue"))

p21 <- ggplot(data = stats2[which(stats2$variable == "Beta"),], aes(x = (Time), y = value, col = WeekDuration, alpha=WeekDuration)) + geom_line(size = 1.02, stat = "identity") +
  labs(x ="Time Relative to Lockdown Date (Days)", y = "??") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text=element_text(size=13), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=14), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))  +
  scale_alpha_manual(values = c(0.4, 0.4, 1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1)) +
  scale_color_manual(values=c("#A6CEE3", "#1F78B4","darkred", "#B2DF8A", "#33A02C", "#FB9A99",
                              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A","royalblue"))

plot_grid(p11,  p21, align = "v", nrow = 2, rel_heights = c(0.3, 0.25))

#### SIS Model Tangent Calculation - the 2nd Peak ####

stats3 <- data.frame(matrix(nrow = 0, ncol = 5))

for(i in 1:length(duration)) {
  temp <- data.frame(matrix(nrow = length(seq(0,730)), ncol = 5))
  parms = c(mu = 1/(GenTime(4.6,2.4)), 
            int_timestart = 98, 
            int_timeend = 98 + (duration[i]*7))
  out1 <- data.frame(ode(y = init, func = SIS, times = times, parms = parms))
  temp[,1] <- paste0(as.character(duration[i]), " Weeks")
  temp[,2] <- out1$time #Time
  temp[,3] <- out1$S #Time
  temp[,4] <- out1$I #Time
  temp[,5] <- betastatdecrease(times, 
                               as.numeric(parms[2]), 
                               as.numeric(parms[3]))
  stats3 <- rbind(stats3, temp)
}

colnames(stats3) <- c("WeekDuration", "Time","Susc", "Infec", "Beta")

stats3$WeekDuration <- factor(stats3$WeekDuration, levels = unique(stats3$WeekDuration))

stats3$slope <- stats3$Beta * stats3$S * stats3$I- (1/(GenTime(4.6,2.4))) * stats3$I
stats3$slope1 <- (stats3$Infec[98]-((stats3$Time[98]*stats3$Infec-stats3$Time*stats3$Infec[98])/
                                     (stats3$Time[98]-stats3$Time)))/stats3$Infec[98]

tipping_points <- stats3 %>% filter(Time>110) %>% group_by(WeekDuration) %>% filter(slope1 > slope) %>% top_n(n=1)
tipping_points$dur = c(seq(1:12))

ggplot(tipping_points, aes(x=dur, y=Time-100)) +geom_line(stat="identity", col="darkblue") + 
  scale_x_continuous(name="Duration of intervention (weeks)", breaks = c(1:12)) +ylab(label = "Interval time (days)")
