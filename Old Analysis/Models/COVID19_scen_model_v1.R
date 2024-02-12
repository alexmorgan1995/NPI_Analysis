rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/nCoV Work/Models")
library("deSolve"); library("ggplot2")

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### Baseline Model #### 

SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - beta*S*I
    dI = beta*S*I - mu*I
    dR = mu*I 
    dC = beta*S*I
    return(list(c(dS,dI,dR,dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,300,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))))
outbase <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
outbase[,6] <- as.numeric(parms[2])


#Plotting
outdata <- rbind(data.frame("Compartment" = "Infec", "Times" = outbase[,1], "Prev" = outbase[,3]),
                  data.frame("Compartment" = "Susc", "Times" = outbase[,1], "Prev" = outbase[,2]),
                  data.frame("Compartment" = "Rec", "Times" = outbase[,1], "Prev" = outbase[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = outbase[,1], "Prev" = outbase[,5]),
                  data.frame("Compartment" = "Beta", "Times" = outbase[,1], "Prev" = outbase[,6]))

ggplot(data = outdata[outdata$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.2) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata[outdata$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") +
  scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata[outdata$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

#### Initialisation - Up until 41 Days ####
#Normal SIR Model - with visible and non visible compartments

#C compartment is to count the cumulative infections
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - beta*S*I
    dI = beta*S*I - mu*I
    dR = mu*I 
    dC = beta*S*I
    return(list(c(dS,dI,dR,dC)))
  })
}

# Time = 41, I = 0.0109214845
init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,41,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))))
outbase <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
#Beta Parameter
outbase[,6] <- as.numeric(parms[2])

#Initial Conditions for "continuing" the model - after 41 days
Snew <- as.numeric(outbase[nrow(outbase),2]); Inew <- as.numeric(outbase[nrow(outbase),3]); Rnew <- as.numeric(outbase[nrow(outbase),4])
Cnew <- as.numeric(outbase[nrow(outbase),5])

#### SCENARIO 1 - Static Shift to 0.625*Beta ####

#SIR Model
SIR1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- approxfun(x=c(84,85),y=c(beta_base,beta_int),method="constant", rule  = 2)
    dS = - beta1(time)*S*I
    dI = beta1(time)*S*I - mu*I
    dR = mu*I 
    dC = beta1(time)*S*I
    return(list(c(dS,dI,dR,dC)))
  })
}

init <- c(S = Snew, I = Inew, R = Rnew, C = Cnew)
times <- seq(0,259,by = 1)

#The model (after 41 days) is initiated at 0.625* Beta - After "84 days" the Beta returns to baseline (approxfun)
parms = c(mu = 1/(GenTime(6,2)),
          beta_base = (2*(1/(GenTime(6,2)))*0.625),
          beta_int = (2*(1/(GenTime(6,2)))))
out <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
out <- out[-1,]

#Beta Function - so I can include in the dataframe
beta1 <- approxfun(x=c(84,85),y=c(as.numeric(parms[2]),as.numeric(parms[3])),method="constant", rule  = 2)
out[,6] <- beta1(seq(1,259))
outstat <- rbind(outbase, out)
#Binding pre and post 41 day dataframes together + Renaming the Times - 0 -> 300
outstat[43:301,1] <- seq(42, 300, by = 1)

#Plotting
outdata1 <- rbind(data.frame("Compartment" = "Infec", "Times" = outstat[,1], "Prev" = outstat[,3]),
                 data.frame("Compartment" = "Susc", "Times" = outstat[,1], "Prev" = outstat[,2]),
                 data.frame("Compartment" = "Rec", "Times" = outstat[,1], "Prev" = outstat[,4]),
                 data.frame("Compartment" = "Cumulative", "Times" = outstat[,1], "Prev" = outstat[,5]),
                 data.frame("Compartment" = "Beta", "Times" = outstat[,1], "Prev" = outstat[,6]))

ggplot(data = outdata1[outdata1$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata1[outdata1$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") +
  scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata1[outdata1$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

#### SCENARIO 2 - Linear Increase #### 
#SIR Model
SIR2 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- approxfun(x=c(0,84),y=c(beta_base,beta_int),method="linear", rule  =2)
    dS = - beta1(time)*S*I
    dI = beta1(time)*S*I - mu*I
    dR = mu*I 
    dC = beta1(time)*S*I
    return(list(c(dS,dI,dR,dC)))
  })
}

init <- c(S = Snew, I = Inew, R = Rnew, C = Cnew)
times <- seq(0,259,by = 1)

#The model (after 41 days) becomes 0.25*Beta - then linearly increases back to baseline within 12 weeks
parms = c(mu = 1/(GenTime(6,2)),
          beta_base = (2*(1/(GenTime(6,2)))*0.25),
          beta_int = (2*(1/(GenTime(6,2)))))
#Beta Function - so I can include in the dataframe
beta1 <- approxfun(x=c(0,84),y=c(as.numeric(parms[2]),as.numeric(parms[3])),method="linear", rule  =2)
out <- data.frame(ode(y = init, func = SIR2, times = times, parms = parms))
out <- out[-1,]
out[,6] <- beta1(seq(1,259))
#Binding pre and post 41 day dataframes together + Renaming the Times - 0 -> 300
outstat <- rbind(outbase, out); outstat[43:301,1] <- seq(42, 300, by = 1)

#Plotting
outdata2 <- rbind(data.frame("Compartment" = "Infec", "Times" = outstat[,1], "Prev" = outstat[,3]),
                  data.frame("Compartment" = "Susc", "Times" = outstat[,1], "Prev" = outstat[,2]),
                  data.frame("Compartment" = "Rec", "Times" = outstat[,1], "Prev" = outstat[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = outstat[,1], "Prev" = outstat[,5]),
                  data.frame("Compartment" = "Beta", "Times" = outstat[,1], "Prev" = outstat[,6]))

ggplot(data = outdata2[outdata2$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata2[outdata2$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") +
  scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata2[outdata2$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

#### SCENARIO 3 - Linear Decrease #### 
SIR3 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- approxfun(x=c(0,84),y=c(beta_base,beta_int),method="linear", rule  =2)
    dS = - beta1(time)*S*I
    dI = beta1(time)*S*I - mu*I
    dR = mu*I 
    dC = beta1(time)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = Snew, I = Inew, R = Rnew, C = Cnew)
times <- seq(0,84,by = 1)

#The model (after 41 days) remains Beta - then linearly decreases to 0.25*Beta for 12 weeks - then increases back to baseline after 12 weeks
parms = c(mu = 1/(GenTime(6,2)),
          beta_base = (2*(1/(GenTime(6,2)))),
          beta_int = (2*(1/(GenTime(6,2))))*0.25)
#Beta Function - so I can include in the dataframe
beta1 <- approxfun(x=c(0,84),y=c(as.numeric(parms[2]),as.numeric(parms[3])),method="linear", rule  =2)

out <- data.frame(ode(y = init, func = SIR3, times = times, parms = parms))
out <- out[-1,]
out[,6] <- beta1(seq(1,84))
#Binding pre and post 41 day dataframes together + Renaming the Times - 0 -> 300
outstat <- rbind(outbase, out); outstat[43:126,1] <- seq(42, 125, by = 1)

Snew3 <- as.numeric(outstat[nrow(outstat),2])
Inew3 <- as.numeric(outstat[nrow(outstat),3])
Rnew3 <- as.numeric(outstat[nrow(outstat),4])
Cnew3 <- as.numeric(outstat[nrow(outstat),5])

#Modelling the Return back to Baseline Beta after 12 weeks
init <- c(S = Snew3, I = Inew3, R = Rnew3, C = Cnew3)
times <- seq(0,175,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))))
outadd <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
outadd <- outadd[-1,]
outadd[,6] <- (2*(1/(GenTime(6,2))))
outstatnew <- rbind(outstat, outadd); outstatnew[127:301,1] <- seq(126, 300, by = 1)

#Plotting
outdata3 <- rbind(data.frame("Compartment" = "Infec", "Times" = outstatnew[,1], "Prev" = outstatnew[,3]),
                  data.frame("Compartment" = "Susc", "Times" = outstatnew[,1], "Prev" = outstatnew[,2]),
                  data.frame("Compartment" = "Rec", "Times" = outstatnew[,1], "Prev" = outstatnew[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = outstatnew[,1], "Prev" = outstatnew[,5]),
                  data.frame("Compartment" = "Beta", "Times" = outstatnew[,1], "Prev" = outstatnew[,6]))

ggplot(data = outdata3[outdata3$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata3[outdata3$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") +
  scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata3[outdata3$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

#### SCENARIO 4 - Linear Decrease - Then Increase  #### 

#Normal SIR Model - with visible and non visible compartments
SIR4 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- approxfun(x=c(0,42),y=c(beta_base,beta_int),method="linear", rule  =2)
    dS = - beta1(time)*S*I
    dI = beta1(time)*S*I - mu*I
    dR = mu*I 
    dC = beta1(time)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = Snew, I = Inew, R = Rnew, C = Cnew)
times <- seq(0,42,by = 1)

#The model (after 41 days) remains Beta - then linearly decreases to 0.25*Beta for 6 weeks - then increases back to baseline during weeks 6-12 
parms = c(mu = 1/(GenTime(6,2)),
          beta_base = (2*(1/(GenTime(6,2)))),
          beta_int = (2*(1/(GenTime(6,2))))*0.25)
#Beta Function - so I can include in the dataframe
beta1 <- approxfun(x=c(0,42),y=c(as.numeric(parms[2]),as.numeric(parms[3])),method="linear", rule  =2)

out <- data.frame(ode(y = init, func = SIR4, times = times, parms = parms))
out <- out[-1,]
out[,6] <- beta1(seq(1,42))
#Binding pre and post 41 day (up till week 6) dataframes together + Renaming the Times - 0 -> week 6
outstat <- rbind(outbase, out); outstat[43:84,1] <- seq(42, 83, by = 1)

Snew4 <- as.numeric(outstat[nrow(outstat),2])
Inew4 <- as.numeric(outstat[nrow(outstat),3])
Rnew4 <- as.numeric(outstat[nrow(outstat),4])
Cnew4 <- as.numeric(outstat[nrow(outstat),5])

#Modelling the linear increase back to baseline during week 6 - 12
SIR4_5 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- approxfun(x=c(0,42),y=c(beta_base,beta_int),method="linear", rule  =2)
    dS = - beta1(time)*S*I
    dI = beta1(time)*S*I - mu*I
    dR = mu*I 
    dC = beta1(time)*S*I
    print(beta1(time))
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = Snew4, I = Inew4, R = Rnew4, C = Cnew4)
times <- seq(0,42,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta_base = (2*(1/(GenTime(6,2))))*0.25,
          beta_int = (2*(1/(GenTime(6,2)))))
outadd <- data.frame(ode(y = init, func = SIR4_5, times = times, parms = parms))
outadd <- outadd[-1,]
beta1 <- approxfun(x=c(0,42),y=c(as.numeric(parms[2]),as.numeric(parms[3])),method="linear", rule  =2)
outadd[,6] <- beta1(seq(1,42))
outstatnew <- rbind(outstat, outadd); outstatnew[85:126,1] <- seq(84, 125, by = 1)

#Modelling the baseline Beta after week 12 

init <- c(S = as.numeric(outstatnew[nrow(outstatnew),2]), 
          I = as.numeric(outstatnew[nrow(outstatnew),3]), 
          R = as.numeric(outstatnew[nrow(outstatnew),4]),
          C = as.numeric(outstatnew[nrow(outstatnew),5]))
times <- seq(0,175,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))))
outaddnew <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
outaddnew <- outaddnew[-1,]
outaddnew[,6] <- (2*(1/(GenTime(6,2))))
outstatnewer <- rbind(outstatnew, outaddnew); outstatnewer[127:301,1] <- seq(126, 300, by = 1)

#Plotting
outdata4 <- rbind(data.frame("Compartment" = "Infec", "Times" = outstatnewer[,1], "Prev" = outstatnewer[,3]),
                  data.frame("Compartment" = "Susc", "Times" = outstatnewer[,1], "Prev" = outstatnewer[,2]),
                  data.frame("Compartment" = "Rec", "Times" = outstatnewer[,1], "Prev" = outstatnewer[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = outstatnewer[,1], "Prev" = outstatnewer[,5]),
                  data.frame("Compartment" = "Beta", "Times" = outstatnewer[,1], "Prev" = outstatnewer[,6]))

ggplot(data = outdata4[outdata4$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata4[outdata4$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") +
  scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata4[outdata4$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

#### SCENARIO 5 - Beta Pulse  #### 

#Beta Function to "Pulse" every other week
betapulse <- function(time) {
  ifelse((time >= 7 & time <= 21) | (time >= 35 & time <= 49) | (time >= 63 & time <= 77), 2*(1/(GenTime(6,2)))*0.25, 2*(1/(GenTime(6,2))))
}

SIR4 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betapulse(time)*S*I
    dI = betapulse(time)*S*I - mu*I
    dR = mu*I 
    dC = betapulse(time)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}


init <- c(S = Snew, I = Inew, R = Rnew, C = Cnew)
times <- seq(0,259,by = 1)
parms = c(mu = 1/(GenTime(6,2)))
out <- data.frame(ode(y = init, func = SIR4, times = times, parms = parms))
out <- out[-1,]
#Beta Function - so I can include in the dataframe
out[,6] <- betapulse(seq(1,259))
#Binding pre and post 41 day (up till week 6) dataframes together + Renaming the Times - 0 -> week 6
outstat <- rbind(outbase, out); outstat[43:301,1] <- seq(42, 300, by = 1)

plot(seq(300), betapulse(seq(300)))

#Plotting
outdata5 <- rbind(data.frame("Compartment" = "Infec", "Times" = outstat[,1], "Prev" = outstat[,3]),
                  data.frame("Compartment" = "Susc", "Times" = outstat[,1], "Prev" = outstat[,2]),
                  data.frame("Compartment" = "Rec", "Times" = outstat[,1], "Prev" = outstat[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = outstat[,1], "Prev" = outstat[,5]),
                  data.frame("Compartment" = "Beta", "Times" = outstat[,1], "Prev" = outstat[,6]))

ggplot(data = outdata5[outdata5$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata5[outdata5$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") +
  scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata5[outdata5$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 