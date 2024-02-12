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
times <- seq(0,35,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))))
outbase <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
#Beta Parameter
outbase[,6] <- as.numeric(parms[2])

#Initial Conditions for "continuing" the model - after 41 days
Snew <- as.numeric(outbase[nrow(outbase),2]); Inew <- as.numeric(outbase[nrow(outbase),3]); Rnew <- as.numeric(outbase[nrow(outbase),4])
Cnew <- as.numeric(outbase[nrow(outbase),5])

#### SCENARIO 3 - Linear Decrease #### 
SIR3 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- approxfun(x=c(0,112),y=c(beta_base,beta_int),method="linear", rule  =2)
    dS = - beta1(time)*S*I
    dI = beta1(time)*S*I - mu*I
    dR = mu*I 
    dC = beta1(time)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}

init <- c(S = Snew, I = Inew, R = Rnew, C = Cnew)
times <- seq(0,112,by = 1)

#The model (after 41 days) remains Beta - then linearly decreases to 0.25*Beta for 12 weeks - then increases back to baseline after 12 weeks
parms = c(mu = 1/(GenTime(6,2)),
          beta_base = (2*(1/(GenTime(6,2)))),
          beta_int = (2*(1/(GenTime(6,2))))*0.25)
#Beta Function - so I can include in the dataframe
beta1 <- approxfun(x=c(0,112),y=c(as.numeric(parms[2]),as.numeric(parms[3])),method="linear", rule  =2)
out <- data.frame(ode(y = init, func = SIR3, times = times, parms = parms))
out <- out[-1,]
out[,6] <- beta1(seq(1,112))
#Binding pre and post 41 day dataframes together + Renaming the Times - 0 -> 300
outstat <- rbind(outbase, out)

outstat[37:148,1] <- seq(36, 147, by = 1)

Snew3 <- as.numeric(outstat[nrow(outstat),2])
Inew3 <- as.numeric(outstat[nrow(outstat),3])
Rnew3 <- as.numeric(outstat[nrow(outstat),4])
Cnew3 <- as.numeric(outstat[nrow(outstat),5])

#Modelling the Return back to Baseline Beta after 12 weeks
init <- c(S = Snew3, I = Inew3, R = Rnew3, C = Cnew3)
times <- seq(0,253,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))))
outadd <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
outadd <- outadd[-1,]
outadd[,6] <- (2*(1/(GenTime(6,2))))
outstatnew <- rbind(outstat, outadd); outstatnew[149:401,1] <- seq(148, 400, by = 1)

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

ggplot(data = outdata3[outdata3$Compartment == "Susc",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkgreen") +
  labs(x ="Time (Days)", y = "Proportion of Susceptibles") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 