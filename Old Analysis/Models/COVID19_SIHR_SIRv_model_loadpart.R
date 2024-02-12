rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/nCoV Work/Models")
library("deSolve"); library("ggplot2"); library("bayestestR"); library("svMisc")

#### Model Functions ####

#Normal SIR Model - with visible and non visible compartments
SIVR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - v*(beta*S*(I+Inv)) - (1-v)*(beta*S*(I+Inv))
    dI = v*(beta*S*(I+Inv)) - mu*I
    dInv = (1-v)*(beta*S*(I+Inv)) - mu*Inv
    dR = mu*I + mu*Inv
    dCI = v*(beta*S*(I+Inv))
    dCInv = (1-v)*(beta*S*(I+Inv))
    return(list(c(dS,dI,dInv,dR,dCI,dCInv)))
  })
}

#SIHR Model 
SIHR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = -beta1*S*I -beta2*S*I
    dI = beta1*S*I - mu*I
    dH = beta2*S*I - mu*H
    dR = mu*I + mu*H
    dC1 = beta1*S*I
    dC2 = beta2*S*I
    return(list(c(dS,dI,dH,dR,dC1,dC2)))
  })
}

#Function for Generation time
#Log function is actually the natural log function 
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### Obtaining SIHR Output ####

init <- c(S = 0.9999, I = 0.0001, H = 0, R = 0, C1= 0, C2 = 0)

times <- seq(0,200,by = 1)

parms = c(mu = 1/(GenTime(6,2)),
          beta1 = 2*(1/(GenTime(6,2))),
          beta2 = (2*(1/(GenTime(6,2))))*5
          )

out <- ode(y = init, func = SIHR, times = times, parms = parms)

outdata <- rbind(data.frame("Compartment" = "Infec", "Times" = times, "Prev" = out[,3]),
                   data.frame("Compartment" = "Hidden_Inf", "Times" = times, "Prev" = out[,4]),
                   data.frame("Compartment" = "Cum_Inf", "Times" = times, "Prev" = out[,6]),
                   data.frame("Compartment" = "Cum_Hid_Inf", "Times" = times, "Prev" = out[,7]))

#### SIHR Model Plotting ####

ggplot(data = outdata[outdata$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.03), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata[outdata$Compartment == "Hidden_Inf",], aes(x = Times, y = Prev)) + geom_line(size = 1.05) +
  labs(x ="Time (Days)", y = "Proportion of Hidden Infected Humans") +
  scale_y_continuous(limits = c(0,0.25), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata[outdata$Compartment == "Cum_Inf",], aes(x = Times, y = Prev)) + geom_line(size = 1.05) +
  labs(x ="Time (Days)", y = "Cumulative Number of Human Infections") +
  scale_y_continuous(limits = c(0,0.25), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata[outdata$Compartment == "Cum_Hid_Inf",], aes(x = Times, y = Prev)) + geom_line(size = 1.05) +
  labs(x ="Time (Days)", y = "Cumulative Number of Hidden Human Infections") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

#### ABC Model Fitting - Fitting a SIR model to the SIHR model ####

#SIHR Initialisation
init <- c(S = 0.9999, I = 0.0001, H = 0, R = 0, C1= 0, C2 = 0)

times <- seq(0,200,by = 1)

parms = c(mu = 1/(GenTime(6,2)),
          beta1 = 2*(1/(GenTime(6,2))),
          beta2 = (2*(1/(GenTime(6,2))))*5)

outSIHR <- ode(y = init, func = SIHR, times = times, parms = parms)

#SIVR Initialisation
init1 <- c(S = 0.9999, I = 0.0001, Inv = 0, R = 0, CI = 0, CInv = 0)

times <- seq(0,200,by = 1)

parms1 = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))),
          v = 0.05)

outSIVR <- ode(y = init1, func = SIVR, times = times, parms = parms1)

#SUMMARY STATISTICS 
#Timing of Outbreak Peak

sum_time_max <- function(prev) {
  return(prev[,1][which(prev[,3] == max(prev[,3]))])
}

#Cases at Outbreak Peak
sum_prev_max <- function(prev) {
  return(prev[,3][which(prev[,3] == max(prev[,3]))])
}

#Dataobs must be a vector of the parameters we observe 

dist_function <- function(sum.stats, data.obs, model.obs) {
  peak_time <- sapply(sum.stats[1], function(x) {abs(data.obs[1] - x(model.obs))})
  peak_max <- sapply(sum.stats[2], function(x) {abs(data.obs[2] - x(model.obs))})
  return(c(peak_time, peak_max))
}

#This I guess is very similar to the for loop I programmed before - which will output a distance for any given parameter set or initial conditions. 

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, thetaparm, init.state, times, data) {
  temp <- matrix(NA, nrow = length(times), ncol = 3)
  parms2 = c(beta = thetaparm[["beta"]], 
             mu = thetaparm[["mu"]], 
             v = thetaparm[["v"]])
  out <- ode(y = init.state, func = fitmodel, times = times, parms = parms2)
  temp[,1] <- out[,1]
  temp[,2] <- out[,3]
  temp[,3] <- out[,4]
  colnames(temp) <- c("Times", "VisInf", "NonVisInf")
  return(distanceABC(sum.stats, data, temp))
}

computeDistanceABC_ALEX(sum.stats = list(sum_time_max, sum_prev_max), 
                        distanceABC = dist_function, 
                        fitmodel = SIVR, 
                        thetaparm = c(mu = 1/(GenTime(6,2)),
                                      beta = (2*(1/(GenTime(6,2)))),
                                      v = 0.05), 
                        init.state = c(S = 0.9999, I = 0.0001, Inv = 0, R = 0, CI = 0, CInv = 0), 
                        times = times <- seq(0,100,by = 1), 
                        data = c(63, 0.02566227))

# use the ABC rejection algorithm to find population 1 in the ABC-SMC algorithm
start_time <- Sys.time()

ABC_algorithm <- function(N, epsilon, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data) {
  dump1 <- matrix(nrow = 0, ncol = 5)
  i <- 0
  while(i < N) {
    thetaparm <- c(beta = runif(1, min = 0, max = 1),
                   mu = runif(1, min = 0, max = 1),
                   v = runif(1, min = 0, max = 1))
    dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, thetaparm, init.state, times, data)
    if((dist[1] <= epsilon[1]) && (dist[2] <= epsilon[2]) && (!is.na(dist))) {
      thetaparm <- c(thetaparm, dist)
      dump1 <- rbind(dump1, thetaparm)
    }
    i <- dim(dump1)[1]
    print(i/N)
  }
  return(dump1)
}

data_ABC <- ABC_algorithm(N = 1000, 
                          epsilon = c((63*0.05), (0.02566227*0.05)),
                          sum.stats = list(sum_time_max, sum_prev_max), 
                          distanceABC = dist_function, 
                          fitmodel = SIVR, 
                          init.state = c(S = 0.9999, I = 0.0001, Inv = 0, R = 0, CI = 0, CInv = 0), 
                          times = seq(0,100,by = 1), 
                          data = c(63, 0.02566227))

end_time <- Sys.time(); end_time - start_time

#### Particle Import - N = 6000 ####

combpart <- rbind(data.frame(read.csv("COVID_fit_v1.csv")),
                  data.frame(read.csv("COVID_fit_v2.csv")),
                  data.frame(read.csv("COVID_fit_v3.csv")),
                  data.frame(read.csv("COVID_fit_v4.csv")),
                  data.frame(read.csv("COVID_fit_v5.csv")))
                
#### Plotting ####
#Finding the Posterior Distribution Mode - Maximum a posteriori  
plot(density(runif(2500, min = 0, max = 1)), xlim = c(0, 1), lwd = 1.2, xlab = "Beta") #BetaAA
plot(density(runif(2500, min = 0, max = 1)), xlim = c(0, 1), lwd = 1.2, xlab = "Mu") #Phi
plot(density(runif(2500, min = 0, max = 1)), xlim = c(0, 1), lwd = 1.2, xlab = "V") #Theta

#Posterior Distributions
#Beta_AA
betaAA <- data.frame("value" = combpart[,2])
betaAA_MAP <- data.frame("MAP" = map_estimate(combpart[,2]))

ggplot(betaAA, aes(x=value)) + geom_density(alpha=.3, fill = "red")  + 
  geom_vline(data = betaAA_MAP, aes(xintercept = MAP), size=1.05, linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + scale_y_continuous(expand = c(0, 0), limits = c(0,7))+
  geom_text(data = betaAA_MAP, aes(x = MAP + 0.1, label= signif(MAP, digits = 3), y= 5.5), show.legend=FALSE, size = 5) +
  labs(x ="Beta_AA Parameter Value", y = "Density") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

#1/Generation Time - Mu 
mu <- data.frame("value" = combpart[,3])
mu_MAP <- data.frame("MAP" = map_estimate(combpart[,3]))

ggplot(mu, aes(x=value)) + geom_density(alpha=.3, fill = "red")  + 
  geom_vline(data = mu_MAP, aes(xintercept = MAP), size=1.05, linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + scale_y_continuous(expand = c(0, 0), limits = c(0,6))+
  geom_text(data = mu_MAP, aes(x = MAP + 0.1, label= signif(MAP, digits = 3), y= 5.5), show.legend=FALSE, size = 5) +
  labs(x ="Mu (1/Generation Time)", y = "Density") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

#Proportion of non-Visible Infections - v
vparm <- data.frame("value" = combpart[,4])
vparm_MAP <- data.frame("MAP" = map_estimate(combpart[,4]))

ggplot(vparm, aes(x=value)) + geom_density(alpha=.3, fill = "red")  + 
  geom_vline(data = vparm_MAP, aes(xintercept = MAP), size=1.05, linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1.75))+
  geom_text(data = vparm_MAP, aes(x = MAP + 0.07, label= signif(MAP, digits = 3), y= 1.5), show.legend=FALSE, size = 5) +
  labs(x ="V (Proportion of non-Visible Infections)", y = "Density") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))





#### Model Fit Test ####

#SIHR Initialisation 

init <- c(S = 0.9999, I = 0.0001, H = 0, R = 0, C1= 0, C2 = 0)
  
times <- seq(0,200,by = 1)
  
parms = c(mu = 1/(GenTime(6,2)),
            beta1 = 2*(1/(GenTime(6,2))),
            beta2 = (2*(1/(GenTime(6,2))))*5)
  
outSIHR <- ode(y = init, func = SIHR, times = times, parms = parms)

outSIHRdata <- rbind(data.frame("Compartment" = "Infec", "Times" = times, "Prev" = outSIHR[,3]),
                 data.frame("Compartment" = "Hidden_Inf", "Times" = times, "Prev" = outSIHR[,4]),
                 data.frame("Compartment" = "Cum_Inf", "Times" = times, "Prev" = outSIHR[,6]),
                 data.frame("Compartment" = "Cum_Hid_Inf", "Times" = times, "Prev" = outSIHR[,7]))

#SIVR Initialisation
init1 <- c(S = 0.9999, I = 0.0001, Inv = 0, R = 0, CI = 0, CInv = 0)
  
times <- seq(0,200,by = 1)
  
parms1 = c(mu = 0.371, beta = 0.482, v = 0.115)

outSIVR <- ode(y = init1, func = SIVR, times = times, parms = parms1)

outSIVRdata <- rbind(data.frame("Compartment" = "Infec", "Times" = times, "Prev" = outSIVR[,3]),
                 data.frame("Compartment" = "nonVis_Inf", "Times" = times, "Prev" = outSIVR[,4]),
                 data.frame("Compartment" = "Cum_Inf", "Times" = times, "Prev" = outSIVR[,6]),
                 data.frame("Compartment" = "Cum_nonVis_Inf", "Times" = times, "Prev" = outSIVR[,7]))

#Comparison Plots 
SIVR_SIHR_inf <- rbind(data.frame("State" = "SIHR_inf", "Value" = outSIHR[,3], "Time" = outSIHR[,1]),
                       data.frame("State" = "SIVR_inf", "Value" = outSIVR[,4], "Time" = outSIVR[,1]))

ggplot(data = SIVR_SIHR_inf, aes(x = Time, y = Value, color = State)) + geom_line(size = 1.05) +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.03)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + 
  scale_color_manual(values=c("darkblue", "darkred"), labels=c("SIHR non-Hidden Infecteds", "SIVR Visible Infecteds")) +
  scale_x_continuous(expand = c(0, 0)) 
