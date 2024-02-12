rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/nCoV Work"); library("deSolve"); library("ggplot2")

#### NORMAL SIR ####
#### Model Functions ####
#Normal SIR Model

SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = -beta*S*I
    dI = beta*S*I - mu*I
    dR = mu*I
    dC = beta*S*I
    return(list(c(dS,dI,dR,dC)))
  })
}

#### Integration Steps + Dataframe ####
#Initial COnditions

init <- c(S = 0.9999, I = 0.0001, R = 0, C= 0 )

times <- seq(0,100,by = 1)

#Need to Specify Model Parameters

parms = c(beta = 4*(1/2),
          mu = 1/8)


out <- ode(y = init, func = SIR, times = times, parms = parms)

plot(x = out[,1], y = out[,3], type= "l")

#Multiple 

parmframe <- data.frame("mu" = c(1/4, 1/5, 1/7, 1/10), "beta" = c((2*(1/4)), (2.5*(1/5)), (3*(1/7)), (4*(1/10))), 
                        "Scenario" = c("One","Two","Three","Four"))

out <- ode(y = init, func = SIR, times = times, parms = parms)
output <- data.frame()

for (i in 1:nrow(parmframe)) {
  temp <- data.frame(matrix(ncol = 4, nrow = length(times)))
  parms1 = c(beta = parmframe$beta[i], mu = parmframe$mu[i])
  out <- ode(y = init, func = SIR, times = times, parms = parms1)
  temp[,1] <- times
  temp[,2] <- out[,3]
  temp[,3] <- parmframe$Scenario[i]
  temp[,4] <- out[,5]
  output <- rbind(output, temp)
}

colnames(output) <- c("Time","InfFrac","Scen", "CumInf") 

#### Plotting ####
ggplot(data = output, aes(x = Time, y = InfFrac, color = Scen)) + 
  geom_line(size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.5), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_colour_discrete(name  ="Scenarios",
                        breaks=c("One","Two","Three","Four"),
                        labels=c(expression(paste("R"[0],"= 2,  ", frac(1, mu), "  = 4")),
                                 expression(paste("R"[0],"= 2.5, ", frac(1, mu), "  = 5")),
                                 expression(paste("R"[0],"= 3, ", frac(1, mu), "  = 7")),
                                 expression(paste("R"[0],"= 4, ", frac(1, mu), "  = 10")))) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = output, aes(x = Time, y = CumInf, color = Scen)) + 
  geom_line(size = 1.1) +
  labs(x ="Time (Days)", y = "Cumulative Number of Infections") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_colour_discrete(name  ="Scenarios",
                        breaks=c("One","Two","Three","Four"),
                        labels=c(expression(paste("R"[0],"= 2,  ", frac(1, mu), "  = 4")),
                                 expression(paste("R"[0],"= 2.5, ", frac(1, mu), "  = 5")),
                                 expression(paste("R"[0],"= 3, ", frac(1, mu), "  = 7")),
                                 expression(paste("R"[0],"= 4, ", frac(1, mu), "  = 10")))) 

#### ADAPTED SIR - HIDDEN COMPARTMENT ####
#w/ "Hidden" Compartment - SIR Model

SIRHid <- function(time, state, parameters) {
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

#### Integration + Dataframe ####

init <- c(S = 0.9999, I= 0.0001, H= 0, R = 0, C1 = 0, C2 = 0)
times <- seq(0,100,by = 1)

#Need to Specify Model Parameters

output <- data.frame()

betascen <- c(5,10,20)
parmframe <- data.frame("mu" = c(1/4, 1/5, 1/7, 1/10), "beta" = c((2*(1/4)), (2.5*(1/5)), (3*(1/7)), (4*(1/10))), 
                        "Scenario" = c("One","Two","Three","Four"))

for (j in 1:length(betascen)) {
  outputdump <- data.frame() 
  for (i in 1:nrow(parmframe)) {
    temp <- data.frame(matrix(ncol = 7, nrow = length(times)))
    parms1 = c(beta1 = parmframe$beta[i], 
               beta2 = parmframe$beta[i]*betascen[j], 
               mu = parmframe$mu[i])
    out <- ode(y = init, func = SIRHid, times = times, parms = parms1)
    temp[,1] <- times
    temp[,2] <- out[,3]
    temp[,3] <- out[,4]
    temp[,4] <- parmframe$Scenario[i]
    temp[,5] <- out[,6]
    temp[,6] <- out[,7]
    temp[,7] <- betascen[j]
    outputdump <- rbind(outputdump, temp)
  }
  output <- rbind(output, outputdump)
}

colnames(output) <- c("Time","InfFrac","HidInf", "Scen", "CumInfI", "CumInfH", "BetaScen") 

#### Normal Infections Plotting ####

ggplot(data = output[output$BetaScen == "20",], aes(x = Time, y = InfFrac, color = Scen)) + 
  geom_line(size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.1), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_colour_discrete(name  ="Scenarios",
                        breaks=c("One","Two","Three","Four"),
                        labels=c(expression(paste("R"[0],"= 2,  ", frac(1, mu), "  = 4")),
                                 expression(paste("R"[0],"= 2.5, ", frac(1, mu), "  = 5")),
                                 expression(paste("R"[0],"= 3, ", frac(1, mu), "  = 7")),
                                 expression(paste("R"[0],"= 4, ", frac(1, mu), "  = 10")))) +
  scale_x_continuous(expand = c(0, 0)) 


ggplot(data = output[output$BetaScen == "20",], aes(x = Time, y = CumInfI, color = Scen)) + 
  geom_line(size = 1.1) +
  labs(x ="Time (Days)", y = "Cumulative Number of Infections") +
  scale_y_continuous(limits = c(0,0.25), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_colour_discrete(name  ="Scenarios",
                        breaks=c("One","Two","Three","Four"),
                        labels=c(expression(paste("R"[0],"= 2,  ", frac(1, mu), "  = 4")),
                                 expression(paste("R"[0],"= 2.5, ", frac(1, mu), "  = 5")),
                                 expression(paste("R"[0],"= 3, ", frac(1, mu), "  = 7")),
                                 expression(paste("R"[0],"= 4, ", frac(1, mu), "  = 10"))))


#### Hidden Compartment ####

ggplot(data = output[output$BetaScen == "20",], aes(x = Time, y = HidInf, color = Scen)) + 
  geom_line(size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.5), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_colour_discrete(name  ="Scenarios",
                        breaks=c("One","Two","Three","Four"),
                        labels=c(expression(paste("R"[0],"= 2,  ", frac(1, mu), "  = 4")),
                                 expression(paste("R"[0],"= 2.5, ", frac(1, mu), "  = 5")),
                                 expression(paste("R"[0],"= 3, ", frac(1, mu), "  = 7")),
                                 expression(paste("R"[0],"= 4, ", frac(1, mu), "  = 10")))) +
  scale_x_continuous(expand = c(0, 0)) 


ggplot(data = output[output$BetaScen == "20",], aes(x = Time, y = CumInfH, color = Scen)) + 
  geom_line(size = 1.1) +
  labs(x ="Time (Days)", y = "Cumulative Number of Infections") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_colour_discrete(name  ="Scenarios",
                        breaks=c("One","Two","Three","Four"),
                        labels=c(expression(paste("R"[0],"= 2,  ", frac(1, mu), "  = 4")),
                                 expression(paste("R"[0],"= 2.5, ", frac(1, mu), "  = 5")),
                                 expression(paste("R"[0],"= 3, ", frac(1, mu), "  = 7")),
                                 expression(paste("R"[0],"= 4, ", frac(1, mu), "  = 10"))))