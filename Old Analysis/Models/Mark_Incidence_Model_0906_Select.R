library("deSolve"); library("ggplot2"); library("reshape2"); library("Cairo"); library("ggpubr")
rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures") # This is where the plots Output

#### Model FUnctions ####

#Generation Time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Beta Parameter
beta <- function(time,tstart1) {
  gamma <- 1/(GenTime(3.3,2.8))
  ifelse((time <= tstart1), 
         2.8*gamma, #Phase 2
         beta1_2 <- 0.8*gamma)
}

plot(beta(seq(0,730), 71))

#Model FUnction
SIRS2 <- function(time, y, parms) {
  
  beta <- beta(time, parms["tstart1"])

  dS = - beta*y["I"]*y["S"] + parms["zeta"]*y["R"] 
  dI = beta*y["I"]*y["S"] - parms["gamma"]*y["I"]
  dR = parms["gamma"]*y["I"] - parms["zeta"]*y["R"]
  dN = beta*y["I"]*y["S"]*parms["N"]
  
  return(list(c(dS, dI, dR, dN)))
}


N <- 5.5*10^6
init <- c(S = 1-0.0001, I = 0.0001, R = 0, N = 0)
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          N = 5.5*10^6,
          zeta = 1/365,
          tstart1 =27)

times <- seq(0, 365, by = 1)
out1 <- data.frame(rk(y = init, func = SIRS2, times = times, parms = parms, method = "rk4"))
out1$Beta1 <- beta(times, 27)
out1$N[nrow(out1)]

ggplot(data = out1, aes(x = as.numeric(time), y = I)) + geom_line(size = 1.02, stat = "identity")


#### Run the Model + Obtain Incidences ####

N <- 5.5*10^6
init <- c(S = 1-0.0001, I = 0.0001, R = 0, N = 0)

startdates <- c(27-3.3, 27, 27+3.3)

outdata <- data.frame(matrix(nrow = 0, ncol = 6))
cuminf <- numeric(0)

for(i in 1:3) {
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            N = 5.5*10^6,
            zeta = 1/365,
            tstart1 = startdates[i])
  
  times <- seq(0, 365, by = 1)
  out1 <- data.frame(ode(y = init, func = SIRS2, times = times, parms = parms))
  out1$Beta1 <- beta(times, startdates[i])
  cuminf <- c(cuminf, out1$N[nrow(out1)])
  out1$group <- as.character(startdates[i])
  outdata <- rbind.data.frame(outdata, out1)
}
  
out1$N[nrow(out1)]

(1-(cuminf[3]/cuminf[2]))*100

#### Manipulate Data + Plotting ####

pinf <- ggplot(data = outdata, aes(x = as.numeric(time), y = as.numeric(I*parms["N"]), col = group))  + theme_bw() +
  labs(x ="Time (Days)", y = "Prevalence", col = "Lockdown Trigger Day") + scale_y_continuous(limits = c(0, 300000),expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=20), legend.text=element_text(size=20),  axis.text=element_text(size=20),
        axis.title.y=element_text(size=20),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  geom_line(size = 1.02, stat = "identity")

#Beta Plot
pbeta <- ggplot(data = outdata, aes(x = (time), y = Beta1, col = group)) + theme_bw() +
  labs(x ="Time (Days)", y = "Beta Value", col = "Lockdown Trigger Day") + scale_y_continuous(limits = c(0,0.35),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=20), legend.text=element_text(size=20),  axis.text=element_text(size=20),
        axis.title.y=element_text(size=20),axis.title.x= element_text(size=20), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  geom_line(size = 1.02, stat = "identity")

#Anti-Aliasing in Plots - Can Ignore 
plot <- ggarrange(pinf, NULL, pbeta, nrow = 3, ncol =1, align = "v", 
                  heights = c(0.9, -0.03, 0.5), 
                  labels = c("A","", "B"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom")

ggsave(plot, filename = "Selec_Com_fig.png", dpi = 300, type = "cairo", width = 8.5, height = 11, units = "in")
