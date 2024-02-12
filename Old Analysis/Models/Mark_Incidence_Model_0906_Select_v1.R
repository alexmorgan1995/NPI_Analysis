library("deSolve"); library("ggplot2"); library("reshape2"); library("Cairo"); library("ggpubr"); library("scales")
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
  gamma <- 1/6.5
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
parms = c(gamma = 1/6.5, 
          N = 5.5*10^6,
          zeta = 1/365,
          tstart1 = 21)

times <- seq(0, 365, by = 1)
out1 <- data.frame(ode(y = init, func = SIRS2, times = times, parms = parms))
out1$Beta1 <- beta(times, 2)
out1$N[nrow(out1)]

ggplot(data = out1, aes(x = as.numeric(time), y = I)) + geom_line(size = 1.02, stat = "identity")


#### Run the Model + Obtain Incidences - 1 Week ####

N <- 5.5*10^6
init <- c(S = 1-0.0001, I = 0.0001, R = 0, N = 0)

startdates <- c(21-7, 21, 21+7)

outdata <- data.frame(matrix(nrow = 0, ncol = 6))
cuminf <- numeric(0)

for(i in 1:3) {
  parms = c(gamma = 1/6.5, 
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

(cuminf[3]/cuminf[2])
(1-(cuminf[1]/cuminf[2]))


pinf <- ggplot(data = outdata, aes(x = as.numeric(time), y = as.numeric(I*parms["N"]), col = group))  + theme_bw() +
  labs(x ="Time (Days)", y = "Prevalence", col = "Lockdown Trigger Day") + scale_y_continuous(limits = c(0, 850000),expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=20), legend.text=element_text(size=20),  axis.text=element_text(size=20),
        axis.title.y=element_text(size=20),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  geom_line(size = 1.02, stat = "identity") + scale_colour_manual(
    values = c("red", "black", "blue"),
    breaks = c("14", "21", "28"),
    labels = c("14", "21 (Baseline)", "28"))

#Beta Plot
pbeta <- ggplot(data = outdata, aes(x = (time), y = Beta1, col = group)) + theme_bw() +
  labs(x ="Time (Days)", y = "Beta Value", col = "Lockdown Trigger Day") + scale_y_continuous(limits = c(0,0.5),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=20), legend.text=element_text(size=20),  axis.text=element_text(size=20),
        axis.title.y=element_text(size=20),axis.title.x= element_text(size=20), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  geom_line(size = 1.02, stat = "identity") + scale_colour_manual(
    values = c("red", "black", "blue"),
    breaks = c("14", "21", "28"),
    labels = c("14", "21 (Baseline)", "28"))

#Anti-Aliasing in Plots - Can Ignore 
plot <- ggarrange(pinf, NULL, pbeta, nrow = 3, ncol =1, align = "v", 
                  heights = c(0.9, -0.03, 0.5), 
                  labels = c("A","", "B"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom")

ggsave(plot, filename = "Selec_Com_fig_v1_1week.png", dpi = 300, type = "cairo", width = 8.5, height = 11, units = "in")

#### Run the Model + Obtain Incidences - 1 Day ####

N <- 5.5*10^6
init <- c(S = 1-0.0001, I = 0.0001, R = 0, N = 0)

startdates <- c(21-1, 21, 21+1)

outdata <- data.frame(matrix(nrow = 0, ncol = 6))
cuminf <- numeric(0)

for(i in 1:3) {
  parms = c(gamma = 1/6.5, 
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

(cuminf[3]/cuminf[2])
(1-(cuminf[1]/cuminf[2]))
#### Manipulate Data + Plotting ####

pinf <- ggplot(data = outdata, aes(x = as.numeric(time), y = as.numeric(I*parms["N"]), col = group))  + theme_bw() +
  labs(x ="Time (Days)", y = "Prevalence", col = "Lockdown Trigger Day") + scale_y_continuous(limits = c(0, 850000),expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=20), legend.text=element_text(size=20),  axis.text=element_text(size=20),
        axis.title.y=element_text(size=20),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  geom_line(size = 1.02, stat = "identity") + scale_colour_manual(
    values = c("red", "black", "blue"),
    breaks = c("20", "21", "22"),
    labels = c("20", "21 (Baseline)", "22"))

#Beta Plot
pbeta <- ggplot(data = outdata, aes(x = (time), y = Beta1, col = group)) + theme_bw() +
  labs(x ="Time (Days)", y = "Beta Value", col = "Lockdown Trigger Day") + scale_y_continuous(limits = c(0,0.5),  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=20), legend.text=element_text(size=20),  axis.text=element_text(size=20),
        axis.title.y=element_text(size=20),axis.title.x= element_text(size=20), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) +
  geom_line(size = 1.02, stat = "identity") + scale_colour_manual(
    values = c("red", "black", "blue"),
    breaks = c("20", "21", "22"),
    labels = c("20", "21 (Baseline)", "22"))

#Anti-Aliasing in Plots - Can Ignore 
plot <- ggarrange(pinf, NULL, pbeta, nrow = 3, ncol =1, align = "v", 
                  heights = c(0.9, -0.03, 0.5), 
                  labels = c("A","", "B"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom")

ggsave(plot, filename = "Selec_Com_fig_v1_1Day.png", dpi = 300, type = "cairo", width = 8.5, height = 11, units = "in")


# Plotting R0, Gen Time and Doubling Time  --------------------------------

#DoubTime
DoubTime <- function(R0) {
  T2 =  6.5 / ((R0-1)/log(2))
  return(T2)
}

G <- seq(1.1, 4,by = 0.05)

data <- data.frame("r0" = G, "doubltime" = sapply(G, DoubTime))

R0 <-ggplot(data = data, aes(x = (r0), y = doubltime)) + theme_bw() + geom_line(size = 1.02, stat = "identity") + 
  labs(x = expression(R[0]), y = "Doubling Time (Days)") +
  theme(legend.position = "bottom", legend.title = element_text(size=20), legend.text=element_text(size=20),  axis.text=element_text(size=20),
        axis.title.y=element_text(size=20),axis.title.x= element_text(size=20), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + 
  scale_y_continuous(trans="log10", breaks=c(0,1,10,20,30,40,50), limits = c(1,50),expand = c(0,0),
                     minor_breaks = seq(1, 50, 1)) +
  scale_x_continuous(trans="log10", breaks=c(0,1,2,3,4), limits = c(1,4), expand = c(0,0),
                     minor_breaks = seq(1, 4, 0.1))

ggsave(R0, filename = "Double_Time.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in")
