library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/18_09_20")

# Generic Model Functions ----------------------------------------------------------
#Function for the generation time - a function of R0 and the doubling time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - (beta*S*I)/N
    dI = (beta*S*I)/N - gamma*I + c
    dR = gamma*I 
    
    dC = (beta*S*I)/N + c
    return(list(c(dS, dI, dR, dC)))
  })
} 

# Epidemic Trajectory Plots for the 5 NPI scenarios -----------------------------------

#Initial Conditions and Parameters
init <- c(S = 5000000-2500, I = 2500, R = 0, C = 0)
times <- seq(0,720,by = 1)
gamma <- (1/GenTime(3, 2.8))
parms = c(gamma = 1/GenTime(3, 2.8),
          beta = 1.8*gamma, #R should be from 1.1 to 2 in steps of 0.1
          c = 5,
          N = 5000000)

#For Loop to run the 5 NPI scenarios and create plots with I(t), Beta(t) and Re(t)

data <- data.frame(matrix(nrow = 9, ncol = 0))

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

for(i in 1:2) { #Run the model for unmitigated and NPI scenarios
  
  parms["c"] <- c(0,5)[i]
  
  #Run the model and store/calculate important run characteristics
  out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)), 
               "group" =  c("No Importation", paste0(as.character(parms["c"]), " cases per day"))[i])
  data <- rbind(data, out)
}

#Convert the dataframe into a suitable format for the model plotting
plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
plotC <- melt(data, id.vars = c("time", "group"), measure.vars = ("C"))

p1 <- ggplot(data = plotdata, aes(x = time, y = value, col= group)) + theme_bw() +
  scale_y_continuous(limits = c(0, max(data$I)+10000), expand = c(0,0)) + scale_x_continuous(limits = c(0, out$time[out$I < 2500][1]), expand = c(0, 0)) + 
  scale_color_manual(values = c("darkred", "darkblue")) + 
  labs(x = "Time", y = "Prevalence", col = "Importation", title = bquote(italic(R)*" = 1.8")) +
  theme(axis.title.y=element_text(size=16), legend.text = element_text(size=14), legend.title =  element_text(size=14),  axis.title.x=element_text(size=16), plot.title = element_text(size = 22, vjust = 3, hjust = 0.5, face = "bold"), 
        axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
  geom_line(size = 1.1, stat = "identity")

p2 <- ggplot(data = plotC, aes(x = time, y = value, col= group)) + theme_bw() +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(limits = c(0, out$time[out$I < 2500][1]), expand = c(0, 0)) + 
  scale_color_manual(values = c("darkred", "darkblue")) + 
  labs(x = "Time", y = "Cumulative Incidence", col = "Importation", title = NULL) +
  theme(axis.title.y=element_text(size=16), legend.text = element_text(size=14), legend.title =  element_text(size=14),  axis.title.x=element_text(size=16), plot.title = element_text(size = 22, vjust = 3, hjust = 0.5, face = "bold"), 
        axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
  geom_line(size = 1.1, stat = "identity")

#Final Size at 2500
cat2500No <- data[data$group == "No Importation",][data$I < 2500,][1,5]
cat2500Import <- data[data$group == paste0(as.character(parms["c"]), " cases per day"),][data$I < 2500,][1,5]

#Cases Due to Importation
1-(cat2500No/cat2500Import)

combplot <- ggarrange(p1,p2,nrow = 2, ncol = 1 , common.legend = TRUE, legend = "bottom")

ggsave(combplot, filename = paste0(parms["beta"]/gamma,"_",parms["c"],"_import.png"), dpi = 300, type = "cairo", width = 7, height = 10, units = "in")

