library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures")

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#FUnction to Obtain Betas for the Model 
combbeta <- function(time, tstart, t_dur, R0Dec) {
  gamma <- 1/GenTime(6, 2)
  R0lin <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(2, R0Dec), method="linear", rule =2) 
  output <- ifelse((time >= tstart & time <= tstart + t_dur),
                   R0lin(time)*gamma, 
                   2*gamma)
  return(output)
}

plot(seq(0,125),combbeta(seq(0,125), 41, 12*7, 0.5))

#SEIR set of ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(time, tstart, t_dur, R0Dec)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 


# Baseline Model Analysis - 5 Scenarios -----------------------------------

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,125,by = 1)
parms = c(gamma = 1/GenTime(6, 2),
          tstart = 41,
          t_dur = 12*7,
          R0Dec = 0.5)

out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)), 
             "r0" = combbeta(times, parms[["tstart"]], parms[["t_dur"]], parms[["R0Dec"]])/parms[["gamma"]])

out$re <- out$r0*out$S

test <- ggplot(data = out, aes(x = time, y = I)) + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
        plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
        axis.title.y=element_text(size=18), axis.title.x = element_text(size=18), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence") +
  scale_x_continuous(expand = c(0,0)) + 
  geom_hline(yintercept = 0.0070940951	, size = 1.02, lty = 2, col = "red") +
  geom_vline(xintercept = 76, size = 1.02, lty = 2, col = "black", alpha = 0.5) + 
  annotate(geom = "rect", xmin = 41, ymin = 0, xmax = 41+12*7, ymax = Inf, fill = "darkblue", alpha = 0.2)  + 
  scale_y_continuous(trans="log10", breaks=c(0.001, 0.01, 0.1), limits = c(0.001,0.15),expand = c(0,0)) 


ggsave(test, filename = "Scen3_Request.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in")

