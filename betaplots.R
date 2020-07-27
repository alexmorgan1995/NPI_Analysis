library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#FUnction to Obtain Betas for the Model 
combbeta <- function(scen, time, tstart, t_dur, R0Dec) {
  gamma <- 1/GenTime(3, 2.8)
  if(scen == 0) {
    output <-  1.96*gamma
  }
  if(scen == 1) {
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     R0Dec*gamma, 1.96*gamma)
  }
  if(scen == 2) {
    R0lin <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(R0Dec, 1.96), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     R0lin(time)*gamma, 1.96*gamma)
  }
  if(scen == 3) {
    R0lin <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(1.96, R0Dec), method="linear", rule =2) 
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     R0lin(time)*gamma, 1.96*gamma)
  }
  if(scen == 4) {
    R0incdec <-approxfun(x=c(tstart, tstart+(t_dur/2), tstart+t_dur), y= c(1.96, R0Dec, 1.96), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     R0incdec(time)*gamma, 1.96*gamma)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart & time <= tstart+(t_dur*0.16667)) | (time >= tstart+(t_dur*0.3333) & time <= tstart+(t_dur*0.5)) |
                       (time >= tstart+(t_dur*0.6667) & time <= tstart+(t_dur*0.83333)), 
                     R0Dec*gamma, 1.96*gamma)
  }
  return(output)
}

plot(seq(0,365),combbeta(4, seq(0,365), 41, 12*7, 2.1*0.4))

#SEIR set of ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, R0Dec)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 


# Model Analysis - 5 Scenarios -----------------------------------
init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 3),
          scen = 3,
          tstart = 52,
          t_dur = 12*7,
          R0Dec = 0.784)

betaplots <- list()

for(i in 1:5) {
  betaplots[[i]] <- local({
    parms["scen"] = seq(1:5)[i]
    out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)),
                 "beta" = combbeta(parms["scen"], times, parms[["tstart"]], parms[["t_dur"]], parms[["R0Dec"]]))
    out$re <- out$beta/parms[["gamma"]]*out$S    
    
    p1 <- ggplot(data = out, aes(x = time, y = beta)) + theme_bw() + geom_line(size = 1.1, stat = "identity") + 
      scale_y_continuous(limits = c(0 , 0.3), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=18),axis.title.x = element_text(size=15), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = expression(beta[(t)]))
    
    return(p1)
    
  })
  ggsave(betaplots[[i]], filename = paste0("Beta_Scenario_", seq(1:5)[i],".png"), dpi = 300, type = "cairo", width = 6, height = 3, units = "in")
}


# Cmin analysis -----------------------------------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 3),
          scen = 3,
          tstart = 52,
          t_dur = 12*7,
          R0Dec = 0.784)

betaplots <- list()

for(i in 1:5) {
  betaplots[[i]] <- local({
    parms["scen"] = seq(1:5)[i]
    out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)),
                 "beta" = combbeta(parms["scen"], times, parms[["tstart"]], parms[["t_dur"]], parms[["R0Dec"]]))
    out$re <- out$beta/parms[["gamma"]]*out$S    
    
    p1 <- ggplot(data = out, aes(x = time, y = beta)) + theme_bw() + geom_line(size = 1.1, stat = "identity") + 
      scale_y_continuous(limits = c(0 , 0.3), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=18),axis.title.x = element_text(size=15), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = expression(beta[(t)]))
    
    return(p1)
    
  })
  ggsave(betaplots[[i]], filename = paste0("Beta_Scenario_", seq(1:5)[i],".png"), dpi = 300, type = "cairo", width = 6, height = 3, units = "in")
}

