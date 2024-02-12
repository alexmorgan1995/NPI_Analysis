library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")
start_time <- Sys.time()

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

combbetamult <- function(scen, time, tstart1, t_dur1, tstart2, t_dur2, R0Dec1, R0Dec2) {
  gamma <- 1/GenTime(3.3, 2.8)
  if(scen == 0) {
    output <-  1.7*gamma
  }
  if(scen == 1) {
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse((time >= (tstart1 + t_dur1 + tstart2)),
                            R0Dec2*gamma,
                            R0Dec1*gamma),
                     1.7*gamma)
  }
  if(scen == 2) {
    R0lin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(R0Dec1, 1.7), method="linear", rule =2)
    R0lin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(R0Dec2, 1.7), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0lin1(time)*gamma,
                            R0lin2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 3) {
    R0lin1 <- approxfun(x=c(tstart1, (tstart1 + t_dur1)), y= c(1.7, R0Dec1), method="linear", rule =2)
    R0lin2 <- approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 + t_dur2)), y= c(1.7, R0Dec2), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0lin1(time)*gamma,
                            R0lin2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 4) {
    R0incdec1 <-approxfun(x=c(tstart1, tstart1+(t_dur1/2), tstart1+t_dur1), y= c(1.7, R0Dec1, 1.7), method="linear", rule =2)
    R0incdec2 <-approxfun(x=c((tstart1 + t_dur1 + tstart2), (tstart1 + t_dur1 + tstart2 +(t_dur2/2)), (tstart1 + t_dur1 + tstart2 +t_dur2)), y= c(1.7, R0Dec2, 1.7), method="linear", rule =2)
    output <- ifelse(((time >= (tstart1) & time <= (tstart1 + t_dur1)) | (time >= (tstart1 + t_dur1 + tstart2) & time <= (tstart1 + t_dur1 + tstart2 + t_dur2))),
                     ifelse(time >= (tstart1) & time <= (tstart1 + t_dur1),
                            R0incdec1(time)*gamma,
                            R0incdec2(time)*gamma),
                     1.7*gamma)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart1 & time <= tstart1+(t_dur1*0.16667)) | (time >= tstart1+(t_dur1*0.3333) & time <= tstart1+(t_dur1*0.5)) |
                       (time >= tstart1+(t_dur1*0.6667) & time <= tstart1+(t_dur1*0.83333)) |
                       (time >= tstart1 + t_dur1 + tstart2 & time <= tstart1 + t_dur1 + tstart2+(t_dur2*0.16667)) | 
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.3333) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.5)) |
                       (time >= tstart1 + t_dur1 + tstart2 + (t_dur2*0.6667) & time <= tstart1 + t_dur1 + tstart2 + (t_dur2*0.83333)),
                     ifelse((time >= tstart1 + t_dur1 + tstart2),
                            R0Dec2*gamma,
                            R0Dec1*gamma),
                     1.7*gamma)
  }
  return(output)
}

plot(seq(0,365),combbetamult(5, seq(0,365), 71, 12*7, 20, 12*7, 0.8, 0.6))

#SEIR set of ODEs
SIRmulti <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbetamult(scen, time, tstart1, t_dur1, tstart2, t_dur2, R0Dec1, R0Dec2)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

# Multiple Intervention Case Study - Trajectory ----------------------------------------

# Run the Model

r0optim <- expand.grid("r01" = seq(0, 1, by = 0.5), "r01" = seq(0, 1, by = 0.5))

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,500,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 1,
          tstart1 = 71,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          R0Dec1 = 0.8,
          R0Dec2 = 0.8)

datalist <- list()

for(j in 1:nrow(r0optim)) {
  
  datalist[[j]] <- local({
    
    j=j
    
    parms["R0Dec1"] = r0optim[j,1]
    parms["R0Dec2"] = r0optim[j,2]
    
    out <- cbind(data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms)), 
                 "r0" = combbetamult(parms[["scen"]], times, parms[["tstart1"]], parms[["t_dur1"]], parms[["tstart2"]], parms[["t_dur2"]],
                                     parms[["R0Dec1"]], parms[["R0Dec2"]])/parms[["gamma"]])
    out$re <- out$r0*out$S

    shade <- data.frame(xmin =  c(parms[["tstart1"]], parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]]), 
                        xmax = c(parms[["tstart1"]]+parms[["t_dur1"]], parms[["tstart1"]]+parms[["t_dur1"]]+parms[["tstart2"]]+parms[["t_dur2"]]), 
                        ymin = 0, ymax = Inf)
    
    p1 <- ggplot(data = out, aes(x = time, y = I)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.075),expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
            plot.title = element_text(size = 18, vjust = 2, hjust = 0.5, face = "bold"),
            axis.title.y=element_text(size=15), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
      scale_x_continuous(expand = c(0, 0))  +
      geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", 
                                                                                     title = paste0(r0optim[j,1], " / ", r0optim[j,2]))
    
    p2 <- ggplot(out, aes(x = time, y = r0))  + theme_bw() +
      scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
            axis.title.y=element_text(size=15),axis.title.x = element_blank(), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_line(size = 1.1, stat = "identity")  + labs(x ="", y = expression(R[0]))
    
    p3 <- ggplot(out, aes(x = time, y = re)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=12),
            axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_hline(yintercept = 1, size = 1.1, lty = 2, col = "black") + geom_line(size = 1.02, stat = "identity") + labs(x ="Time (Days)", y = expression(R[e]))
    
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.45, 0.55))
    return(combplot)
  })
}

combplotmulti_TRAJ <- ggarrange(datalist[[1]], datalist[[2]], datalist[[3]], 
                           datalist[[4]], datalist[[5]], datalist[[6]],
                           datalist[[7]], datalist[[8]], datalist[[9]],
                           nrow = 3, ncol = 3, legend = "none", align = "h")

ggsave(combplotmulti_TRAJ, filename = "1_scenarios_multi_TRAJ.png", dpi = 300, type = "cairo", width = 10, height = 15, units = "in")

# Multiple Intervention Case Study - Trajectory ----------------------------------------

# Run the Model

r0optim <- expand.grid("r01" = seq(0, 1, by = 0.5), "r01" = seq(0, 1, by = 0.5))

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,500,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 2,
          tstart1 = 71,
          t_dur1 = 12*7,
          tstart2 = 42,
          t_dur2 = 12*7,
          R0Dec1 = 0.8,
          R0Dec2 = 0.8)

datalist <- list()

for(j in 1:nrow(r0optim)) {
  
  datalist[[j]] <- local({
    
    j=j
    
    parms["R0Dec1"] = r0optim[j,1]
    parms["R0Dec2"] = r0optim[j,2]
    
    out <- cbind(data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms)), 
                 "r0" = combbetamult(parms[["scen"]], times, parms[["tstart1"]], parms[["t_dur1"]], parms[["tstart2"]], parms[["t_dur2"]],
                                     parms[["R0Dec1"]], parms[["R0Dec2"]])/parms[["gamma"]])
    out$re <- out$r0*out$S
    
    shade <- data.frame(xmin =  c(parms[["tstart1"]], parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]]), 
                        xmax = c(parms[["tstart1"]]+parms[["t_dur1"]], parms[["tstart1"]]+parms[["t_dur1"]]+parms[["tstart2"]]+parms[["t_dur2"]]), 
                        ymin = 0, ymax = Inf)
    
    p1 <- ggplot(data = out, aes(x = time, y = I)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.075),expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
            plot.title = element_text(size = 18, vjust = 2, hjust = 0.5, face = "bold"),
            axis.title.y=element_text(size=15), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
      scale_x_continuous(expand = c(0, 0))  +
      geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", 
                                                                                     title = paste0(r0optim[j,1], " / ", r0optim[j,2]))
    
    p2 <- ggplot(out, aes(x = time, y = r0))  + theme_bw() +
      scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=12),
            axis.title.y=element_text(size=15),axis.title.x = element_blank(), 
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_line(size = 1.1, stat = "identity")  + labs(x ="", y = expression(R[0]))
    
    p3 <- ggplot(out, aes(x = time, y = re)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text=element_text(size=12),
            axis.title.y=element_text(size=18),axis.title.x = element_text(size=18),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
      geom_hline(yintercept = 1, size = 1.1, lty = 2, col = "black") + geom_line(size = 1.02, stat = "identity") + labs(x ="Time (Days)", y = expression(R[e]))
    
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.45, 0.55))
    return(combplot)
  })
}

combplotmulti_TRAJ2 <- ggarrange(datalist[[1]], datalist[[2]], datalist[[3]], 
                                datalist[[4]], datalist[[5]], datalist[[6]],
                                datalist[[7]], datalist[[8]], datalist[[9]],
                                nrow = 3, ncol = 3, legend = "none", align = "h")

ggsave(combplotmulti_TRAJ2, filename = "2_scenarios_multi_TRAJ.png", dpi = 300, type = "cairo", width = 10, height = 15, units = "in")

