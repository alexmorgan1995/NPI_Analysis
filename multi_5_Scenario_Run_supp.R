library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis/supplementary")
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


# Multiple Optimisations - Trigger + Length --------------------------------------------------

lengthdata <- expand.grid("int1length" = seq(3,9, by = 3), "int2length" = seq(3,9, by = 3))
optimdata <- expand.grid("tstart1" = seq(0,125, by = 5), "tstart2" = seq(0,125, by = 5))

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,420,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 2,
          tstart1 = 71,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          R0Dec1 = 0.8,
          R0Dec2 = 0.8)

outcomelist <- list()

for (z in 1:5) {
  
  parms["scen"] = z
  
  for(j in 1:nrow(lengthdata)) {
    j = j
    
    parms["t_dur1"] <- lengthdata[j,1] * 7
    parms["t_dur2"] <- lengthdata[j,2] * 7 
    
    outcomelist[[j]] <- local({
      
      optim <- data.frame(matrix(ncol = 6, nrow = nrow(optimdata)))
      
      for(i in 1:nrow(optimdata)) {
        parms["tstart1"] <- optimdata[i,1]
        parms["tstart2"] <- optimdata[i,2]
        
        out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
        optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                       "tstart1" = parms[["tstart1"]], "tstart2" = parms[["tstart2"]],
                       "realstart2" = parms[["tstart1"]] + parms[["t_dur1"]] + parms[["tstart2"]])
      }
      
      colnames(optim) <- c("peak", "cum", "scen", "tstart1", "tstart2", "realstart2")
      
      p1 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= peak))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + labs(x = "Intervention 1 Trigger", y = "Intervention 2 Trigger", fill = "Peak I(t)", 
                                                          title = paste0("Int 1 & 2 Duration: ", lengthdata[j,1], "/", lengthdata[j,2], " weeks")) + 
        scale_fill_viridis_c(direction = -1)
      
      p2 <- ggplot(optim, aes(x = tstart1, y = tstart2, fill= cum))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + 
        scale_fill_viridis_c(direction = -1, option = "magma") + labs(x = "Intervention 1 Trigger", y = "Intervention 2 Trigger", fill = "Cumulative\nIncidence", 
                                                                      title = paste0("Int 1 & 2 Duration: ", lengthdata[j,1], "/", lengthdata[j,2], " wks"))
      
      print(paste0("Scenario: ", z, " | Length Analysis: ", (j/nrow(lengthdata))*100, "%"))
      
      return(list(p1,p2))
    })   
  }
  
  peak_multilength <- ggarrange(outcomelist[[1]][[1]], outcomelist[[2]][[1]], outcomelist[[3]][[1]],
                                outcomelist[[4]][[1]], outcomelist[[5]][[1]], outcomelist[[6]][[1]],
                                outcomelist[[7]][[1]], outcomelist[[8]][[1]], outcomelist[[9]][[1]],
                                nrow = 3, ncol =3, align = "hv")
  
  ggsave(peak_multilength, filename = paste0("heatpeak_multilength",z,".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
  
  cum_multilength <- ggarrange(outcomelist[[1]][[2]], outcomelist[[2]][[2]], outcomelist[[3]][[2]],
                               outcomelist[[4]][[2]], outcomelist[[5]][[2]], outcomelist[[6]][[2]],
                               outcomelist[[7]][[2]], outcomelist[[8]][[2]], outcomelist[[9]][[2]],
                               nrow = 3, ncol =3, align = "hv")
  
  ggsave(cum_multilength, filename = paste0("heatcum_multilength", z, ".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
}


# Multiple Optimisations - R0 + Length --------------------------------------------------

lengthdata <- expand.grid("int1length" = seq(3,9, by = 3), "int2length" = seq(3,9, by = 3))
r0optim <- expand.grid("r01" = seq(0, 1.7, by = 0.1), "r01" = seq(0, 1.7, by = 0.1))

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,1000,by = 1)

parms = c(gamma = 1/GenTime(3.3, 2.8),
          scen = 2,
          tstart1 = 71,
          t_dur1 = 6*7,
          tstart2 = 42,
          t_dur2 = 6*7,
          R0Dec1 = 0.8,
          R0Dec2 = 0.8)

outcomelistr0 <- list()

for (z in 1:5) {
  
  parms["scen"] = z
  
  for(j in 1:nrow(lengthdata)) {
    j = j
    
    parms["t_dur1"] <- lengthdata[j,1] * 7
    parms["t_dur2"] <- lengthdata[j,2] * 7 
    
    outcomelistr0[[j]] <- local({
      
      optim <- data.frame(matrix(ncol = 5, nrow = nrow(optimdata)))
      
      for(i in 1:nrow(r0optim)) {
        parms["R0Dec1"] <- r0optim[i,1]
        parms["R0Dec2"] <- r0optim[i,2]
        
        out <- data.frame(ode(y = init, func = SIRmulti, times = times, parms = parms))
        optim[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                       "R0Dec1" = parms[["R0Dec1"]], "R0Dec2" = parms[["R0Dec2"]])
      }
      
      colnames(optim) <- c("peak", "cum", "scen", "R0Dec1", "R0Dec2")
      
      p1 <- ggplot(optim, aes(x = R0Dec1, y = R0Dec2, fill= peak))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + labs(x = bquote(.(Intervention ~ 1 ~ R[0])), y = bquote(.(Intervention ~ 2 ~ R[0])), fill = "Peak I(t)", 
                                                          title = paste0("Int 1 & 2 Duration: ", lengthdata[j,1], "/", lengthdata[j,2], " weeks")) + 
        scale_fill_viridis_c(direction = -1)
      
      p2 <- ggplot(optim, aes(x = R0Dec1, y = R0Dec2, fill= cum))  + geom_tile()  +
        scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
        theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
              axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
              legend.key.width =  unit(0.5, "cm")) + scale_fill_viridis_c(direction = -1, option = "magma") + 
        labs(x = bquote(.(Intervention ~ 1 ~ R[0])), y = bquote(.(Intervention ~ 2 ~ R[0])), fill = "Cumulative\nIncidence", 
                                                                      title = paste0("Int 1 & 2 Duration: ", lengthdata[j,1], "/", lengthdata[j,2], " wks"))
      
      print(paste0("Scenario: ", z, " | Length Analysis: ", (j/nrow(lengthdata))*100, "%"))
      
      return(list(p1,p2))
    })   
  }
  
  peak_multilength_r0 <- ggarrange(outcomelistr0[[1]][[1]], outcomelistr0[[2]][[1]], outcomelistr0[[3]][[1]],
                                outcomelistr0[[4]][[1]], outcomelistr0[[5]][[1]], outcomelistr0[[6]][[1]],
                                outcomelistr0[[7]][[1]], outcomelistr0[[8]][[1]], outcomelistr0[[9]][[1]],
                                nrow = 3, ncol =3, align = "hv")
  
  ggsave(peak_multilength_r0, filename = paste0("heatpeakr0_multilength",z,".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
  
  cum_multilength_r0 <- ggarrange(outcomelistr0[[1]][[2]], outcomelistr0[[2]][[2]], outcomelistr0[[3]][[2]],
                                  outcomelistr0[[4]][[2]], outcomelistr0[[5]][[2]], outcomelistr0[[6]][[2]],
                                  outcomelistr0[[7]][[2]], outcomelistr0[[8]][[2]], outcomelistr0[[9]][[2]],
                               nrow = 3, ncol =3, align = "hv")
  
  ggsave(cum_multilength_r0, filename = paste0("heatcumr0_multilength", z, ".png"), dpi = 300, type = "cairo", width = 16, height = 16, units = "in")
}


# Multiple Intervention Case Study - R0 Trajectory - Scen 1----------------------------------------

# Run the Model

r0optim <- expand.grid("r01" = seq(0.2, 1, by = 0.4), "r01" = seq(0.2, 1, by = 0.4))

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

# Multiple Intervention Case Study - R0 Trajectory - Scen 2 ----------------------------------------

# Run the Model

r0optim <- expand.grid("r01" = seq(0.2, 1, by = 0.4), "r01" = seq(0.2, 1, by = 0.4))

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


end_time <- Sys.time()
end_time - start_time