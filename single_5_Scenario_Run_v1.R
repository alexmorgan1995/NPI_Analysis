library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis")

# Model Functions ----------------------------------------------------------
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Need to rescale the model for c(t) and not beta(t) 

#FUnction to Obtain Betas for the Model 
combbeta <- function(scen, time, tstart, t_dur, cmin) {
  gamma <- 1/GenTime(3, 2.8)
  betascale <- (2.8*gamma)*0.7
  if(scen == 0) {
    output <- betascale
  }
  if(scen == 1) {
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     betascale*cmin, betascale)
  }
  if(scen == 2) {
    cminfun <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(cmin, 1), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     cminfun(time)*betascale, 
                     betascale)
  }
  if(scen == 3) {
    cminfun <- approxfun(x=c(tstart, (tstart + t_dur)), y= c(1, cmin), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     cminfun(time)*betascale, 
                     betascale)
  }
  if(scen == 4) {
    cminfun <-approxfun(x=c(tstart, tstart+(t_dur/2), tstart+t_dur), y= c(1, cmin, 1), method="linear", rule =2)
    output <- ifelse((time >= tstart & time <= tstart + t_dur),
                     cminfun(time)*betascale, betascale)
  }
  if(scen == 5) {
    output <- ifelse((time >= tstart & time <= tstart+(t_dur*0.16667)) | (time >= tstart+(t_dur*0.3333) & time <= tstart+(t_dur*0.5)) |
                       (time >= tstart+(t_dur*0.6667) & time <= tstart+(t_dur*0.83333)), 
                     betascale*cmin, betascale)
  }
  return(output)
}

plot(seq(0,400),combbeta(5, seq(0,400), 12, 24*7, 0.4))

#SIR set of ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, cmin)
    
    dS = - beta*S*(I)
    dI = beta*S*(I)- gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dI, dR, dC)))
  })
} 

# Baseline Model Analysis - 5 Scenarios -----------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

datalist <- list()

for(j in 1:length(seq(1,5))) {
  datalist[[j]] <- local({
    j=j
    data <- data.frame(matrix(nrow = 9, ncol = 0))
    explor_scen <- c(0, seq(1,5)[j])
    
    
    for(i in 1:2) {
      
      parms["scen"] <- explor_scen[i]
      
      if(parms["scen"] != 0 && parms["scen"] != 1) {
        parms["t_dur"] = 24*7
      }
      
      out <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], 
                   "beta" = combbeta(explor_scen[i], times, parms[["tstart"]], parms[["t_dur"]], parms[["cmin"]]))
      out$re <- out$beta/parms[["gamma"]]*out$S
      data <- rbind(data, out)
    }
    
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    plotbeta <- melt(data, id.vars = c("time", "group"), measure.vars = ("beta"))
    plotre <- melt(data, id.vars = c("time", "group"), measure.vars = ("re"))
    peak <- round(max(out$I), 3); cum <- round(max(out$C), 3)
    
    datatext <- data.frame(x = c(200, 200), y = c(0.1875, 0.163), label = c( paste0("italic(I)[italic(max)]", " ==", peak), 
                                                                             paste0("italic(I)[italic(c)](italic(t)[italic(max)])", " ==", cum)))
    
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group, alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0, 0.2),expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) + 
      scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkred", "darkred"))
    
    p2 <- ggplot(plotbeta, aes(x = time, y = value, col = group, alpha= group)) + theme_bw() + scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(limits = c(0 , 0.3), expand = c(0,0)) + geom_line(size = 1.1, stat = "identity") + 
      scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    p3 <- ggplot(plotre, aes(x = time, y = value, col = group, alpha= group)) + theme_bw() + scale_x_continuous(expand = c(0,0)) +
      geom_hline(yintercept = 1, size = 1.1, lty = 2, col = "black") + geom_line(size = 1.02, stat = "identity") + 
      scale_alpha_manual(values = c(0.35,1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    shade <- data.frame(xmin =  parms[["tstart"]], xmax = parms[["tstart"]]+(parms[["t_dur"]]), ymin = 0, ymax = Inf)
    
    if(parms[["scen"]]  == 1) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x =NULL, y = "Prevalence", col = "", title = "Scenario 1") +
        theme(axis.title.y=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
              axis.text.x=element_blank(),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      
      p2 <- p2 + theme(axis.title.y=element_text(size=18), axis.text.x=element_blank(), axis.text.y=element_text(size=15), 
                       plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x =NULL, y = expression(beta[(t)]), col = "") 
      
      p3 <- p3 + theme(axis.text=element_text(size=15), axis.title.y=element_text(size=18), axis.title.x = element_text(size=18), 
                       plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x ="Time (Days)", y = expression(R[e]), col = "")
    }
    
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x =NULL, y = NULL, col = "", title = "Scenario 2") +
        theme(plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
              axis.text.x=element_blank(), axis.text.y=element_blank(), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      
      p2 <- p2 + theme(axis.text=element_blank(), plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x =NULL, y = NULL, col = "") 
      
      p3 <- p3 + theme(axis.text.x =element_text(size=15), axis.text.y=element_blank(), axis.title.x = element_text(size=18), 
                       plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x ="Time (Days)", y = NULL, col = "")
    }
    
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x =NULL, y = NULL, col = "", title = "Scenario 3") +
        theme(plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
              axis.text.x=element_blank(), axis.text.y=element_blank(), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      
      p2 <- p2 + theme(axis.text=element_blank(), plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x =NULL, y = NULL, col = "") 
      
      p3 <- p3 + theme(axis.text.x =element_text(size=15), axis.text.y=element_blank(), axis.title.x = element_text(size=18), 
                       plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x ="Time (Days)", y = NULL, col = "")
    }
    
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x =NULL, y = NULL, col = "", title = "Scenario 3") +
        theme(plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
              axis.text.x=element_blank(), axis.text.y=element_blank(), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      
      p2 <- p2 + theme(axis.text=element_blank(), plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x =NULL, y = NULL, col = "") 
      
      p3 <- p3 + theme(axis.text.x =element_text(size=15), axis.text.y=element_blank(), axis.title.x = element_text(size=18), 
                       plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x ="Time (Days)", y = NULL, col = "")
    }
    
    if(parms[["scen"]]  == 5) {
      shade5 <- data.frame(xmin =  c(parms[["tstart"]], parms[["tstart"]]+(parms[["t_dur"]]*0.333), parms[["tstart"]]+(parms[["t_dur"]]*0.667)),
                           xmax = c(parms[["tstart"]]+(parms[["t_dur"]]*0.1667), parms[["tstart"]]+(parms[["t_dur"]]*0.53), 
                                    parms[["tstart"]]+(parms[["t_dur"]]*0.833)), ymin = 0, ymax = Inf)
      
      p1 <- p1 + geom_rect(data = shade5, inherit.aes = F, aes(ymin = ymin,  ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2, fill = "darkblue") + 
        geom_line(size = 1.1, stat = "identity")  + labs(x =NULL, y = NULL, col = "", title = "Scenario 5") + 
        theme(plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), axis.text.x=element_blank(),axis.text.y=element_blank(),
              plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      
      p2 <- p2 + theme(axis.text=element_blank(), plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x =NULL, y = NULL, col = "") 
      
      p3 <- p3 + theme(axis.text.x =element_text(size=15), axis.text.y=element_blank(), axis.title.x = element_text(size=18), 
                       plot.margin=unit(c(0.4,0.6,0.4,0),"cm")) + labs(x ="Time (Days)", y = NULL, col = "")
    }
    
    combplot <- ggarrange(p1,NULL, p2, NULL, p3, nrow = 5, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1,-0.02, 0.4,-0.02, 0.5))
    dump <- list(combplot, plotdata, max(out$C))
    return(dump)
  })
}

combplot <- ggarrange(NULL,datalist[[1]][[1]], datalist[[2]][[1]],  datalist[[3]][[1]], datalist[[4]][[1]], datalist[[5]][[1]], NULL, nrow = 1, ncol = 7, 
          legend = "none", widths = c(0.1,1.2,0.9,0.9,0.9,0.9,0.1))

for(i in 1:5) {
  print(datalist[[i]][[2]][datalist[[i]][[2]]$group == "scenario",]
        [which.max(datalist[[i]][[2]]$value[datalist[[i]][[2]]$group == "scenario"]),])
  print(datalist[[i]][[3]])
}

# Sensitivity Analysis ----------------------------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,800,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

sens <- list("tstart" = seq(1,150, by = 1),
             "cmin" = seq(0.01,1, by = 0.01),
             "t_dur" = seq(1,400, by = 1))

senslist <- list()

for(i in 1:3) {
  senslist[[i]] <- local({
    sensitivity <- sens[[i]]
    datasens <- data.frame(matrix(nrow = 0, ncol = 4))
    
    for(k in 1:5) {
      data <- cbind(data.frame(matrix(ncol = 4, nrow = length(sensitivity))), "group" = as.factor(seq(1,5)[k]))
      parms["scen"] <- seq(1,5)[k]
      
      if(parms["scen"] != 0 && parms["scen"] != 1 && i != 3) {
        parms["t_dur"] = 24*7
      }
      
      for(j in 1:length(sensitivity)) {
        parms[names(sens)[i]] <- sensitivity[j]
        #print(parms[names(sens)[i]])
        out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
        data[j,1] <- max(out$I)
        data[j,2] <- max(out$C)
        data[j,3] <- sensitivity[j]
        data[j,4] <- names(sens)[i]
      }
      datasens <- rbind(datasens, data)
      print(paste0("Sensitivity Test: ", i,", Scenario:", k))
    }
    
    colnames(datasens) <- c("peak", "cum", "sensval", "sens","group")
    
    p1 <- ggplot(datasens, aes(x = sensval, y = peak, col = group)) + geom_line(size = 1.02) + theme_bw() +
      scale_y_continuous(limits = c(0,0.150), expand = c(0,0))  
    
    p2 <- ggplot(datasens, aes(x = sensval, y = cum, col = group)) + geom_line(size = 1.02) + theme_bw() + 
      scale_y_continuous(limits = c(0,1),expand = c(0,0))
    
    if(names(sens)[i] == "tstart") {
      p1 <- p1 + labs(x = NULL, y = "I(t) Peak", col = "Scenario") + scale_x_continuous(limits = c(0,150),expand = c(0, 0))  + 
        theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text.x=element_blank(), axis.text.y=element_text(size=15),
              axis.title.y=element_text(size=18),legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.1,0.3,0.1,0.3),"cm"))
      p2 <- p2 + labs(x = bquote("Intervention Trigger ("*italic(t[p])*")"), y = "Cumulative Incidence", col = "Scenario") +
        theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.1,0.3,0.1,0.3),"cm")) +
      scale_x_continuous(limits = c(0,150), expand = c(0, 0))
    }
    
    if(names(sens)[i] == "cmin") {
      p1 <- p1 + labs(x = NULL, y = NULL, col = "Scenario") + scale_x_continuous(limits = c(0,1),expand = c(0, 0))  + 
        theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18), axis.text.x=element_blank(), axis.text.y=element_blank(),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.1,0.3,0.1,0.3),"cm"))
      p2 <- p2 + labs(x = bquote("Lockdown Scaling Factor ("*italic(c[min])*")"), y = NULL, col = "Scenario") +
        theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18), axis.text.x=element_text(size=15), axis.text.y=element_blank(),
              axis.title.x = element_text(size=18), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.1,0.3,0.1,0.3),"cm")) +
        scale_x_continuous(limits = c(0,1),expand = c(0, 0))
    }
    
    if(names(sens)[i] == "t_dur") {
      p1 <- p1 + labs(x = NULL, y = NULL, col = "Scenario") + scale_x_continuous(limits = c(0,400),expand = c(0, 0))  + 
        theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text.x=element_blank(), axis.text.y=element_blank(),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.1,0.3,0.1,0.3),"cm"))
      p2 <- p2 + labs(x = bquote("Intervention Duration ("*italic(d[t])*")"), y = NULL, col = "Scenario") + 
        theme(legend.position = "bottom", legend.title = element_text(size=18), legend.text=element_text(size=18),  axis.text.x=element_text(size=15), axis.text.y=element_blank(),
              axis.title.x = element_text(size=18), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.1,0.3,0.1,0.3),"cm")) +
        scale_x_continuous(limits = c(0,400),expand = c(0, 0)) 
    }
    dump <- list(p1,p2, datasens)
    return(dump)
  })
}

combsensplot <- ggarrange(NULL, senslist[[1]][[1]], NULL, senslist[[2]][[1]],  NULL, senslist[[3]][[1]],NULL,
          NULL, NULL, NULL, NULL,NULL,NULL, NULL,
          NULL, senslist[[1]][[2]], NULL, senslist[[2]][[2]],NULL, senslist[[3]][[2]],NULL,
          nrow = 3, ncol = 7, common.legend = TRUE,
          legend = "bottom", align = "h", heights = c(1,-0.12,1), 
          widths = c(0.1,1,0.04,0.9,0.04,0.9,0.1))

for(j in 1:3) {
  print(c("tstart1", "cmin", "t_dur")[j])
  for(i in 1:5) {
    print(senslist[[j]][[3]][senslist[[j]][[3]]$group == i,][which.max(senslist[[j]][[3]]$peak[senslist[[j]][[3]]$group == i]),])
  }
}

# CombPlot Test -----------------------------------------------------------

test <- ggarrange(NULL, combplot, NULL, combsensplot, nrow = 4, ncol = 1, labels =  c("","A","","B") ,
          font.label = c(size = 35), vjust = -0.7, hjust = -0.1, heights = c(0.1,1,0.1,0.9))
ggsave(test, filename = "sens_traj_combplot.png", dpi = 300, type = "cairo", width = 13, height = 14, units = "in")

# Multi-Parameter Optimisation --------------------------------------------

init <- c(S = 0.99999, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)

parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          t_dur = 12*7,
          cmin = 0.4)

parameterspace <- expand.grid("trigday" = seq(0,100, by =5), "length" = seq(1,250, by =5))

scensens <- list()

for(j in 1:5) {
  
  scensens[[j]] = local({
    i = 0
    scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 5))
    
    for(i in 1:nrow(parameterspace)) {
      
      print(paste0("Scenario ", j," - ", round(i/nrow(parameterspace), digits = 2)))
      parms["tstart"] <- parameterspace[i,1]
      parms["t_dur"] <- parameterspace[i,2] 
      parms["scen"] <- j
      
      out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
      scendata[i,] <- c("peak" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                        "tstart" = parms[["tstart"]], "t_dur" = parms[["t_dur"]])
    }
    
    colnames(scendata) <- c("peak", "cum", "scen", "tstart", "t_dur")
    
    p1 <- ggplot(scendata, aes(x = tstart, y = t_dur, fill= peak))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + 
      labs(x = bquote("Trigger Date ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "I(t) Peak", title = paste("Scenario", j)) + 
      scale_fill_viridis_c(direction = -1) #, breaks= seq(0.01, 0.15, by = (0.15-0.01)/4), limits = c(0.01, 0.15)
    
    p2<- ggplot(scendata, aes(x = tstart, y = t_dur, fill = cum))  + geom_tile() +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 2, hjust = -0.2),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.3,0.4,0.3,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + 
      labs(x = bquote("Trigger Date ("*italic(t[p])*")"), y = bquote("Duration ("*italic(d[t])*")"), fill = "Total\nCumulative\nIncidence", title = "") + 
      scale_fill_viridis_c(direction = -1, option = "magma") #, breaks= seq(0.1, 0.8, by = (0.8-0.1)/4), limits = c(0.1, 0.8) 
    
    combplot <- ggarrange(p1,p2, ncol = 2, nrow = 1, widths = c(1,1.05), align = "h")
    print(combplot)
    dump <- list(combplot, scendata)
    return(dump)
  })
}

for(j in 1:5) {
    print(scensens[[j]][[2]]$tstart[which.min(scensens[[j]][[2]]$peak)])
}

combplot <- ggarrange(scensens[[1]][[1]],scensens[[2]][[1]],scensens[[3]][[1]],scensens[[4]][[1]],scensens[[5]][[1]],
                      nrow = 5, ncol = 1)

ggsave(combplot, filename = "Heat_5_scenarios_sensitivity.png", dpi = 300, type = "cairo", width = 10, height = 16, units = "in")
