library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/WriteUpAnalysis/supplementary")

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

#SEIR set of ODEs
SEIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- combbeta(scen, time, tstart, t_dur, cmin)
    
    dS = - beta*S*(I)
    dE = beta*S*(I)- sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I 
    
    dC = beta*S*I
    return(list(c(dS, dE, dI, dR, dC)))
  })
} 

# Baseline Model Analysis - 5 Scenarios -----------------------------------

init <- c(S = 0.99999, E = 0, I = 0.00001, R = 0, C = 0)
times <- seq(0,400,by = 1)
parms = c(gamma = 1/GenTime(3, 2.8),
          scen = 0,
          tstart = 52,
          sigma = 1/3,
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
      
      out <- cbind(data.frame(ode(y = init, func = SEIR, times = times, parms = parms)), 
                   "group" =  c("baseline", "scenario")[i], 
                   "beta" = combbeta(explor_scen[i], times, parms[["tstart"]], parms[["t_dur"]], parms[["cmin"]]))
      out$re <- out$beta/parms[["gamma"]]*out$S
      data <- rbind(data, out)
    }
 
    plotdata <- melt(data, id.vars = c("time", "group"), measure.vars = ("I"))
    plotbeta <- melt(data, id.vars = c("time", "group"), measure.vars = ("beta"))
    plotre <- melt(data, id.vars = c("time", "group"), measure.vars = ("re"))
    peak <- round(max(out$I), 3)
    cum <- round(max(out$C), 3)
    
    datatext <- data.frame(x = c(200, 200), y = c(0.15, 0.135), label = c( paste0("italic(I)[italic(max)]", " ==", peak), 
                                                                            paste0("italic(I)[italic(c)](italic(t)[italic(max)])", " ==", cum)))
    
    p1 <- ggplot(data = plotdata, aes(x = time, y = value, color = group , alpha= group)) + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.16),expand = c(0,0)) + 
      scale_x_continuous( expand = c(0, 0)) + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkred", "darkred"))

    p2 <- ggplot(plotbeta, aes(x = time, y = value, col = group, alpha= group))  + theme_bw() +
      scale_y_continuous(limits = c(0 , 0.3), expand = c(0,0)) + 
      geom_line(size = 1.1, stat = "identity") + scale_alpha_manual(values = c(0.35, 1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    p3 <- ggplot(plotre, aes(x = time, y = value, col = group, alpha= group)) + theme_bw() + 
      geom_hline(yintercept = 1, size = 1.1, lty = 2, col = "black") + geom_line(size = 1.02, stat = "identity") + 
      scale_alpha_manual(values = c(0.35,1)) + scale_color_manual(values = c("darkblue", "darkblue"))
    
    shade <- data.frame(xmin =  parms[["tstart"]], xmax = parms[["tstart"]]+(parms[["t_dur"]]), ymin = 0, ymax = Inf)
    
    if(parms[["scen"]]  == 1) {
     p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                          fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="Time (Days)", y = "Prevalence", col = "", title = "Scenario 1") +
       theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
             plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x=element_blank(),axis.text.y=element_text(size=15),
             axis.title.y=element_text(size=18), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
       geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
     p2 <- p2 + labs(x ="", y = expression(beta[(t)]), col = "") +
       theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
             axis.title.y=element_text(size=18),axis.title.x = element_blank(), axis.text.x=element_blank(),axis.text.y=element_text(size=15),
             legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) 
     p3 <- p3 + labs(x ="Time (Days)", y = expression(R[e]), col = "")  +
       theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text=element_text(size=15),
             axis.title.y=element_text(size=18),axis.title.x = element_text(size=18), 
             legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) +
       scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) 
    }
    if(parms[["scen"]] == 2) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 2") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),
              axis.title.y=element_text(size=18), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "")+
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              axis.title.y=element_text(size=18),axis.title.x = element_blank(), axis.text.x=element_blank(),axis.text.y=element_blank(),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) 
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  axis.text.x=element_text(size=15),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18), axis.text.y=element_blank(),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) 
    }
    if(parms[["scen"]]  == 3) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 3") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),
              axis.title.y=element_text(size=18), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              axis.title.y=element_text(size=18),axis.title.x = element_blank(), axis.text.x=element_blank(),axis.text.y=element_blank(),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) 
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),axis.text.x=element_text(size=15),axis.text.y=element_blank(),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18), 
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) 
    }
    if(parms[["scen"]]  == 4) {
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.2,
                           fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + labs(x ="", y = "", col = "", title = "Scenario 4")  +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x =element_blank(),axis.text.y=element_blank(),
              axis.title.y=element_text(size=18), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              axis.title.y=element_text(size=18),axis.title.x = element_blank(), axis.text.x=element_blank(),axis.text.y=element_blank(),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) 
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),axis.text.x=element_text(size=15),axis.text.y=element_blank(),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18), 
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) 
    }
    if(parms[["scen"]]  == 5) {
      
      shade <- data.frame(xmin =  c(parms[["tstart"]], parms[["tstart"]]+(parms[["t_dur"]]*0.333), 
                             parms[["tstart"]]+(parms[["t_dur"]]*0.667)),
                   xmax = c(parms[["tstart"]]+(parms[["t_dur"]]*0.1667), parms[["tstart"]]+(parms[["t_dur"]]*0.53), 
                            parms[["tstart"]]+(parms[["t_dur"]]*0.833)),
                   ymin = 0, ymax = Inf)
      p1 <- p1 + geom_rect(data = shade, inherit.aes = F, aes(ymin = ymin,  ymax = ymax, 
                     xmin = xmin, xmax = xmax), alpha = 0.2, fill = "darkblue") + geom_line(size = 1.1, stat = "identity") +
      labs(x ="", y = "", col = "",  title = "Scenario 5") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),
              axis.title.y=element_text(size=18), axis.title.x = element_blank(), legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
        geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", parse = TRUE, fontface = "bold", fill = "white")
      p2 <- p2 + labs(x ="", y = "", col = "") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),  
              axis.title.y=element_text(size=18),axis.title.x = element_blank(), axis.text.x=element_blank(),axis.text.y=element_blank(),
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) 
      p3 <- p3 + labs(x ="Time (Days)", y = "", col = "") +
        theme(legend.position = "bottom", legend.title = element_text(size=15), legend.text=element_text(size=18),axis.text.x=element_text(size=15),axis.text.y=element_blank(),
              axis.title.y=element_text(size=18),axis.title.x = element_text(size=18), 
              legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0 , 2), expand = c(0,0)) 
    }
    
    combplot <- ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "none", align = "v",heights = c(1, 0.4, 0.5))
    dump <- list(combplot, plotdata, max(out$C))
    return(dump)
  })
}

combplot <- ggarrange(NULL,datalist[[1]][[1]], datalist[[2]][[1]], datalist[[3]][[1]], datalist[[4]][[1]], datalist[[5]][[1]], NULL,nrow = 1, ncol = 7, 
          legend = "none", align = "h", widths = c(0.1,1,0.9,0.9,0.9,0.9,0.1))


for(i in 1:5) {
print(datalist[[i]][[2]][datalist[[i]][[2]]$group == "scenario",]
      [which.max(datalist[[i]][[2]]$value[datalist[[i]][[2]]$group == "scenario"]),])
print(datalist[[i]][[3]])
}

ggsave(combplot, filename = "5_scenarios_SEIR.png", dpi = 300, type = "cairo", width = 14, height = 10, units = "in")

