library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures")

# WAIFW Matrix ------------------------------------------------------------

#Function for the generation time - a function of R0 and the doubling time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

gamma <- 1/GenTime(3, 2.8)
beta <- (2.8*gamma)
1/2.8


#Beta_BASELINE
beta_base <- matrix(data = beta, nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k"))) 
beta_base_HALF <- beta_base; beta_base_HALF[c(3,7)] <- beta_base_HALF[c(3,7)]/2

#Beta_ALL

scale_ALL <- c(0.1,0.1,0.1,
               0.1,0.1,0.1,
               0.1,0.1,0.1)
beta_ALL <- beta_base*scale_ALL
beta_ALL_HALF <- beta_ALL; beta_ALL_HALF[c(3,7)] <- beta_ALL_HALF[c(3,7)]/2

#Just_K

scale_K <- c(1,1,0.1,
             1,1,0.1,
             0.1,0.1,0.1)
beta_K <- beta_base*scale_K
beta_K_HALF <- beta_K; beta_K_HALF[c(3,7)] <- beta_K_HALF[c(3,7)]/2

#Just_JK

scale_JK <- c(1,0.1,0.1,
              0.1,0.1,0.1,
              0.1,0.1,0.1)
beta_JK <- beta_base*scale_JK
beta_JK_HALF <- beta_JK; beta_JK_HALF[c(3,7)] <- beta_JK_HALF[c(3,7)]/2


# Time-Varying Beta Function ----------------------------------------------

gamma <- 1/GenTime(3, 2.8)
beta <- (2.8*gamma)
1/2.8


#Beta_BASELINE
beta_base <- matrix(data = beta, nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k"))) 
beta_base_HALF <- beta_base; beta_base_HALF[c(3,7)] <- beta_base_HALF[c(3,7)]/2

#Beta_BASELINE
beta_base <- matrix(data = beta, nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k"))) 
beta_base_HALF <- beta_base; beta_base_HALF[c(3,7)] <- beta_base_HALF[c(3,7)]/2

#Beta_ALL

scale_ALL <- c(0.1,0.1,0.1,
               0.1,0.1,0.1,
               0.1,0.1,0.1)
beta_ALL <- beta_base*scale_ALL
beta_ALL_HALF <- beta_ALL; beta_ALL_HALF[c(3,7)] <- beta_ALL_HALF[c(3,7)]/2

#if t == time when intervention starts - then set beta matrix to have lower value 

betafunc <- function(time, T_i, dt_i, T_j, dt_j, T_k, dt_k, socdist) {
  beta <- matrix(data = (1.5*gamma), nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k"))) 
  
  if(time > (T_i + dt_i)) {
    beta[1,] <- beta[1,]*2.8
    beta[2:3,1] <- beta[2:3,1]*2.8
  }
  
  if(time > (T_j + dt_j)) {
    beta[2,2:3] <- beta[2,2:3]*2.8
    beta[3,2] <- beta[3,2]*2.8
  }
  
  if(time > (T_k + dt_k)) {
    beta[3,3] <- beta[3,3]*2.8
  }
  
  if(socdist == "No") {
    beta <- matrix(data = (2.8*gamma), nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k")))
  }  
  
  return(beta)
}


# Metapopulation Model ----------------------------------------------------

#ODE equations - SIR model with Beta(t) defined
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    beta <- betafunc(time, T_i, dt_i, T_j, dt_j, T_k, dt_k, socdist)
    r_i <- ifelse(time < T_i | time > T_i + dt_i, 0, (-log(1-P_i)/dt_i))
    
    dS_i = - beta["i","i"]*S_i*(I_i) - beta["i","j"]*S_i*(I_j) - beta["i","k"]*S_i*(I_k) - r_i*S_i*eff      
    dI_i = beta["i","i"]*S_i*(I_i) + beta["i","j"]*S_i*(I_j) + beta["i","k"]*S_i*(I_k) - gamma*I_i
    dV_i = r_i*S_i*eff    
    dR_i = gamma*I_i
    dC_i = beta["i","i"]*S_i*(I_i) + beta["i","j"]*S_i*(I_j) + beta["i","k"]*S_i*(I_k)
    
    r_j <- ifelse(time < T_j | time > T_j + dt_j, 0, (-log(1-P_j)/dt_j))
    
    dS_j = - beta["j","i"]*S_j*(I_i) - beta["j","j"]*S_j*(I_j) - beta["j","k"]*S_j*(I_k) - r_j*S_j*eff  
    dI_j = beta["j","i"]*S_j*(I_i) + beta["j","j"]*S_j*(I_j) + beta["j","k"]*S_j*(I_k) - gamma*I_j
    dV_j = r_j*S_j*eff   
    dR_j = gamma*I_j 
    dC_j = beta["j","i"]*S_j*(I_i) + beta["j","j"]*S_j*(I_j) + beta["j","k"]*S_j*(I_k) 

    r_k <- ifelse(time < T_k | time > T_k + dt_k, 0, (-log(1-P_k)/dt_k))
    
    dS_k = - beta["k","i"]*S_k*(I_i) - beta["k","j"]*S_k*(I_j) - beta["k","k"]*S_k*(I_k) - r_k*S_k*eff  
    dI_k = beta["k","i"]*S_k*(I_i) + beta["k","j"]*S_k*(I_j) + beta["k","k"]*S_k*(I_k) - gamma*I_k
    dV_k = r_k*S_k*eff 
    dR_k = gamma*I_k
    dC_k = beta["k","i"]*S_k*(I_i) + beta["k","j"]*S_k*(I_j) + beta["k","k"]*S_k*(I_k)
    
    return(list(c(dS_i, dI_i, dV_i, dR_i, dC_i,
                  dS_j, dI_j, dV_j, dR_j, dC_j,
                  dS_k, dI_k, dV_k, dR_k, dC_k)))
  })
} 

# Epidemic Trajectory Plots for the 5 NPI scenarios -----------------------------------


#Initial Conditions and Parameters
init <- c(S_i = (1-(0.0071+0.1))/3, I_i = 0.0071/3, V_i = 0, R_i = 0.1/3, C_i = 0,
          S_j = (1-(0.0071+0.1))/3, I_j = 0.0071/3, V_j = 0, R_j = 0.1/3, C_j = 0,
          S_k = (1-(0.0071+0.1))/3, I_k = 0.0071/3, V_k = 0, R_k = 0.1/3, C_k = 0)

times <- seq(0,365,by = 1)
parms = list(gamma = 1/GenTime(3, 2.8),
             beta_matrix = beta_base*0.357,
             eff = 0.9,
             P_i = 0.9,
             P_j = 0.9,
             P_k = 0.9,
             dt_i = 90,
             dt_j = 90,
             dt_k = 90,
             T_i = 0,
             T_j = 90,
             T_k = 180,
             socdist = "Yes")

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

data_i <- melt(out, id.vars = c("time"), measure.vars = c("I_i", "V_i", "R_i"), value.name = "I")
data_j <- melt(out, id.vars = c("time"), measure.vars = c("I_j", "V_j", "R_j"), value.name = "I")
data_k <- melt(out, id.vars = c("time"), measure.vars = c("I_k", "V_k", "R_k"), value.name = "I")

#Plot i
shade_i <- data.frame(xmin =  parms[["T_i"]], xmax = parms[["T_i"]]+(parms[["dt_i"]]), ymin = 0, ymax = Inf)

p_i <- ggplot(data = data_i, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.5), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen", "darkblue")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_i, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Population i")

p_i2 <- ggplot(data = data_i[data_i$variable == "I_i",], aes(x = time, y = I)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.005), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=8), plot.margin=unit(c(0.2,0,0,0),"cm")) +
  geom_rect(data = shade_i, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity", color = "darkred") + 
  labs(x ="", y = "")

p_iplot <- p_i + annotation_custom(ggplotGrob(p_i2), xmin = 200, xmax = 360, 
                       ymin = 0.25, ymax = 0.5)

#Plot j
shade_j <- data.frame(xmin =  parms[["T_j"]], xmax = parms[["T_j"]]+(parms[["dt_j"]]), ymin = 0, ymax = Inf)

p_j <- ggplot(data = data_j, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.5), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen", "darkblue")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_j, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_j"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Population j")

p_j2 <- ggplot(data = data_j[data_j$variable == "I_j",], aes(x = time, y = I)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.01), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=8), plot.margin=unit(c(0.2,0,0,0),"cm")) +
  geom_rect(data = shade_j, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity", color = "darkred") + 
  labs(x ="", y = "")

p_jplot <- p_j + annotation_custom(ggplotGrob(p_j2), xmin = 200, xmax = 360, 
                                   ymin = 0.25, ymax = 0.5)

#Plot k
shade_k <- data.frame(xmin =  parms[["T_k"]], xmax = parms[["T_k"]]+(parms[["dt_k"]]), ymin = 0, ymax = Inf)

p_k <- ggplot(data = data_k, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.5), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen", "darkblue")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_k, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_k"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Population k")

p_k2 <- ggplot(data = data_k[data_k$variable == "I_k",], aes(x = time, y = I)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.01), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=8), plot.margin=unit(c(0.2,0,0,0),"cm")) +
  geom_rect(data = shade_k, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity", color = "darkred") + 
  labs(x ="", y = "")

p_kplot <- p_k + annotation_custom(ggplotGrob(p_k2), xmin = 200, xmax = 360, 
                                   ymin = 0.25, ymax = 0.5)

# What happens if there is no vaccination -Entire Population ---------------------------------

parms1 <- parms
parms1[["eff"]] <- 0

parms2 <- parms1
parms2[["socdist"]] <- "No"

outvac <- out
out_int <- data.frame(ode(y = init, func = SIR, times = times, parms = parms1))
out_noint <- data.frame(ode(y = init, func = SIR, times = times, parms = parms2))

combframe <- data.frame(cbind("time" = out$time, "vac_inf" = (out$I_i + out$I_j + out$I_k),
                              "novac_inf" = (out_int$I_i + out_int$I_j + out_int$I_k), 
                              "noint_inf" = (out_noint$I_i + out_noint$I_j + out_noint$I_k)))

data_i_comb <- melt(combframe, id.vars = c("time"), measure.vars = c("vac_inf", "novac_inf", "noint_inf"), value.name = "I")

shade_i <- data.frame(xmin =  c(parms[["T_i"]],parms[["T_j"]],parms[["T_k"]]), 
                      xmax = c(parms[["T_i"]]+(parms[["dt_i"]]),
                               parms[["T_j"]]+(parms[["dt_j"]]),
                               parms[["T_k"]]+(parms[["dt_k"]])), ymin = c(0,0,0), ymax = c(Inf,Inf,Inf))

datatext <- data.frame(x = c(182.5, 182.5, 182.5), y = c(0.24, 0.23, 0.22), 
                       label = c(paste0("Attack Rate (Vac/Int)", " = ", round(tail(out$C_i + out$C_j + out$C_k, 1), digits =3)), 
                                 paste0("Attack Rate (No Vac/Int)", " = ", round(tail(out_int$C_i + out_int$C_j + out_int$C_k, 1), digits =3)),
                                 paste0("Attack Rate (No Vac/No Int)", " = ", round(tail(out_noint$C_i + out_noint$C_j + out_noint$C_k, 1), digits =3))))

p_noint_i <- ggplot(data = data_i_comb, aes(x = time, y = I, lty = variable, alpha = variable, col = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.25), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
 scale_alpha_manual(name = "Scenario", 
                    values = c(1,0.4,0.4), labels = c("Vac + Int", "No Vac + Int" , "No Vac + No Int")) + 
  scale_linetype_manual(name = "Scenario",
                        values = c(1,1,2), labels = c("Vac + Int", "No Vac + Int" , "No Vac + No Int")) +
  scale_color_manual(name = "Scenario",
                        values = c("darkred","darkred","darkred"), labels = c("Vac + Int", "No Vac + Int" , "No Vac + No Int")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18), legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_i, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = c(0.4,0.2,0.4),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Prevalence", alpha = "Scenario", title = "Effects of Vaccination (ALL)") + 
  geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, col = "black", fontface = "bold", fill = "white")

# Final plotting ----------------------------------------------------------

pcomb <- ggarrange(p_iplot, NULL,
                   p_jplot, NULL,
                   p_kplot,NULL,
                   nrow = 3, ncol = 2, widths = c(1,0.05))

ggsave(pcomb, filename = "vacc_scenario_seq.png", dpi = 300, type = "cairo", width = 10, height = 14, units = "in")

ggsave(p_noint_i, filename = "total_vacc_scenario_seq.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in")
