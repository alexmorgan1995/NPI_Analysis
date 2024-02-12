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

# Metapopulation Model ----------------------------------------------------

#ODE equations - SIR model with Beta(t) defined
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    beta <- beta_matrix
    
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
init <- c(S_i = 0.99999/3, I_i = 0.00001/3, V_i = 0, R_i = 0, C_i = 0,
          S_j = 0.99999/3, I_j = 0.00001/3, V_j = 0, R_j = 0, C_j = 0,
          S_k = 0.99999/3, I_k = 0.00001/3, V_k = 0, R_k = 0, C_k = 0)

times <- seq(0,200,by = 1)
parms = list(gamma = 1/GenTime(3, 2.8),
             beta_matrix = beta_base,
             eff = 0.9,
             P_i = 0.8,
             P_j = 0.7,
             P_k = 0.6,
             dt_i = 90,
             dt_j = 90,
             dt_k = 90,
             T_i = 20,
             T_j = 40,
             T_k = 60)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

data_i <- melt(out, id.vars = c("time"), measure.vars = c("I_i", "V_i", "R_i"), value.name = "I")
data_j <- melt(out, id.vars = c("time"), measure.vars = c("I_j", "V_j", "R_j"), value.name = "I")
data_k <- melt(out, id.vars = c("time"), measure.vars = c("I_k", "V_k", "R_k"), value.name = "I")

shade_i <- data.frame(xmin =  parms[["T_i"]], xmax = parms[["T_i"]]+(parms[["dt_i"]]), ymin = 0, ymax = Inf)

p_i <- ggplot(data = data_i, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.25), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen", "darkblue")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_i, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Population i")

shade_j <- data.frame(xmin =  parms[["T_j"]], xmax = parms[["T_j"]]+(parms[["dt_j"]]), ymin = 0, ymax = Inf)

p_j <- ggplot(data = data_j, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.3), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen", "darkblue")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_j, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_j"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Population j")

shade_k <- data.frame(xmin =  parms[["T_k"]], xmax = parms[["T_k"]]+(parms[["dt_k"]]), ymin = 0, ymax = Inf)

p_k <- ggplot(data = data_k, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.3), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen", "darkblue")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_k, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_k"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Population k")


# What happens if there is no vaccination ---------------------------------

parms1 <- parms
parms1[["dt_i"]] <- 0
parms1[["dt_j"]] <- 0
parms1[["dt_k"]] <- 0

out_noint <- cbind(data.frame(ode(y = init, func = SIR, times = times, parms = parms1)),
                   out$I_i,out$I_j, out$I_k)

data_i_noint <- melt(out_noint, id.vars = c("time"), measure.vars = c("I_i", "out$I_i"), value.name = "I")
data_j_noint <- melt(out_noint, id.vars = c("time"), measure.vars = c("I_j", "out$I_j"), value.name = "I")
data_k_noint <- melt(out_noint, id.vars = c("time"), measure.vars = c("I_k", "out$I_k"), value.name = "I")

p_noint_i <- ggplot(data = data_i_noint, aes(x = time, y = I, alpha = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.1), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
 scale_alpha_manual(values = c(0.4,1), labels = c("No Vaccination", "Vaccinated" )) + 
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18), legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_i, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity", color = "darkred") + 
  labs(x ="Time (Days)", y = "Prevalence", alpha = "Scenario", title = "Effects of Vaccination (i)") 

p_noint_j <- ggplot(data = data_j_noint, aes(x = time, y = I, alpha = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.1), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_alpha_manual(values = c(0.4,1), labels = c("No Vaccination", "Vaccinated" )) + 
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_j, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_j"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity", color = "darkred") + 
  labs(x ="Time (Days)", y = "Prevalence", alpha = "Scenario", title = "Effects of Vaccination (j)") 

p_noint_k <- ggplot(data = data_k_noint, aes(x = time, y = I, alpha = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.1), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_alpha_manual(values = c(0.4,1), labels = c("No Vaccination", "Vaccinated" )) + 
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) +
  geom_rect(data = shade_k, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_k"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity", color = "darkred") + 
  labs(x ="Time (Days)", y = "Prevalence", alpha = "Scenario", title = "Effects of Vaccination (k)") 


# Final plotting ----------------------------------------------------------

pcomb <- ggarrange(p_i, p_noint_i,NULL,
          p_j, p_noint_j,NULL,
          p_k, p_noint_k,NULL,
          nrow = 3, ncol = 3, widths = c(1,1,0.05))

ggsave(pcomb, filename = "vacc_scenario.png", dpi = 300, type = "cairo", width = 14, height = 12, units = "in")
