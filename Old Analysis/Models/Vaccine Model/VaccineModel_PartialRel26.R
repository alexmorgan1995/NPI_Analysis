library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Vaccine_Figures")

#Function for the generation time - a function of R0 and the doubling time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

gamma <- 1/GenTime(3, 2.8)

# Time-Varying Beta Function ----------------------------------------------

#if t == time when intervention starts - then set beta matrix to have lower value 

betafunc <- function(time, T_i, dt_i, T_j, dt_j, T_k, dt_k, socdist) {
  beta <- matrix(data = (1*gamma), nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k"))) 
  
  if(time > (T_i + dt_i)) {
    beta[1,] <- beta[1,]*2.6
    beta[2:3,1] <- beta[2:3,1]*2.6
  }
  
  if(time > (T_j + dt_j)) {
    beta[2,2:3] <- beta[2,2:3]*2.6
    beta[3,2] <- beta[3,2]*2.6
  }
  
  if(time > (T_k + dt_k)) {
    beta[3,3] <- beta[3,3]*2.6
  }
  
  if(socdist == "No") {
    beta <- matrix(data = (2.6*gamma), nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k")))
  }  
  
  return(beta)
}


# Simplified --------------------------------------------------------------

betafuncsimp <- function(time, timing) {
  beta <- matrix(data = (1*gamma), nrow = 3, ncol = 3, dimnames = list(c("i", "j", "k"), c("i", "j", "k"))) 
  
  if(time > timing) {
    beta <- beta*2.6
  }
  
  return(beta)
}

# Metapopulation Model ----------------------------------------------------

#ODE equations - SIR model with Beta(t) defined
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    ifelse(simp == "No", 
           beta <- betafunc(time, T_i, dt_i, T_j, dt_j, T_k, dt_k, socdist),
           beta <- betafuncsimp(time, timing))
    
    #Subpop i
    
    deltar_i <- ifelse(time < T_i | time > T_i + dt_i | vac == "No", 0, r_i)
    
    dS_i = sigma1*(V_i+R_i) + sigma2*(Rv_i+Vr_i) - 
      beta["i","i"]*S_i*(I_i) - beta["i","j"]*S_i*(I_j) - beta["i","k"]*S_i*(I_k) - 
      (1-eff2)*beta["i","i"]*S_i*(Vi_i) - (1-eff2)*beta["i","j"]*S_i*(Vi_j) - (1-eff2)*beta["i","k"]*S_i*(Vi_k) - 
      deltar_i*(S_i/(S_i+I_i+R_i))  
    
    dI_i = beta["i","i"]*S_i*(I_i) + beta["i","j"]*S_i*(I_j) + beta["i","k"]*S_i*(I_k) + 
      (1-eff2)*beta["i","i"]*S_i*(Vi_i) + (1-eff2)*beta["i","j"]*S_i*(Vi_j) + (1-eff2)*beta["i","k"]*S_i*(Vi_k) - 
      gamma*I_i - deltar_i*(I_i/(S_i+I_i+R_i))
    
    dR_i = gamma*I_i - deltar_i*(R_i/(S_i+I_i+R_i)) - sigma1*R_i
    
    dRv_i = deltar_i*((I_i/(S_i+I_i+R_i)) + (R_i/(S_i+I_i+R_i))) - sigma2*Rv_i
    
    dV_i = deltar_i*(S_i/(S_i+I_i+R_i)) -
      (1-eff1)*beta["i","i"]*V_i*(I_i) - (1-eff1)*beta["i","j"]*V_i*(I_j) - (1-eff1)*beta["i","k"]*V_i*(I_k) - 
      (1-eff1)*(1-eff2)*beta["i","i"]*V_i*(Vi_i) - (1-eff1)*(1-eff2)*beta["i","j"]*V_i*(Vi_j) - (1-eff1)*(1-eff2)*beta["i","k"]*V_i*(Vi_k) - 
      sigma1*V_i
    
    dVi_i = (1-eff1)*beta["i","i"]*V_i*(I_i) + (1-eff1)*beta["i","j"]*V_i*(I_j) + (1-eff1)*beta["i","k"]*V_i*(I_k) + 
      (1-eff1)*(1-eff2)*beta["i","i"]*V_i*(Vi_i) + (1-eff1)*(1-eff2)*beta["i","j"]*V_i*(Vi_j) + (1-eff1)*(1-eff2)*beta["i","k"]*V_i*(Vi_k) - 
      gamma*Vi_i
    
    dVr_i = gamma*Vi_i - sigma2*Vr_i
    
    dC_i = beta["i","i"]*S_i*(I_i) + beta["i","j"]*S_i*(I_j) + beta["i","k"]*S_i*(I_k) + 
      (1-eff2)*beta["i","i"]*S_i*(Vi_i) + (1-eff2)*beta["i","j"]*S_i*(Vi_j) + (1-eff2)*beta["i","k"]*S_i*(Vi_k) + 
      (1-eff1)*beta["i","i"]*V_i*(I_i) + (1-eff1)*beta["i","j"]*V_i*(I_j) + (1-eff1)*beta["i","k"]*V_i*(I_k) + 
      (1-eff1)*(1-eff2)*beta["i","i"]*V_i*(Vi_i) + (1-eff1)*(1-eff2)*beta["i","j"]*V_i*(Vi_j) + (1-eff1)*(1-eff2)*beta["i","k"]*V_i*(Vi_k)
    
    
    dCi_i = beta["i","i"]*S_i*(I_i) + beta["i","j"]*S_i*(I_j) + beta["i","k"]*S_i*(I_k) + 
      (1-eff2)*beta["i","i"]*S_i*(Vi_i) + (1-eff2)*beta["i","j"]*S_i*(Vi_j) + (1-eff2)*beta["i","k"]*S_i*(Vi_k)
    
    #Subpop j
    
    deltar_j <- ifelse(time < T_j | time > T_j + dt_j | vac == "No", 0, r_j)
    #r_j <- ifelse(time < T_j | time > T_j + dt_j | vac == "No", 0, (P_i)/dt_i)
    
    dS_j = sigma1*(V_j+R_j) + sigma2*(Rv_j+Vr_j) - 
      beta["j","i"]*S_j*(I_i) - beta["j","j"]*S_j*(I_j) - beta["j","k"]*S_j*(I_k) - 
      (1-eff2)*beta["j","i"]*S_j*(Vi_i) - (1-eff2)*beta["j","j"]*S_j*(Vi_j) - (1-eff2)*beta["j","k"]*S_j*(Vi_k) - 
      deltar_j*(S_j/(S_j+I_j+R_j))  
    
    dI_j = beta["j","i"]*S_j*(I_i) + beta["j","j"]*S_j*(I_j) + beta["j","k"]*S_j*(I_k) + 
      (1-eff2)*beta["j","i"]*S_j*(Vi_i) + (1-eff2)*beta["j","j"]*S_j*(Vi_j) + (1-eff2)*beta["j","k"]*S_j*(Vi_k) - 
      gamma*I_j - deltar_j*(I_j/(S_j+I_j+R_j))
    
    dR_j = gamma*I_j - deltar_j*(R_j/(S_j+I_j+R_j)) - sigma1*R_j
    
    dRv_j = deltar_j*((I_j/(S_j+I_j+R_j)) + (R_j/(S_j+I_j+R_j))) - sigma2*Rv_j
    
    dV_j = deltar_j*(S_j/(S_j+I_j+R_j)) -
      (1-eff1)*beta["j","i"]*V_j*(I_i) - (1-eff1)*beta["j","j"]*V_j*(I_j) - (1-eff1)*beta["i","k"]*V_j*(I_k) - 
      (1-eff1)*(1-eff2)*beta["j","i"]*V_j*(Vi_i) - (1-eff1)*(1-eff2)*beta["j","j"]*V_j*(Vi_j) - (1-eff1)*(1-eff2)*beta["j","k"]*V_j*(Vi_k) - 
      sigma1*V_j
    
    dVi_j = (1-eff1)*beta["j","i"]*V_j*(I_i) + (1-eff1)*beta["j","j"]*V_j*(I_j) + (1-eff1)*beta["j","k"]*V_j*(I_k) + 
      (1-eff1)*(1-eff2)*beta["j","i"]*V_j*(Vi_i) + (1-eff1)*(1-eff2)*beta["j","j"]*V_j*(Vi_j) + (1-eff1)*(1-eff2)*beta["j","k"]*V_j*(Vi_k) - 
      gamma*Vi_j
    
    dVr_j = gamma*Vi_j - sigma2*Vr_j
    
    dC_j = beta["j","i"]*S_j*(I_i) + beta["j","j"]*S_j*(I_j) + beta["j","k"]*S_j*(I_k) + 
      (1-eff2)*beta["j","i"]*S_j*(Vi_i) + (1-eff2)*beta["j","j"]*S_j*(Vi_j) + (1-eff2)*beta["j","k"]*S_j*(Vi_k) +
      (1-eff1)*beta["j","i"]*V_j*(I_i) + (1-eff1)*beta["j","j"]*V_j*(I_j) + (1-eff1)*beta["j","k"]*V_j*(I_k) + 
      (1-eff1)*(1-eff2)*beta["j","i"]*V_j*(Vi_i) + (1-eff1)*(1-eff2)*beta["j","j"]*V_j*(Vi_j) + (1-eff1)*(1-eff2)*beta["j","k"]*V_j*(Vi_k)
    
    #Subpop k
    deltar_k <- ifelse(time < T_k | time > T_k + dt_k | vac == "No", 0, r_k)
    #r_k <- ifelse(time < T_k | time > T_k + dt_k | vac == "No", 0, (P_i)/dt_i)
    
    dS_k = sigma1*(V_k+R_k) + sigma2*(Rv_k+Vr_k) - 
      beta["k","i"]*S_k*(I_i) - beta["k","j"]*S_k*(I_j) - beta["k","k"]*S_k*(I_k) - 
      (1-eff2)*beta["k","i"]*S_k*(Vi_i) - (1-eff2)*beta["k","j"]*S_k*(Vi_j) - (1-eff2)*beta["k","k"]*S_k*(Vi_k) - 
      deltar_k*(S_k/(S_k+I_k+R_k))  
    
    dI_k = beta["k","i"]*S_k*(I_i) + beta["k","j"]*S_k*(I_j) + beta["k","k"]*S_k*(I_k) + 
      (1-eff2)*beta["k","i"]*S_k*(Vi_i) + (1-eff2)*beta["k","j"]*S_k*(Vi_j) + (1-eff2)*beta["k","k"]*S_k*(Vi_k) - 
      gamma*I_k - deltar_k*(I_k/(S_k+I_k+R_k))
    
    dR_k = gamma*I_k - deltar_k*(R_k/(S_k+I_k+R_k)) - sigma1*R_k
    
    dRv_k = deltar_k*((I_k/(S_k+I_k+R_k)) + (R_k/(S_k+I_k+R_k))) - sigma2*Rv_k
    
    dV_k = deltar_k*(S_k/(S_k+I_k+R_k)) -
      (1-eff1)*beta["k","i"]*V_k*(I_i) - (1-eff1)*beta["k","j"]*V_k*(I_j) - (1-eff1)*beta["k","k"]*V_k*(I_k) - 
      (1-eff1)*(1-eff2)*beta["k","i"]*V_k*(Vi_i) - (1-eff1)*(1-eff2)*beta["k","j"]*V_k*(Vi_j) - (1-eff1)*(1-eff2)*beta["k","k"]*V_k*(Vi_k) - 
      sigma1*V_k
    
    dVi_k = (1-eff1)*beta["k","i"]*V_k*(I_i) + (1-eff1)*beta["k","j"]*V_k*(I_j) + (1-eff1)*beta["k","k"]*V_k*(I_k) + 
      (1-eff1)*(1-eff2)*beta["k","i"]*V_k*(Vi_i) + (1-eff1)*(1-eff2)*beta["k","j"]*V_k*(Vi_j) + (1-eff1)*(1-eff2)*beta["k","k"]*V_k*(Vi_k) - 
      gamma*Vi_k
    
    dVr_k = gamma*Vi_k - sigma2*Vr_k
    
    dC_k = beta["k","i"]*S_k*(I_i) + beta["k","j"]*S_k*(I_j) + beta["k","k"]*S_k*(I_k) + 
      (1-eff2)*beta["k","i"]*S_k*(Vi_i) + (1-eff2)*beta["k","j"]*S_k*(Vi_j) + (1-eff2)*beta["k","k"]*S_k*(Vi_k) + 
      (1-eff1)*beta["k","i"]*V_k*(I_i) + (1-eff1)*beta["k","j"]*V_k*(I_j) + (1-eff1)*beta["k","k"]*V_k*(I_k) + 
      (1-eff1)*(1-eff2)*beta["k","i"]*V_k*(Vi_i) + (1-eff1)*(1-eff2)*beta["k","j"]*V_k*(Vi_j) + (1-eff1)*(1-eff2)*beta["k","k"]*V_k*(Vi_k)
    
    return(list(c(dS_i, dI_i, dR_i, dRv_i, dV_i, dVi_i, dVr_i, dC_i, dCi_i,
                  dS_j, dI_j, dR_j, dRv_j, dV_j, dVi_j, dVr_j, dC_j,
                  dS_k, dI_k, dR_k, dRv_k, dV_k, dVi_k, dVr_k, dC_k)))
  })
} 

# Epidemic Trajectory Plots for the 5 NPI scenarios -----------------------------------

#Initial Conditions and Parameters
init <- c(S_i = (1-(0.0079+0.073))/3, I_i = 0.0079/3, R_i = 0.073/3, Rv_i = 0, V_i = 0, Vi_i = 0, Vr_i = 0, C_i = 0, Ci_i = 0,
          S_j = (1-(0.0079+0.073))/3, I_j = 0.0079/3, R_j = 0.073/3, Rv_j = 0, V_j = 0, Vi_j = 0, Vr_j = 0, C_j = 0,
          S_k = (1-(0.0079+0.073))/3, I_k = 0.0079/3, R_k = 0.073/3, Rv_k = 0, V_k = 0, Vi_k = 0, Vr_k = 0, C_k = 0)

parms = list(gamma = 1/GenTime(3, 2.8),
             sigma1 = 0,
             sigma2 = 0,
             eff1 = 0.9,
             eff2 = 0.9,
             P_i = 0.9,
             P_j = 0.9,
             P_k = 0.9,
             dt_i = 90,
             dt_j = 90,
             dt_k = 90,
             T_i = 0,
             T_j = 90,
             T_k = 180,
             socdist = "Yes",
             vac = "Yes",
             r_i = 0,
             r_j = 0,
             r_k = 0,
             simp = "Yes",
             timing = 270)

#Run the first iteration for (pop i)  - baseline parameters 
  #run the 2nd iteration for (pop j) - using the latter model run
    #run the 3rd iteration for (pop k) - using the latter model run
#Run the finalised model all together with everything that it needs 

#Run to identify Pop i conditions
times <- seq(0,parms[["T_i"]]+1,by = 1)
out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

propSIR_i <- out$S_i[out$time == parms["T_i"]] + out$I_i[out$time == parms["T_i"]] + out$R_i[out$time == parms["T_i"]]
parms["r_i"] <- (propSIR_i*parms[["P_i"]])/parms[["dt_i"]]

#Run to identify Pop j conditions
times <- seq(0,parms[["T_j"]]+1,by = 1)
out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

propSIR_j <- out$S_j[out$time == parms["T_j"]] + out$I_j[out$time == parms["T_j"]] + out$R_j[out$time == parms["T_j"]]
parms["r_j"] <- (propSIR_j*parms[["P_j"]])/parms[["dt_j"]]

#Run to identify Pop k conditions
times <- seq(0,parms[["T_k"]]+1,by = 1)
out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

propSIR_k <- out$S_k[out$time == parms["T_k"]] + out$I_k[out$time == parms["T_k"]] + out$R_k[out$time == parms["T_k"]]
parms["r_k"] <- (propSIR_k*parms[["P_k"]])/parms[["dt_k"]]

#Final Run

times <- seq(0,365,by = 1)
out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out$vac_total_i <- out$Rv_i + out$V_i + out$Vi_i + out$Vr_i
out$vac_total_j <- out$Rv_j + out$V_j + out$Vi_j + out$Vr_j
out$vac_total_k <- out$Rv_k + out$V_k + out$Vi_k + out$Vr_k

data_i <- melt(out, id.vars = c("time"), measure.vars = c("I_i", "vac_total_i", "Vi_i", "R_i"), value.name = "I")
data_j <- melt(out, id.vars = c("time"), measure.vars = c("I_j", "vac_total_j", "Vi_j", "R_j"), value.name = "I")
data_k <- melt(out, id.vars = c("time"), measure.vars = c("I_k", "vac_total_k", "Vi_k", "R_k"), value.name = "I")

# Plotting ----------------------------------------------------------------


max(c(out$I_i, out$I_j, out$I_k))*0.1


#Plot i
shade_i <- data.frame(xmin =  parms[["T_i"]], xmax = parms[["T_i"]]+(parms[["dt_i"]]), ymin = 0, ymax = Inf)

datatexti <- data.frame(x = c(100), y = c(0.45), 
                       label = paste0("Attack Rate ","(i)" , " = ", round(tail(out$C_i, 1), digits =3)))
                                
p_i <- ggplot(data = data_i, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.5), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen", "orange", "darkblue"),
                     labels = c("Nat Inf","Total Vacc","Vacc Inf","Recov")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=18), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.6,0.4,0.4),"cm")) +
  geom_rect(data = shade_i, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Sub-population i") + 
  geom_label(data= datatexti, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, alpha = 0.7, col = "black", fontface = "bold", fill = "white")

p_i2 <- ggplot(data = subset(data_i, variable == "I_i" | variable == "Vi_i"), aes(x = time, y = I, col = variable)) +
  scale_color_manual(values = c("darkred","orange")) + 
  theme_bw() + scale_y_continuous(limits = c(0, max(c(out$I_i, out$I_j, out$I_k))*1.1), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=8), plot.margin=unit(c(0.2,0.4,0,0),"cm"),
        legend.position = "none") +
  geom_rect(data = shade_i, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="", y = "")

p_iplot <- p_i + annotation_custom(ggplotGrob(p_i2), xmin = 200, xmax = 365, 
                       ymin = 0.31, ymax = 0.5)

#Plot j
shade_j <- data.frame(xmin =  parms[["T_j"]], xmax = parms[["T_j"]]+(parms[["dt_j"]]), ymin = 0, ymax = Inf)

datatextj <- data.frame(x = c(100), y = c(0.45), 
                        label = paste0("Attack Rate ","(j)" , " = ", round(tail(out$C_j, 1), digits =3)))

p_j <- ggplot(data = data_j, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.5), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen","orange", "darkblue"),
                     labels = c("Nat Inf","Total Vacc","Vacc Inf","Recov")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=18), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.6,0.4,0.4),"cm")) +
  geom_rect(data = shade_j, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_j"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Sub-population j") + 
  geom_label(data= datatextj, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, alpha = 0.7, col = "black", fontface = "bold", fill = "white")


p_j2 <- ggplot(data = subset(data_j, variable == "I_j" | variable == "Vi_j"), aes(x = time, y = I, col = variable)) +
  scale_color_manual(values = c("darkred","orange")) + 
  theme_bw() + scale_y_continuous(limits = c(0, max(c(out$I_i, out$I_j, out$I_k))*1.1), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=8), plot.margin=unit(c(0.2,0.4,0,0),"cm"),
        legend.position = "none") +
  geom_rect(data = shade_j, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="", y = "")

p_jplot <- p_j + annotation_custom(ggplotGrob(p_j2), xmin = 200, xmax = 365, 
                                   ymin = 0.31, ymax = 0.5)

#Plot k
shade_k <- data.frame(xmin =  parms[["T_k"]], xmax = parms[["T_k"]]+(parms[["dt_k"]]), ymin = 0, ymax = Inf)

datatextk <- data.frame(x = c(100), y = c(0.45), 
                        label = paste0("Attack Rate ","(k)" , " = ", round(tail(out$C_k, 1), digits =3)))

p_k <- ggplot(data = data_k, aes(x = time, y = I, color = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.5), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  scale_color_manual(values = c("darkred","darkgreen","orange", "darkblue"),
                     labels = c("Nat Inf","Total Vacc","Vacc Inf","Recov")) +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=18), legend.title = element_text(size=18),legend.position='bottom',
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.6,0.4,0.4),"cm")) +
  geom_rect(data = shade_k, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_k"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="Time (Days)", y = "Fraction of Population", col = "SubPop", title = "Sub-population k") + 
  geom_label(data= datatextk, inherit.aes = F, aes(x = x, y = y, label = label), size = 5.5, alpha = 0.7, col = "black", fontface = "bold", fill = "white")


p_k2 <- ggplot(data = subset(data_k, variable == "I_k" | variable == "Vi_k"), aes(x = time, y = I, col = variable)) +
  scale_color_manual(values = c("darkred","orange")) + 
  theme_bw() + scale_y_continuous(limits = c(0, max(c(out$I_i, out$I_j, out$I_k))*1.1), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=8), plot.margin=unit(c(0.2,0.4,0,0),"cm"),
        legend.position = "none") +
  geom_rect(data = shade_k, inherit.aes = F, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), alpha = 0.4*as.numeric(parms[["P_i"]]),
            fill = "darkblue") + geom_line(size = 1.1, stat = "identity") + 
  labs(x ="", y = "")

p_kplot <- p_k + annotation_custom(ggplotGrob(p_k2), xmin = 200, xmax = 365, 
                                   ymin = 0.31, ymax = 0.5)

# Attack Rate Analysis ---------------------------------

times <- seq(0,365,by = 1)
parms1 <- parms
parms1[["simp"]] <- "Yes" 

#Baseline
outbase <- out
#After i vacc
parms1[["timing"]] <- 90 
out_i <- data.frame(ode(y = init, func = SIR, times = times, parms = parms1))

#After j vacc
parms1[["timing"]] <- 180 
out_j <- data.frame(ode(y = init, func = SIR, times = times, parms = parms1))

#After k vacc
parms1[["timing"]] <- 270 
out_k <- data.frame(ode(y = init, func = SIR, times = times, parms = parms1))

combframe <- data.frame(cbind("time" = out_i$time, 
                              "outbase" = (outbase$I_i + outbase$I_j + outbase$I_k + outbase$Vi_i + outbase$Vi_j + outbase$Vi_k),
                              "out_i" = (out_i$I_i + out_i$I_j + out_i$I_k + out_i$Vi_i + out_i$Vi_j + out_i$Vi_k), 
                              "out_j" = (out_j$I_i + out_j$I_j + out_j$I_k + out_j$Vi_i + out_j$Vi_j + out_j$Vi_k),
                              "out_k" = (out_k$I_i + out_k$I_j + out_k$I_k + out_k$Vi_i + out_k$Vi_j + out_k$Vi_k)))

data_i_comb <- melt(combframe, id.vars = c("time"), measure.vars = c("outbase", "out_i", "out_j", "out_k"), value.name = "I")

shade_i <- data.frame(xmin =  c(parms[["T_i"]],parms[["T_j"]],parms[["T_k"]]), 
                      xmax = c(parms[["T_i"]]+(parms[["dt_i"]]),
                               parms[["T_j"]]+(parms[["dt_j"]]),
                               parms[["T_k"]]+(parms[["dt_k"]])), ymin = c(0,0,0), ymax = c(Inf,Inf,Inf))

datatext <- data.frame(x = c(182.5, 182.5, 182.5, 182.5), y = c(0.148, 0.140, 0.132, 0.124), 
                       label = c(paste0("Attack Rate (Baseline)", " = ", round(tail(outbase$C_i + outbase$C_j + outbase$C_k, 1), digits =3)), 
                                 paste0("Attack Rate (After i vac, t > 90)", " = ", round(tail(out_i$C_i + out_i$C_j + out_i$C_k, 1), digits =3)),
                                 paste0("Attack Rate (After j vac, t > 180)", " = ", round(tail(out_j$C_i + out_j$C_j + out_j$C_k, 1), digits =3)),
                                 paste0("Attack Rate (After k vac, t > 270)", " = ", round(tail(out_k$C_i + out_k$C_j + out_k$C_k, 1), digits =3))))

p_noint_i <- ggplot(data = data_i_comb, aes(x = time, y = I, col = variable)) + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.15), expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) + 
  scale_color_manual(name = "Scenario",
                        values = c("darkred","darkblue","darkgreen", "purple"), labels = c("Seq Release", "After i" , "After j", "After k")) +
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
                   nrow = 3, ncol = 2, widths = c(1,0.05), common.legend = TRUE, legend = "bottom")

if(parms[["simp"]] == "No") {
  ggsave(pcomb, filename = paste0("vacc_scenario_partial26_seq.png"), dpi = 300, type = "cairo", width = 10, height = 14, units = "in")
} else{
  ggsave(pcomb, filename = paste0("vacc_scenario_partial26_simple_", parms[["timing"]] ,".png"), dpi = 300, type = "cairo", width = 10, height = 14, units = "in")
}

if(parms[["simp"]] == "No") {
  ggsave(p_noint_i, filename = "total_vacc_scenario_partial26_seq.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in")
}


tail(out$Ci_i, 1)
