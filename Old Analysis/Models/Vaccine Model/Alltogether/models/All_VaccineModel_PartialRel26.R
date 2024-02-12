library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/Vaccine Model/Alltogether")

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

if(parms[["simp"]] == "No") {
  data_i$group <- "partial26_seq"
} else{
  data_i$group <- paste0("partial26_simple_", parms[["timing"]])
  paste0("partial26_simple_", parms[["timing"]])
}

if(parms[["simp"]] == "No") {
  write.csv(data_i,'vacc_scenario_partial26_seq.csv')
} else{
  write.csv(data_i,paste0("partial26_simple_", parms[["timing"]] ,".csv"))
}

tail(out$Ci_i, 1)


