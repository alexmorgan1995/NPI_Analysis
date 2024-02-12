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
    
    dS_i = - beta["i","i"]*S_i*(I_i) - beta["i","j"]*S_i*(I_j) - beta["i","k"]*S_i*(I_k)
    dI_i = beta["i","i"]*S_i*(I_i) + beta["i","j"]*S_i*(I_j) + beta["i","k"]*S_i*(I_k) - gamma*I_i
    dR_i = gamma*I_i
    dC_i = beta["i","i"]*S_i*(I_i) + beta["i","j"]*S_i*(I_j) + beta["i","k"]*S_i*(I_k)
    
    dS_j = - beta["j","i"]*S_j*(I_i) - beta["j","j"]*S_j*(I_j) - beta["j","k"]*S_j*(I_k)
    dI_j = beta["j","i"]*S_j*(I_i) + beta["j","j"]*S_j*(I_j) + beta["j","k"]*S_j*(I_k) - gamma*I_j
    dR_j = gamma*I_j 
    dC_j = beta["j","i"]*S_j*(I_i) + beta["j","j"]*S_j*(I_j) + beta["j","k"]*S_j*(I_k) 
    
    dS_k = - beta["k","i"]*S_k*(I_i) - beta["k","j"]*S_k*(I_j) - beta["k","k"]*S_k*(I_k) 
    dI_k = beta["k","i"]*S_k*(I_i) + beta["k","j"]*S_k*(I_j) + beta["k","k"]*S_k*(I_k) - gamma*I_k
    dR_k = gamma*I_k
    dC_k = beta["k","i"]*S_k*(I_i) + beta["k","j"]*S_k*(I_j) + beta["k","k"]*S_k*(I_k)
    
    return(list(c(dS_i, dI_i, dR_i, dC_i,
                  dS_j, dI_j, dR_j, dC_j,
                  dS_k, dI_k, dR_k, dC_k)))
  })
} 

# Epidemic Trajectory Plots for the 5 NPI scenarios -----------------------------------


#Initial Conditions and Parameters
init <- c(S_i = 0.99999/3, I_i = 0.00001/3, R_i = 0, C_i = 0,
          S_j = 0.99999/3, I_j = 0.00001/3, R_j = 0, C_j = 0,
          S_k = 0.99999/3, I_k = 0.00001/3, R_k = 0, C_k = 0)

times <- seq(0,400,by = 1)
parms = list(gamma = 1/GenTime(3, 2.8),
             beta_matrix = beta_K)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

plot(out$S_i, out$I_i)

plot(out$time, out$I_i)

data <- melt(out, id.vars = c("time"), measure.vars = c("I_i", "I_k"), value.name = "I")

out$C_i[nrow(out)]
out$C_k[nrow(out)]

p1 <- ggplot(data = data, aes(x = time, y = I, color = variable)) + geom_line(size = 1.1, stat = "identity") + 
  theme_bw() + scale_y_continuous(limits = c(0, 0.1),expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) + 
  labs(x = "Time (Days)", y = "Prevalence", col = "SubPop", title = "Beta_K") +
  theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
        legend.text = element_text(size=15), legend.title = element_text(size=18),
        axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) 

# For Loop for All Matrices -----------------------------------------------

listmatrices <- list(beta_base, beta_ALL, beta_K, beta_JK,
                     beta_base_HALF, beta_ALL_HALF, beta_K_HALF, beta_JK_HALF)

storagelist <- list()

for(i in 1:length(listmatrices)) {
  parms = list(gamma = 1/GenTime(3, 2.8),
               beta_matrix = listmatrices[[i]])
  out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
  data <- melt(out, id.vars = c("time"), measure.vars = c("I_i", "I_j", "I_k"), value.name = "I")
  final_i <- round(out$C_i[nrow(out)], digits = 4)
  final_j <- round(out$C_j[nrow(out)], digits = 4)
  final_k <- round(out$C_k[nrow(out)], digits = 4)
  
  datatext <- data.frame(x = c(260, 260, 260), y = c(0.095, 0.085, 0.075), 
                         label = c( paste0("Final~Size~italic(I)[italic(i)]", " ==", final_i), 
                                    paste0("Final~Size~italic(I)[italic(j)]", " ==", final_j),
                                    paste0("Final~Size~italic(I)[italic(k)]", " ==", final_k)))
  
  p1 <- ggplot(data = data, aes(x = time, y = I, color = variable)) + geom_line(size = 1.1, stat = "identity") + 
    theme_bw() + scale_y_continuous(limits = c(0, 0.1),expand = c(0,0)) + scale_x_continuous( expand = c(0, 0)) + 
    labs(x = "Time (Days)", y = "Prevalence", col = "Sub-Population", 
         title = c("beta_base", "beta_ALL", "beta_K", "beta_JK",
                  "beta_base_HALF", "beta_ALL_HALF", "beta_K_HALF", "beta_JK_HALF")[i]) +
    scale_color_manual(values = c("darkred","darkgreen", "darkblue")) +
    theme(axis.title=element_text(size=18), plot.title = element_text(size = 20, vjust = 3, hjust = 0.5, face = "bold"), 
          legend.text = element_text(size=18), legend.title = element_text(size=18),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm")) + 
    geom_label(data= datatext, inherit.aes = F, aes(x = x, y = y, label = label), size = 6, col = "black", 
               parse = TRUE, fontface = "bold", fill = "white")
  
  dump <- list(p1, final_i, final_k)
  storagelist[[i]] <- dump
}


figure <- ggarrange(storagelist[[1]][[1]], storagelist[[3]][[1]], storagelist[[4]][[1]], storagelist[[2]][[1]],
          storagelist[[5]][[1]], storagelist[[7]][[1]], storagelist[[8]][[1]], storagelist[[6]][[1]],
          nrow = 2, ncol = 4, common.legend = TRUE, legend = "bottom")

ggsave(figure, filename = "sub_pops_0.1.png", dpi = 300, type = "cairo", width = 20, height = 10, units = "in")
 