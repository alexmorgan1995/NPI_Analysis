rm(list=ls())

library("deSolve"); library("ggplot2"); library("cowplot")

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#### SCENARIO 1 No Intervention ####

SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - beta*S*I
    dI = beta*S*I - gamma*I 
    dR = gamma*I 
    dC = beta*S*I 
    return(list(c(dS,dI,dR,dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001,  R = 0, C = 0)
times <- seq(0,365,by = 1)

parms = c(gamma = 1/(GenTime(4.6,2.4)),
          beta = 1.5*(1/(GenTime(4.6,2.4))))

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

ggplot(out, aes(x = time, y =I)) + geom_line(size = 1.05) + 
  scale_y_continuous(limits = c(0, 0.25), expand = c(0,0)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,365))

ggplot(out, aes(x = time, y =C)) + geom_line(size = 1.05) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,365))

write.csv(comb10data, "delaydata3day.csv", row.names = FALSE)
