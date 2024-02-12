library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/01_10_20")

# Generic Model Functions ----------------------------------------------------------
#Function for the generation time - a function of R0 and the doubling time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

beta <- function(time, int_timestart, r0) {
  gamma <- (1/GenTime(3, 2.8))
  ifelse((time >= 31 + int_timestart),
         r0*gamma,
         1.4*gamma)
}

beta(seq(0,365), 0 , 0.7)

inc <- data.frame(matrix(nrow = 0, ncol = 2))

#ODEs
SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- beta(time, parms[["int_timestart"]], parms[["r0"]])
    
    inc <<- rbind(inc ,c(time,(beta*S*I)/N))
    inc <<- inc[!duplicated(inc$X0), ]
    
    dS = - (beta*S*I)/N
    dI = (beta*S*I)/N - gamma*I 
    dR = gamma*I 
    
    dC = (beta*S*I)/N
    
    return(list(c(dS, dI, dR, dC)))
  })
} 

# Epidemic Trajectory Plots for the 5 NPI scenarios -----------------------------------

#Initial Conditions and Parameters
init <- c(S = 66650000-8400-1509, I = 8400, R = 0, C = 1509)
times <- seq(0,365,by = 1)
gamma <- (1/GenTime(3, 2.8))
parms = c(gamma = 1/GenTime(3, 2.8),
          beta = 1.4*gamma, #R should be from 1.1 to 2 in steps of 0.1,
          N = 66650000,
          int_timestart = 0,
          r0 = 0.9)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms, method = rk4))
out$beta <- beta(seq(0,365, by = 1), parms[["int_timestart"]], parms[["r0"]])

colnames(inc) <- c("time", "inc")

incalter <- inc[-c(seq(2,times[length(times)]*2, by =2)),]
incalter$date <- seq(as.Date("2020-09-01"), as.Date("2020-09-01") + times[length(times)], by="days")
incalter$beta <- beta(seq(0,365, by = 1), parms[["int_timestart"]], parms[["r0"]])

time1500 <- incalter$time[incalter$inc < 1500][1]

# Model Runs --------------------------------------------------------------

times <- seq(0,time1500,by = 1)

inc <- data.frame(matrix(nrow = 0, ncol = 2))
#Run the model and store/calculate important run characteristics
out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms, method = rk4))
colnames(inc) <- c("time", "inc")
incalter <- inc[-c(seq(2,times[length(times)]*2, by =2)),]
incalter$date <-  seq(as.Date("2020-09-01"), as.Date("2020-09-01") + times[length(times)], by="days")
incalter$inchalf <- incalter$inc/10
incalter$cum <- cumsum(incalter$inchalf)

vec <- c("maxCstart" = incalter$cum[nrow(incalter)], "maxCint" = incalter$cum[nrow(incalter)] -  incalter$cum[incalter$time == 30])
vec[3] <- c("b_a" =  vec[2]/vec[1])

vec[1]
vec[2]
vec[3]


# Plotting ----------------------------------------------------------------

output <- ggplot(incalter, aes(x = date, y = inchalf)) + theme_bw() +
  scale_y_continuous(limits = c(0, 1500), expand = c(0,0)) + scale_x_date(limits = c(as.Date("2020-09-01"), as.Date("2020-09-01") + times[length(times)]), expand = c(0, 0)) + 
  labs(x = "Time", y = "Incidence", col = "", title = "") +
  theme(axis.title.y=element_text(size=18), legend.text = element_text(size=18), legend.title =  element_text(size=18),  axis.title.x=element_text(size=18), plot.title = element_text(size = 22, vjust = 3, hjust = 0.5, face = "bold"), 
        axis.text.x=element_text(size=18), axis.text.y=element_text(size=18), plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm"), legend.position = "bottom") + 
  geom_line(size = 1.1, stat = "identity", col = "darkblue")

ggsave(output, filename = paste0("Scotland_R",parms["r0"],"_intdelay",parms["int_timestart"],".png"), dpi = 300, type = "cairo", width = 7, height = 5, units = "in")

#### Identifying Initial Conditions ####

testing <- expand.grid("r0" = seq(1, 1.5, by = 0.01), "initI" = seq(8400, 8600, by = 1))

data <- data.frame(matrix(nrow = 0, ncol = 4))

for(i in 1:nrow(testing)) {
  print(i/nrow(testing))
  init <- c(S = 66650000-testing$initI[i], I = testing$initI[i], R = 0, C = 0)
  parms = c(gamma = 1/GenTime(3, 2.8),
            beta = testing$r0[i]*gamma, #R should be from 1.1 to 2 in steps of 0.1,
            N = 66650000,
            int_timestart = 0,
            r0 = testing$r0[i])
  inc <- data.frame(matrix(nrow = 0, ncol = 2))
  out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms, method = rk4))
  colnames(inc) <- c("time", "inc")
  data <- rbind(data, c(inc$inc[inc$time == 1], inc$inc[nrow(inc)], testing$r0[i],testing$initI[i]))
}

colnames(data) <- c("inc109", "inc110", "r0", "init")
