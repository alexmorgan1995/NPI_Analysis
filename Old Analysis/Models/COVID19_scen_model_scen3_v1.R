rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/nCoV Work/Models")
library("deSolve"); library("ggplot2"); library("gridExtra"); library("cowplot") 

#### Functions #### 

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

betadecrease <- function(time, int_timestart, int_timeend, betastart, betaend) {
  betalin <- approxfun(x=c(int_timestart, int_timeend),y= c(betastart, betaend), method="linear", rule  =2)
  ifelse((time >= int_timestart & time <= int_timeend),
         betalin(time),
         (2*(1/(GenTime(6,2)))))
}

#plot(betadecrease(seq(100), 32, 90, 2*(1/(GenTime(6,2))), 2*(1/(GenTime(6,2)))*0.25))

SIR3 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - betadecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I
    dI = betadecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I - mu*I
    dR = mu*I 
    dC = betadecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I
    return(list(c(dS,dI,dR, dC)))
  })
}


#### Baseline Model #### 

SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - beta*S*I
    dI = beta*S*I - mu*I
    dR = mu*I 
    dC = beta*S*I
    return(list(c(dS,dI,dR,dC)))
  })
}

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,300,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          beta = (2*(1/(GenTime(6,2)))))
outbase <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
outbase[,6] <- as.numeric(parms[2])


#Plotting
outdata <- rbind(data.frame("Compartment" = "Infec", "Times" = outbase[,1], "Prev" = outbase[,3]),
                  data.frame("Compartment" = "Susc", "Times" = outbase[,1], "Prev" = outbase[,2]),
                  data.frame("Compartment" = "Rec", "Times" = outbase[,1], "Prev" = outbase[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = outbase[,1], "Prev" = outbase[,5]),
                  data.frame("Compartment" = "Beta", "Times" = outbase[,1], "Prev" = outbase[,6]))

ggplot(data = outdata[outdata$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0,0.2) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata[outdata$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") +
  scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

ggplot(data = outdata[outdata$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") +
  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

#### SCENARIO 3 - Linear Decrease #### 

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,300,by = 1)
parms = c(mu = 1/(GenTime(6,2)),
          int_timestart = 35, 
          int_timeend = 35+(12*7),
          beta_base = (2*(1/(GenTime(6,2)))),
          beta_int = (2*(1/(GenTime(6,2))))*0.25)

out <- data.frame(ode(y = init, func = SIR3, times = times, parms = parms))
out$beta <- betadecrease(times, as.numeric(parms[2]), as.numeric(parms[3]), as.numeric(parms[4]), as.numeric(parms[5]))

#Plotting 

outdata3 <- rbind(data.frame("Compartment" = "Infec", "Times" = out[,1], "Prev" = out[,3]),
                  data.frame("Compartment" = "Susc", "Times" = out[,1], "Prev" = out[,2]),
                  data.frame("Compartment" = "Rec", "Times" = out[,1], "Prev" = out[,4]),
                  data.frame("Compartment" = "Cumulative", "Times" = out[,1], "Prev" = out[,5]),
                  data.frame("Compartment" = "Beta", "Times" = out[,1], "Prev" = out[,6]))

p1 <- ggplot(data = outdata3[outdata3$Compartment == "Infec",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

p2 <- ggplot(data = outdata3[outdata3$Compartment == "Beta",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

p3 <- ggplot(data = outdata3[outdata3$Compartment == "Cumulative",], aes(x = Times, y = Prev)) + geom_line(size = 1.05, col = "darkred") +
  labs(x ="Time (Days)", y = "Cumulative Infections") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  scale_x_continuous(expand = c(0, 0)) 

p4 <- ggplot(data = outdata3[1:903,], aes(x = Times, y = Prev, col = Compartment)) + geom_line(size = 1.05) +
  labs(x ="Time (Days)", y = "Proportion of Population") +  scale_y_continuous(limits = c(0,1) ,  expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + 
  scale_x_continuous(expand = c(0, 0)) 

plot_grid(p1, NULL, p2, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))
plot_grid(p3, NULL, p2, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))
plot_grid(p4, NULL, p2, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))

#### Testing Different Interventions  ####

#Create dataframe with the weeks to test 
#Each row has a start and a end date 
#Each row has a the magnistude of the intervention at the start and at the end

weektest <- rbind(data.frame("Week" = "Five", "Start" = 35, "End" = (35+(5*7))),  data.frame("Week" = "Six", "Start" = 35, "End" = (35+(6*7))),
                  data.frame("Week" = "Seven", "Start" = 35, "End" = (35+(7*7))), data.frame("Week" = "Eight", "Start" = 35, "End" = (35+(8*7))),
                  data.frame("Week" = "Nine", "Start" = 35, "End" = (35+(9*7))), data.frame("Week" = "Ten", "Start" = 35, "End" = (35+(10*7))),
                  data.frame("Week" = "Eleven", "Start" = 35, "End" = (35+(11*7))), data.frame("Week" = "Twelve", "Start" = 35, "End" = (35+(12*7))),
                  data.frame("Week" = "Thirteen", "Start" = 35, "End" = (35+(13*7))), data.frame("Week" = "Fourteen", "Start" = 35, "End" = (35+(14*7))),
                  data.frame("Week" = "Fifteen", "Start" = 35, "End" = (35+(15*7))), data.frame("Week" = "Sixteen", "Start" = 35, "End" = (35+(16*7))), 
                  data.frame("Week" = "Seventeen", "Start" = 35, "End" = (35+(17*7))), data.frame("Week" = "Eighteen", "Start" = 35, "End" = (35+(18*7))), 
                  data.frame("Week" = "Nineteen", "Start" = 35, "End" = (35+(19*7))), data.frame("Week" = "Twenty", "Start" = 35, "End" = (35+(20*7))))

init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
times <- seq(0,300,by = 1)

statsrun <- data.frame(matrix(ncol = 4, nrow = nrow(weektest)))
datarun <- data.frame()

for (i in 1:nrow(weektest)) {
  temp <- numeric(4)
  parms = c(mu = 1/(GenTime(6,2)),
            int_timestart = weektest[i,2], 
            int_timeend = weektest[i,3],
            beta_base = (2*(1/(GenTime(6,2)))),
            beta_int = (2*(1/(GenTime(6,2))))*0.25)
  out <- data.frame(ode(y = init, func = SIR3, times = times, parms = parms))
  out$beta <- betadecrease(times, as.numeric(parms[2]), as.numeric(parms[3]), as.numeric(parms[4]), as.numeric(parms[5]))
  out$name <- as.character(weektest[i,1])
  temp[1] <- out[nrow(out), 5]
  temp[2] <- out[,1][which(out[,3] == max(out[,3]))]
  temp[3] <- out[,3][which(out[,3] == max(out[,3]))]
  temp[4] <- as.character(weektest[i,1])
  statsrun[i,] <- temp
  datarun <- rbind.data.frame(datarun, out)
}

colnames(statsrun) <- c("Tot_Inf", "Time_Peak", "Inf_Peak", "Length")
colnames(datarun) <- c("Times", "Susceptible", "Infected", "Recovered", "Cumulative_Inf", "Beta", "Length")

# create graphing function
COVID_length_graph <- function(df){
  length_list <- unique(df$Length)
  for (i in 1:length(length_list)) { 
    p1 <- ggplot(data = subset(df, df$Length==length_list[i]), aes(x = Times, y = Infected)) + geom_line(size = 1.05, col = "darkblue") +
      labs(x ="Time (Days)", y = "Proportion of Infected Humans") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
      scale_x_continuous(expand = c(0, 0))  +  theme(legend.position="none") + 
      ggtitle(paste("Length Of Intervention - ", length_list[i], " Weeks", sep='')) + 
      annotate("rect", xmin = weektest[i,2], xmax = weektest[i,3], ymin = 0, ymax = 0.15, fill = "darkred", alpha = .2) +
      geom_text(x = 250,
                y = 0.1,
                label = statsrun,
                data = )
    p2 <- ggplot(data = subset(df, df$Length==length_list[i]), aes(x = Times, y = Beta)) + geom_line(size = 1.05, col = "darkred") +
      labs(x ="Time (Days)", y = "Beta Parameter") + scale_y_continuous(limits = c(0,0.3) ,  expand = c(0,0)) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            legend.spacing.x = unit(0.3, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
      scale_x_continuous(expand = c(0, 0)) 
    plot <- plot_grid(p1, NULL, p2, align = "v", nrow = 3, rel_heights = c(2, -0.15, 1))
    print(plot)
  }
}

COVID_length_graph(datarun)
