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

#### Dataframe Modification ####

#Baseline 
outbase <- out

ggplot(outbase, aes(x = time, y = I)) + geom_line(size = 1.05, col = "darkblue") +
  labs(x ="Time (Days)", y = "Proportion of Infected Humans") +
  scale_y_continuous(limits = c(0, 0.1), expand = c(0,0)) +
  scale_x_continuous(expand = c(0, 0))

combdatabase <- data.frame("times" = seq(0,365), "Infec" = (outbase[,3]),
                        "CumInf" = (outbase[,5]),
                        "TimeDelay" = "Baseline")

# +7 Days
out1 <- out; out1[,1] <- out1[,1]+7
outone <- rbind(data.frame("time" = seq(0,6), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out1))
outone <- head(outone, -7)

combdata1 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3])/2),
                        "CumInf" = ((outbase[,5] + outone[,5])/2),
                        "TimeDelay" = "6_Days")

# +14 Days
out2 <- out; out2[,1] <- out2[,1]+13
outtwo <- rbind(data.frame("time" = seq(0,12), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                data.frame(out2))
outtwo <- head(outtwo, -13)

combdata2 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3])/3),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5])/3),
                        "TimeDelay" = "12_Days")

# +21 Days
out3 <- out; out3[,1] <- out3[,1]+19
outthree <- rbind(data.frame("time" = c(seq(0,18)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                data.frame(out3))
outthree <- head(outthree, -19)

combdata3 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3])/4), 
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5])/4),
                        "TimeDelay" = "18_Days")

# +28 Days
out4 <- out; out4[,1] <- out4[,1]+25
outfour <- rbind(data.frame("time" = c(seq(0,24)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                  data.frame(out3))
outfour <- head(outfour, -25)

combdata4 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3])/5),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5])/5),
                        "TimeDelay" = "24_Days")

# +35 Days
out5 <- out; out5[,1] <- out5[,1]+31
outfive <- rbind(data.frame("time" = c(seq(0,30)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out5))
outfive <- head(outfive, -31)

combdata5 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3])/6),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5])/6),
                        "TimeDelay" = "30_Days")

# +42 Days
out6 <- out; out6[,1] <- out6[,1]+37
outsix <- rbind(data.frame("time" = c(seq(0,36)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out6))
outsix <- head(outsix, -37)

combdata6 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3])/7),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5])/7),
                        "TimeDelay" = "36_Days")

# +49 Days
out7 <- out; out7[,1] <- out7[,1]+43
outseven <- rbind(data.frame("time" = c(seq(0,42)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                data.frame(out7))
outseven <- head(outseven, -43)

combdata7 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3])/8),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5])/8),
                        "TimeDelay" = "42_Days")

# +56 Days
out8 <- out; out8[,1] <- out8[,1]+49
outeight <- rbind(data.frame("time" = c(seq(0,48)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                  data.frame(out8))
outeight <- head(outeight, -49)

combdata8 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3] + outeight[,3])/9),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5] + outeight[,5])/9),
                        "TimeDelay" = "48_Days")

# +63 Days
out9 <- out; out9[,1] <- out9[,1]+55
outnine <- rbind(data.frame("time" = c(seq(0,54)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                  data.frame(out9))
outnine <- head(outnine, -55)

combdata9 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3] + outeight[,3] + outnine[,3])/10),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5] + outeight[,5] + outnine[,5])/10),
                        "TimeDelay" = "54_Days")

# +10 Days
out10 <- out; out10[,1] <- out10[,1]+61
outten <- rbind(data.frame("time" = c(seq(0,60)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out10))
outten <- head(outten, -61)

combdata10 <- data.frame("times" = seq(0,365), "Infec" = (outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3] + outeight[,3] + outnine[,3] + outten[,3])/11,
                         "CumInf" = (outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5] + outeight[,5] + outnine[,5] + outten[,5])/11,
                         "TimeDelay" = "60_Days")

comb10data <- rbind(combdatabase,
                    combdata1,
                    combdata2,
                    combdata3,
                    combdata4,
                    combdata5,
                    combdata6,
                    combdata7,
                    combdata8,
                    combdata9,
                    combdata10)

ggplot(comb10data, aes(x = times, y =Infec, col = TimeDelay)) + geom_line(size = 1.05)
ggplot(comb10data, aes(x = times, y =CumInf, col = TimeDelay)) + geom_line(size = 1.05)

write.csv(comb10data, "delaydata6day.csv", row.names = FALSE)
