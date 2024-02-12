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

# +1 Days
out1 <- out; out1[,1] <- out1[,1]+1
outone <- rbind(data.frame("time" = 0, "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out1))
outone <- head(outone, -1)

combdata1 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3])/2),
                        "CumInf" = ((outbase[,5] + outone[,5])/2),
                        "TimeDelay" = "One_Days")

# +2 Days
out2 <- out; out2[,1] <- out2[,1]+2
outtwo <- rbind(data.frame("time" = c(0,1), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                data.frame(out2))
outtwo <- head(outtwo, -2)

combdata2 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3])/3),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5])/3),
                        "TimeDelay" = "Two_Days")

# +3 Days
out3 <- out; out3[,1] <- out3[,1]+3
outthree <- rbind(data.frame("time" = c(seq(0,2)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                data.frame(out3))
outthree <- head(outthree, -3)

combdata3 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3])/4), 
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5])/4),
                        "TimeDelay" = "Three_Days")

# +4 Days
out4 <- out; out4[,1] <- out4[,1]+4
outfour <- rbind(data.frame("time" = c(seq(0,3)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                  data.frame(out3))
outfour <- head(outfour, -4)

combdata4 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3])/5),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5])/5),
                        "TimeDelay" = "Four_Days")

# +5 Days
out5 <- out; out5[,1] <- out5[,1]+5
outfive <- rbind(data.frame("time" = c(seq(0,4)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out5))
outfive <- head(outfive, -5)

combdata5 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3])/6),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5])/6),
                        "TimeDelay" = "Five_Days")

# +6 Days
out6 <- out; out6[,1] <- out6[,1]+6
outsix <- rbind(data.frame("time" = c(seq(0,5)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out6))
outsix <- head(outsix, -6)

combdata6 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3])/7),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5])/7),
                        "TimeDelay" = "Six_Days")

# +7 Days
out7 <- out; out7[,1] <- out7[,1]+7
outseven <- rbind(data.frame("time" = c(seq(0,6)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                data.frame(out7))
outseven <- head(outseven, -7)

combdata7 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3])/8),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5])/8),
                        "TimeDelay" = "Seven_Days")

# +8 Days
out8 <- out; out8[,1] <- out8[,1]+8
outeight <- rbind(data.frame("time" = c(seq(0,7)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                  data.frame(out8))
outeight <- head(outeight, -8)

combdata8 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3] + outeight[,3])/9),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5] + outeight[,5])/9),
                        "TimeDelay" = "Eight_Days")

# +9 Days
out9 <- out; out9[,1] <- out9[,1]+9
outnine <- rbind(data.frame("time" = c(seq(0,8)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                  data.frame(out9))
outnine <- head(outnine, -9)

combdata9 <- data.frame("times" = seq(0,365), "Infec" = ((outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3] + outeight[,3] + outnine[,3])/10),
                        "CumInf" = ((outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5] + outeight[,5] + outnine[,5])/10),
                        "TimeDelay" = "Nine_Days")

# +10 Days
out10 <- out; out10[,1] <- out10[,1]+10
outten <- rbind(data.frame("time" = c(seq(0,9)), "S" = 1, "I" = 0, "R" = 0, "C" = 0),
                 data.frame(out10))
outten <- head(outten, -10)

combdata10 <- data.frame("times" = seq(0,365), "Infec" = (outbase[,3] + outone[,3] + outtwo[,3] + outthree[,3] + outfour[,3] + outfive[,3] + outsix[,3] + 
                                                            outseven[,3] + outeight[,3] + outnine[,3] + outten[,3])/11,
                         "CumInf" = (outbase[,5] + outone[,5] + outtwo[,5] + outthree[,5] + outfour[,5] + outfive[,5] + outsix[,5] + 
                                      outseven[,5] + outeight[,5] + outnine[,5] + outten[,5])/11,
                         "TimeDelay" = "Ten_Days")

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

