library("deSolve"); library("ggplot2"); library("shiny"); library("rsconnect")

Server <- function(input, output) {
  
  GenTime <- function(T2, R0) {
    G = T2 * ((R0-1)/log(2))
    return(G)
  }
  
  betastatdecrease <- function(time, int_timestart, int_timeend, beta_base, beta_int) {
    ifelse((time >= int_timestart & time <= int_timeend),
           beta_int,
           beta_base)
  }
  
  SIR1 <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS = - betastatdecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I
      dI = betastatdecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I - mu*I
      dR = mu*I 
      dC = betastatdecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I
      return(list(c(dS,dI,dR, dC)))
    })
  }
  
  times <- seq(0,365, by = 1)
  init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
  parms = c(mu = 1/(GenTime(input$doublingtime,input$R0)), 
            int_timestart = input$Time, 
            int_timeend = input$Time+(12*7), 
            beta_base = input$R0*(1/(GenTime(input$doublingtime,input$R0))),
            beta_int = (input$R0*(1/(GenTime(input$doublingtime,input$R0)))*(1-input$ReducBeta)))
  out <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
  
  
  stats <- data.frame("stats" = c("TimePeak",  "FracInfPeak", "Infected", "Susceptible"),
                      "value" = as.numeric(c(out[,1][which(out[,3] == max(out[,3]))], out[,3][which(out[,3] == max(out[,3]))],
                                             max(out[,5])*100, (1-max(out[,5]))*100)),
                      "dummy" = c(1))
  
  output$Plot1 <- renderPlot({
    ggplot(out, aes(x = time, y = I)) + geom_line(size = 1.05) +
      labs(x ="Time", y = "Fraction Infected") + scale_y_continuous(limits = c(0,0.25),  expand = c(0,0)) +
      annotate("rect", xmin = as.numeric(parms[2]), 
               xmax = as.numeric(parms[3]), ymin = 0, ymax = 0.25, fill = "darkred", alpha = .3) +
      geom_vline(xintercept= stats$value[stats$stats == "TimePeak"], size = 1.1, color = "darkred") + 
      geom_text(aes(x = stats$value[stats$stats == "TimePeak"] + 40, 
                    label= signif(stats$value[stats$stats == "FracInfPeak"], digits = 2), y= stats$value[stats$stats == "FracInfPeak"] + 0.05),
                colour = "darkred", show.legend=FALSE, size = 5) + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
  })
  
  
  output$Plot2 <- renderPlot({
    betaplot <- data.frame("Times" = seq(0,365, by = 1), 
                           "Beta" = betastatdecrease(seq(0,365, by = 1), 
                                                     input$Time, 
                                                     input$Time+(12*7),
                                                     input$R0*(1/(GenTime(input$doublingtime,input$R0))),
                                                     input$R0*(1/(GenTime(input$doublingtime,input$R0)))*(1-input$ReducBeta)))
    ggplot(betaplot, aes(x = Times, y = Beta)) + geom_line(size = 1.05, color = "darkred") +
      labs(x ="Time", y = "Beta") + scale_y_continuous(limits = c(0,1) ,  expand = c(0,0))
  })
  
  output$Plot3 <- renderPlot({
    ggplot(stats[3:4,], aes(x = dummy, y = value, fill = stats)) + geom_col(position = "stack") +
      geom_text(aes(label = paste0(signif(value, digits = 2), "%")), colour = "white", size = 4,
                position = position_stack(vjust = 0.5)) + 
      theme(legend.position = "right", legend.title = element_blank(), axis.title.y = element_text(margin = margin(r = 20)),
            axis.text.x = element_blank(),  axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      ylab("Percentage") + scale_y_continuous(expand = c(0, 0)) + labs(title = "Population by the End of the Outbreak \n(Cumulative Infections)")
  })
}
