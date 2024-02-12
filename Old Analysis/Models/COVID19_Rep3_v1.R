

#### Shiny App - Scen 3 ####
# Define server logic required to plot various variables

Server <- function(input, output) {
  
  betadecrease <- function(time, int_timestart, int_timeend, betastart, betaend) {
    betalin <- approxfun(x=c(int_timestart, int_timeend),y= c(betastart, betaend), method="linear", rule  =2)
    ifelse((time >= int_timestart & time <= int_timeend),
           betalin(time),
           (2*(1/(GenTime(6,2)))))
  }
  
  SIR3 <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS = - betadecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I
      dI = betadecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I - mu*I
      dR = mu*I 
      dC = betadecrease(time,int_timestart,int_timeend, beta_base, beta_int)*S*I
      return(list(c(dS,dI,dR, dC)))
    }
    )
  }
  
  output$Plot1 <- renderPlot({
    times <- seq(0,365, by = 1)
    init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
    parms = c(mu = 1/(GenTime(6,2)), 
              int_timestart = input$Time, 
              int_timeend = input$Time+(12*7), 
              beta_base = 2*(1/(GenTime(6,2))),
              beta_int = (2*(1/(GenTime(6,2)))*(1-0.375))/2)
    out <- data.frame(ode(y = init, func = SIR3, times = times, parms = parms))
    ggplot(out, aes(x = time, y = I)) + geom_line(size = 1.05) +
      labs(x ="Time", y = "Fraction Infected") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
      annotate("rect", xmin = as.numeric(parms[2]), 
               xmax = as.numeric(parms[3]), ymin = 0, ymax = 0.15, fill = "darkred", alpha = .2)
  })
  output$Plot2 <- renderPlot({
    betaplot <- data.frame("Times" = seq(0,365, by = 1), 
                           "Beta" = betadecrease(times, 
                                                 input$Time, 
                                                 input$Time+(12*7),
                                                 2*(1/(GenTime(6,2))),
                                                 (2*(1/(GenTime(6,2)))*(1-0.375))/2))
    ggplot(betaplot, aes(x = Times, y = Beta)) + geom_line(size = 1.05, color = "darkred") +
      labs(x ="Time", y = "Beta") + scale_y_continuous(limits = c(0,0.25) ,  expand = c(0,0))
  })
}

ui <- pageWithSidebar(
  headerPanel = ("COVID19 Model"),
  sidebarPanel(
    sliderInput("Time", "Time of Intervention Start",value = 41, min = 0, max = 100)),
  mainPanel(plotOutput("Plot1"), plotOutput("Plot2"))
)

shinyApp(ui, Server)


#### Shiny App - Scen 1 ####

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
  
  output$Plot1 <- renderPlot({
    times <- seq(0,365, by = 1)
    init <- c(S = 0.9999, I = 0.0001, R = 0, C = 0)
    parms = c(mu = 1/(GenTime(input$doublingtime,input$R0)), 
              int_timestart = input$Time, 
              int_timeend = input$Time+(12*7), 
              beta_base = input$R0*(1/(GenTime(input$doublingtime,input$R0))),
              beta_int = (input$R0*(1/(GenTime(input$doublingtime,input$R0)))*(1-input$ReducBeta))/2)
    out <- data.frame(ode(y = init, func = SIR1, times = times, parms = parms))
    ggplot(out, aes(x = time, y = I)) + geom_line(size = 1.05) +
      labs(x ="Time", y = "Fraction Infected") + scale_y_continuous(limits = c(0,0.15) ,  expand = c(0,0)) +
      annotate("rect", xmin = as.numeric(parms[2]), 
               xmax = as.numeric(parms[3]), ymin = 0, ymax = 0.15, fill = "darkred", alpha = .2)
  })
  output$Plot2 <- renderPlot({
    betaplot <- data.frame("Times" = seq(0,365, by = 1), 
                           "Beta" = betastatdecrease(times, 
                                                     input$Time, 
                                                     input$Time+(12*7),
                                                     input$R0*(1/(GenTime(input$doublingtime,input$R0))),
                                                     input$R0*(1/(GenTime(input$doublingtime,input$R0)))*(1-input$ReducBeta))/2)
    ggplot(betaplot, aes(x = Times, y = Beta)) + geom_line(size = 1.05, color = "darkred") +
      labs(x ="Time", y = "Beta") + scale_y_continuous(limits = c(0,0.25) ,  expand = c(0,0))
  })
}

ui <- pageWithSidebar(
  headerPanel = ("COVID19 Model"),
  sidebarPanel(
    sliderInput("Time", "Time of Intervention Start",value = 41, min = 0, max = 100),
    sliderInput("R0", "Basic Reproduction Number",value = 2, min = 0, max = 3),
    sliderInput("doublingtime", "Epidemic Doubling Time",value = 6, min = 0, max = 10),
    sliderInput("ReducBeta", "Beta Reduction",value = 0.5, min = 0, max = 1)),
  mainPanel(plotOutput("Plot1"), plotOutput("Plot2"))
)

shinyApp(ui, Server)

