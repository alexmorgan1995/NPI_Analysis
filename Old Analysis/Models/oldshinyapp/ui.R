library("deSolve"); library("ggplot2"); library("shiny"); library("rsconnect")

ui <- pageWithSidebar(
  headerPanel = ("COVID19 Model"),
  sidebarPanel(
    sliderInput("Time", "Time of Intervention Start",value = 41, min = 0, max = 100),
    sliderInput("R0", "Basic Reproduction Number",value = 2, min = 1.1, max = 3, step = 0.1),
    sliderInput("doublingtime", "Epidemic Doubling Time",value = 6, min = 0, max = 10, step = 0.1),
    sliderInput("ReducBeta", "Beta Reduction",value = 0.5, min = 0, max = 1)),
  mainPanel(plotOutput("Plot1"), plotOutput("Plot2"))
)