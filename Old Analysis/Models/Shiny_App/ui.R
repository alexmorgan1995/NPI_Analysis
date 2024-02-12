library("deSolve"); library("ggplot2"); library("shiny"); library("shinydashboard"); library("DT")

ui <- dashboardPage(
  dashboardHeader(title = "COVID19 Model"),
  dashboardSidebar(
    sliderInput("Time", "Time of Intervention Start",value = 41, min = 0, max = 100),
    sliderInput("LenInt", "Length of Intervention (Weeks)",value = 12, min = 0, max = 20),
    sliderInput("ReducBeta", "Strength of Intervention (Beta Reduction)",value = 0.5, min = 0, max = 1),
    sliderInput("R0", "Baseline Basic Reproduction Number",value = 2, min = 1.1, max = 3, step = 0.1),
    sliderInput("doublingtime", "Epidemic Doubling Time",value = 6, min = 0, max = 10, step = 0.1)
    
    ),
  dashboardBody(
    fluidRow(
      column(width = 8,
             box(
               title = "Epidemic Curve and Transmission Rate", width = NULL, status = "primary", solidHeader = TRUE,
               plotOutput("Plot1", height = 600), plotOutput("Plot2", height = 200)
             )
      ),
      
      column(width = 2,
             valueBoxOutput("progressBox",width = NULL)
             
      )
    )
  )
)
 