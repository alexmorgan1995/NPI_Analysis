library("deSolve"); library("ggplot2"); library("shiny"); library("shinydashboard"); library("DT")

ui <- dashboardPage(
  dashboardHeader(title = "Effect of Social Distancing Measures - COVID19 Model"),
  dashboardSidebar(
    sliderInput("Time", "Time of Intervention Start",value = 41, min = 0, max = 100),
    sliderInput("LenInt", "Length of Intervention (Weeks)",value = 12, min = 0, max = 20),
    sliderInput("ReducBeta", "Strength of Intervention (Beta Reduction)",value = 0.5, min = 0, max = 1),
    sliderInput("R0", "Baseline Basic Reproduction Number",value = 2, min = 1.1, max = 3, step = 0.1),
    sliderInput("doublingtime", "Epidemic Doubling Time",value = 6, min = 0, max = 10, step = 0.1)
    
    ),
  dashboardBody(
    fluidRow(
      column(width = 7,
             box(
               title = "Epidemic Curve and Transmission Rate", width = NULL, status = "primary",
               plotOutput("Plot1", height = 600), plotOutput("Plot2", height = 200)
             )
      ),
      
      column(width = 3,
             box(title = "Overall % of Population Infected", 
               status = "warning", width = NULL,
               plotOutput("Plot3")
             ),
             box(
               title = "Summary Statistics", width = NULL, solidHeader = TRUE, status = "warning",
               DT::dataTableOutput("table")
             )
      )
    )
  )
)
           
           
    

#mainPanel(plotOutput("Plot1"), plotOutput("Plot2"), plotOutput("Plot3"))  fluidPage(fluidRow(
#column(
#  width = 3,
#  plotOutput('Plot1', height = 200),
#  plotOutput('Plot2', height = 200)
#),
#column(width = 8, plotOutput('Plot3', height = 400)))
#) 
