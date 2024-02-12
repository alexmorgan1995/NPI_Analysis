library("deSolve"); library("ggplot2"); library("shiny"); library("rsconnect")

ui <- pageWithSidebar(
  headerPanel = ("COVID19 Model"),
  sidebarPanel(
    sliderInput("Time", "Time of Intervention Start",value = 41, min = 0, max = 100),
    sliderInput("LenInt", "Length of Intervention (Weeks)",value = 12, min = 0, max = 20),
    sliderInput("R0", "Basic Reproduction Number",value = 2, min = 1.1, max = 3, step = 0.1),
    sliderInput("doublingtime", "Epidemic Doubling Time",value = 6, min = 0, max = 10, step = 0.1),
    sliderInput("ReducBeta", "Beta Reduction",value = 0.5, min = 0, max = 1)
    ),
  fluidPage(
    fluidRow(
      column(12,
             "",
             fluidRow(
               column(10,
                      "Epidemic Curve",
                      plotOutput('Plot1', height = 200),
                      "COVID-19 Transmission Rate",
                      plotOutput('Plot2', height = 100),
                      ),
               column(width = 2,
                      "Population by the End of the Outbreak \n(Cumulative Infections)",
                      plotOutput('Plot3', height = 300),)
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
