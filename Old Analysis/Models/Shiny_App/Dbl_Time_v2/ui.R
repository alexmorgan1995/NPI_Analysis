library("ggplot2"); library("shiny"); library("shinydashboard"); library("DT")

ui <- dashboardPage(
  dashboardHeader(title = "COVID19 Model"),
  dashboardSidebar(
    selectInput("healthboard", "Scottish Health Board", c("Scotland" = "Scotland",
                                                          "Lothian" = "Lothian",
                                                          "Lanarkshire" = "Lanarkshire",
                                                          "Forth Valley" = "Forth Valley",
                                                          "Fife" = "Fife",
                                                          "Greater Glasgow and Clyde" = "Greater Glasgow and Clyde"
                                                          ))
    ),
  
  
  dashboardBody(
    fluidRow(
      column(width = 6,
             box(
               title = "Epidemic Curve", width = NULL, status = "primary", solidHeader = TRUE,
               plotOutput("Plot1", height = 600)
             )
      ),
      
      column(width = 5,
             valueBoxOutput("progressBox", width = NULL),
             valueBoxOutput("progressBox1", width = NULL),
             box(title = "Instructions", width = NULL, status = "primary", solidHeader = TRUE,
                 htmlOutput( "instructions" ))
             
      )
    )
  )
)
 