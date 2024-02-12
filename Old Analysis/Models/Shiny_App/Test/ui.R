library("ggplot2"); library("shiny"); library("shinydashboard"); library("DT")

ui <- dashboardPage(
  dashboardHeader(title = "COVID19 Model"),
  dashboardSidebar(
    selectInput("healthboard", "Scottish Health Board", c("Ayrshire and Arran" = "Ayrshire",
                                                          "Dumfries and Galloway" = "Dumfries and Galloway",
                                                          "Fife" = "Fife",
                                                          "Forth Valley" = "Forth Valley",
                                                          "Grampian" = "Grampian",
                                                          "Greater Glasgow and Clyde" = "Greater Glasgow and Clyde",
                                                          "Highland" = "Highland",
                                                          "Lanarkshire" = "Lanarkshire",
                                                          "Lothian" = "Lothian",
                                                          "Shetland" = "Shetland",
                                                          "Tayside" = "Tayside")),
    dateRangeInput("date", "Date range",
                   start  = "2020-03-07",
                   end    = "2020-03-19",
                   min    = "2020-03-01",
                   max    = as.character(latest),
                   format = "yyyy/mm/dd",
                   separator = " - ")
    ),
  
  
  dashboardBody(
    fluidRow(
      column(width = 6,
             box(
               title = "Epidemic Curve", width = NULL, status = "primary", solidHeader = TRUE,
               plotOutput("Plot1", height = 500)
             )
      ),
      
      column(width = 4,
             valueBoxOutput("progressBox", width = NULL),
             box(title = "Instructions", width = NULL, status = "primary", solidHeader = TRUE,
                 htmlOutput( "instructions" ))
             
      )
    )
  )
)
 