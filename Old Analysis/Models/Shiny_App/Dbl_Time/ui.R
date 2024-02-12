library("ggplot2");library("plotly"); library("shiny"); library("shinydashboard"); library("DT")

ui <- dashboardPage(
  dashboardHeader(title = "COVID19 Model"),
  dashboardSidebar(
    
    selectInput('healthboard', "Area", choices = list(
      National = c("Scotland" = 'Scotland'),
      Healthboard = c("Ayrshire and Arran" = "Ayrshire",
                      "Borders" = "Borders",
                      "Dumfries and Galloway" = "Dumfries and Galloway",
                      "Fife" = "Fife",
                      "Forth Valley" = "Forth Valley",
                      "Grampian" = "Grampian",
                      "Greater Glasgow and Clyde" = "Greater Glasgow and Clyde",
                      "Highland" = "Highland",
                      "Lanarkshire" = "Lanarkshire",
                      "Lothian" = "Lothian",
                      "Shetland" = "Shetland",
                      "Tayside" = "Tayside")), selectize = FALSE),
    
    dateInput("date", "End Date for 7 Day Doubling Time Window",
                   value  = as.character(Sys.Date()-1),
                   min    = "2020-03-08",
                   max    = as.character(Sys.Date()-1),
                   format = "yyyy/mm/dd")
    ),
  
  
  dashboardBody(
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),
    
    
    fluidRow(
      column(width = 6,
             box(
               title = "Epidemic Curve", width = NULL, status = "primary", solidHeader = TRUE,
               plotOutput("Plot1", height = 400),
               plotlyOutput("Plotcum", height = 300)
             )

      ),
      
      column(width = 6,
             valueBoxOutput("progressBoxsum", width = NULL),
             valueBoxOutput("progressBox", width = NULL),
             valueBoxOutput("progressBoxdeath", width = NULL),
             box(title = "Description", width = NULL, status = "primary", solidHeader = TRUE,
                 htmlOutput( "instructions" ))
             
      )
    )
  )
)
 