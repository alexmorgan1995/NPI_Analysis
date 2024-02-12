library("ggplot2");library("plotly"); library("shiny"); library("shinydashboard"); library("DT") ; library("shinyjs")

ui <- dashboardPage(
  dashboardHeader(title = "ContactTracing"),
  dashboardSidebar(
    
    useShinyjs(),

    radioButtons("model", "Model Structure:",
                 c("SIR (Lifelong Immunity)" = "SIR",
                   "SIS (No Immunity)" = "SIS"),
                 selected = "SIR"),
    
    
    sliderInput("cont_ramp", "Contact Tracing Ramp-up Rate (per day):",
                min = 1, max = 100,
                value = 50, step = 1),
    
    sliderInput("efficacy", "Efficacy of Contact Tracing):",
                min = 0, max = 1,
                value = 0.5, step = 0.1),
    
    sliderInput("target", "Target Incidence",
                min = 1, max = 100,
                value = 10, step = 1),
    
    sliderInput("R0", "R0",
                min = 0.5, max = 1,
                value = 0.8, step = 0.01),
    
    numericInput("iinput", "I(0) Input", value = 10880, 
                 min = 100, max = 100000, step = 1),
    
    numericInput("rinput", "R(0) Input", value = 460000, 
                 min = 100, max = 1000000, step = 1)
    
    ),
  
  
  dashboardBody(
    
    fluidRow(
      column(width = 6,
             box(
               title = "Epidemic Curve", width = NULL, status = "primary", solidHeader = TRUE,
               plotOutput("Plot1", height = 400),
               plotOutput("Plot2", height = 400)
             )

      ),
      
      column(width = 6,
             valueBoxOutput("progressBoxdeath", width = NULL),
             box(title = "Description", width = NULL, status = "primary", solidHeader = TRUE,
                 htmlOutput( "instructions" ))
             
      )
    )
  )
)
 