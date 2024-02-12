library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(rdrop2); library(Rmisc); library(lubridate);library("shiny"); library("rsconnect")


# VARIABLES TO GET FROM USER
# user.input.region
# t1
# t2


# MESSAGES TO DISPLAY
# 1) print('Please enter your region of interest:')
# 2) if free-form input, and that it is wrong: print('Region not in database. Available regions are:'), print(paste(d$Health_Board, collapse = ', '))
# 3) print("Please enter your start date in YYYY-MM-DD format:"), print("Please enter your end date in YYYY-MM-DD format:")
# 4) if free-form input for dates and that it is wrong: print(paste0('Available time window for ', user.input.region, ' are ', min.date, ' to ', max.date))


# CHECKPOINTS
# time of the day (to define latest report to fetch) [present day or day before]
# Is region input by user one of our region (alternatively we could for user to choose in a list instead of free-form input)
# Which dates are selected by user. Cannot compute doubling time over an interval including 0 cases


Server <- function(input, output) {

  
# 1) PLOT CURVE ----  
  # LOG ON DROPBOX AND FETCH LATEST REPORT FILE
  #setwd('Desktop/covid19/COVID_Code/') # The latest report is uploaded from dropbox into this directory, then read.
  
  drop_auth(rdstoken = 'tokenfile.RDS')
  
  # If app run after before 14h, report of 'today' not there yet, so fetches report of today - 1

    latest<- as.character(Sys.Date()-1) 

  
  # Dropbox: gmail account: dailyscraperbox@gmail.com, pw: dsbcid19! (same pw for gmail account and dropbox)
  # The RDS file above should allow you to connect to the dropbox
  # Giles python script must forward those files to that dropbox then
    drop_download(paste0('/Applications/COVID_Scraper/Daily_scraper_reports/SARS-Cov-2-Scotland-', latest, '_raw.xlsx'), overwrite = TRUE)
  d <- read_excel(paste0('SARS-Cov-2-Scotland-', latest, '_raw.xlsx'))
  d[1, 2] <- "Ayrshire" 
  d[13, 2] <- "Scotland"
  #Trying to Force the name to be call'able
  
  # PLOTTING EPICURVE FOR SELECTED REGION
  
  
  
  
  # FUNCTION COMPUTING Td (method 1)
  # Function takes a dataframe as argument, of the form cbind(date, Number of New Cases)
  compute.mu.m1<- function(dat, user.t1, user.t2){
    user.t1<- as.Date(user.t1)
    user.t2<- as.Date(user.t2)
    dat$cumNewCases<- cumsum(dat[,2])
    nb.days<- as.numeric(as.Date(user.t2) - as.Date(user.t1))
    user.n.t1<- dat %>% filter(date == user.t1) %>% select(cumNewCases)
    user.n.t2<- dat %>% filter(date == user.t2) %>% select(cumNewCases)
    Td<- round(nb.days/(log2(user.n.t2[[1]]/user.n.t1[[1]])), 3)
    return(Td)
  }
  # FUNCTION SIMULATING DATASET
  # Takes as argument: its = number of datasets to generate
  # df = dataframe of the form (date, numNewCases)
  sim.epi<- function(df, its){
    # FORMAT DATA< GET CUMULATIVE INCIDENCE
    df0 <- df %>% filter(cumNumCases > 0) # trim data to first reported case
    # MODIF TO MAKE:
    df0$numNewCases<- c(df0$cumNumCases[1], diff(df0$cumNumCases)) # TO CHECK
    #df0$numNewCases<- c(1, diff(df0$cumNumCases)) # TO CHECK
    df0$date<- as.Date(df0$date)
    df0$day_since_start<- c(0, cumsum(as.numeric(diff(df0$date)))) # Get date in terms of "day first case"
    df0$cumIncidence<- cumsum(df0$numNewCases) # TO CHECK
    # SIMULATE DATASETS
    # Each timepoint, draw a number of new cases from a poisson distribution of mean the number of new reported cases for that day in the observed data.
    # Directly append the simulated data to dataframe for easier plotting after
    df0<- cbind(df0, as.data.frame(matrix(NA, nrow = nrow(df0), ncol = its)))
    for(j in 1:length(df0$numNewCases)){ # TO CHECK: this to be modified if we decide on computing mu exclding very early phase
      df0[j,which(substr(colnames(df0), 1, 1) == "V")]<- rpois(n = its, lambda = df0$numNewCases[j])
    }
    return(df0)
  }

  
  
  output$Plot1 <- renderPlot({
    
    
    d2<- as.data.frame(
      d %>%
        filter(Health_Board == input$healthboard) %>% # Get correct row
        select(-c(1,2, ncol(d)-1, ncol(d))))
    
    colnames(d2)[ncol(d2)]<- (latest) # latest day does not have a colname
    
    d2<- d2 %>%
      gather('date', 'cumNumCases', 1:ncol(d2)) %>%
      mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
      arrange(date)
    
    ggplot(d2, aes(x = date, y = cumNumCases, group = region))+
      geom_line(size = 1.1)+ xlab('') + ylab('Cumulative Number of Cases') + 
      scale_y_continuous(expand = c(0,0)) +
      ggtitle(paste0('COVID-19 in ', input$healthboard))+
      theme_bw()+
      theme(legend.position="none",
            panel.border= element_blank(),
            axis.text.y = element_text(face="bold", colour="black", size=10),
            axis.text.x = element_text(face="bold", colour="black", size=10, angle = 45, vjust=1, hjust=1),
            axis.title.y = element_text(face="bold", colour="black", size=11),
            axis.title.x = element_text(face="bold", colour="black", size=11),
            axis.line.y = element_line(color="black", size = 0.5),
            axis.line.x = element_line(color="black", size = 0.5),
            plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5)) +
      annotate("rect", xmin = as.character(input$date-7), 
               xmax = as.character(input$date), ymin = 0, ymax = max(d2$cumNumCases), fill = "lightblue", alpha = .3) + 
      geom_vline(xintercept = "2020-03-13", color = "black", linetype = "dashed", size = 1)
  })
 

# 2) COMPUTE DOUBLING TIME  ----
  
  # GET USER INPUT DATES
    
    output$progressBox <- renderValueBox({
      
      d2<- as.data.frame(
        d %>%
          filter(Health_Board == input$healthboard) %>% # Get correct row
          select(-c(1,2, ncol(d)-1, ncol(d))))
      
      colnames(d2)[ncol(d2)]<- (latest) # latest day does not have a colname

      
      d2<- d2 %>%
        gather('date', 'cumNumCases', 1:ncol(d2)) %>%
        mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
        arrange(date)
      
      
      
      # APPLY FUNCTION TO COMPUTE TD OVER OBSERVED DATASET
      d2.2<- d2 %>%
        mutate(numNewCases = c(cumNumCases[1], diff(cumNumCases))) %>%
        select(date, numNewCases)
      
      td.obs<- compute.mu.m1(d2.2, user.t1 = input$date-7, user.t2 = input$date)
      paste0("Doubling time: ", td.obs, ' for ', input$healthboard," between ", input$date-7, " and ", input$date)
      
      
      # SIMULATE DATA
      d2.3<- sim.epi(d2, 1000)
      Tds<- NULL
      sim.indices<- which(substr(colnames(d2.3), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets
      for(i in 1:length(sim.indices)){
        Tds<- c(Tds, compute.mu.m1(d2.3[,c(1,sim.indices[i])], user.t1 = input$date-7, user.t2 = input$date))
      }
      
      # TEXTBOX OUTPUT
      alpha = 0.1
      ci.low.basic <- NA
      try(ci.low.basic <- round((2*td.obs) - quantile(Tds, probs = c(1 - (alpha/2)))[[1]], 3), silent = TRUE)
      ci.upp.basic <- NA
      try(ci.upp.basic <- round((2*td.obs) - quantile(Tds, probs = c(alpha/2))[[1]], 3), silent = TRUE)
      
      ci.low.basic[ci.low.basic<0] <- 0

      valueBox(tags$p(paste0("Doubling time: ", td.obs, " Days (95% CI: ", ci.low.basic , ' - ', ci.upp.basic , ")"), style = "font-size: 50%;"),
               paste0('For ', input$healthboard," between ", input$date-7, " and ", input$date),
               color = "red"
      )
    })
  
  
  output$progressBox1 <- renderValueBox({
    
    d2<- as.data.frame(
      d %>%
        filter(Health_Board == input$healthboard) %>% # Get correct row
        select(-c(1,2, ncol(d)-1, ncol(d))))
    
    colnames(d2)[ncol(d2)]<- (latest) # latest day does not have a colname
    
    d2<- d2 %>%
      gather('date', 'cumNumCases', 1:ncol(d2)) %>%
      mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
      arrange(date)

    totnumdata <- data.frame(d)
    totnum <- totnumdata[nrow(totnumdata),ncol(totnumdata)-2]
    
    valueBox(
      value = tags$p(paste0("Total Cases in Scotland: ", totnum), style = "font-size: 50%;"),
      paste0("As of ", as.character(d2$date[nrow(d2)])),
      color = "green"
    )
  })
  
  
  output$instructions <- renderUI({
    HTML(paste("<strong>This application calculates the doubling time and the resulting epidemic curve for 12 of Scotland's NHS healthboards.</strong>",
               "",
               "This epidemic doubling time is defined as the time until the number of cases in the population doubles.",
               "",
               "The time window selection for the the doubling time calculations must be made where the number of cases is <strong>AT LEAST 1</strong> on the start date for the specified healthboard. NA 95% CI output occurs due to negative changes in cumulative case data for certain healthboards.",
               "",
               "The <strong>shaded blue region</strong> on the plot describes the range of dates selected for the doubling time calculation. The <strong>dotted line</strong> denotes the shift in Scottish COVID-19 testing strategy on 2020-03-13 to limit testing to hospitalised cases.",
               "",
               "Case data is provided by Health Protection Scotland and the Scottish Government.",
               "",
               "This application was developed by members of the COVID-19 modelling team at the University of Edinburgh: <em>Giles Calder-Gerver, Camille Simonet, Alex Morgan, Bram van Bunnik, Mark Woolhouse and with additional support from Epigroup members.</em>",
               "",
               "Contact Email: alex.morgan@ed.ac.uk",
               sep ="<br/>"))
  })
  
}
  
