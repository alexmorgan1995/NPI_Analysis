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

    latest<- as.character(Sys.Date()-2) 
  

  
  # Dropbox: gmail account: dailyscraperbox@gmail.com, pw: dsbcid19! (same pw for gmail account and dropbox)
  # The RDS file above should allow you to connect to the dropbox
  # Giles python script must forward those files to that dropbox then
  drop_download(paste0('/daily_scraper_reports/SARS-Cov-2-Scotland-', latest, '_raw.xlsx'), overwrite = TRUE)
  d <- read_excel(paste0('SARS-Cov-2-Scotland-', latest, '_raw.xlsx'))
  d[1, 2] <- "Ayrshire" #Trying to Force the name to be call'able
  
  # PLOTTING EPICURVE FOR SELECTED REGION
  
  
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
      geom_line(size = 1.1)+ xlab('') + ylab('Cumulative number of cases') + 
      scale_y_continuous(expand = c(0,0)) +
      ggtitle(paste0('COVID in ', input$healthboard))+
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
      annotate("rect", xmin = as.character(input$date[1]), 
               xmax = as.character(input$date[2]), ymin = 0, ymax = max(d2$cumNumCases), fill = "lightblue", alpha = .3)
      
    
    
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
      
      min.date<- as.Date(d2$date[which(d2$cumNumCases > 0)[1]])
      max.date<- as.Date(d2$date[nrow(d2)])
      
      t1<- input$date[1]
      t2<- input$date[2]
      
      #t1.check<- as.Date(t1) %within% interval(min.date, max.date)
      #t2.check<- as.Date(t2) %within% interval(min.date, max.date)
      
        
        # COMPUTE Td
        nb.days<- as.numeric(as.Date(t2) - as.Date(t1))
        n.t1<- d2 %>% filter(date == t1) %>% select(cumNumCases)
        n.t2<- d2 %>% filter(date == t2) %>% select(cumNumCases)
        Td<- round(nb.days/(log2(n.t2[[1]]/n.t1[[1]])), 2)
      
      valueBox(paste0("Doubling time: ", Td, " Days"),
               paste0('For ', input$healthboard," between ", input$date[1], " and ", input$date[2]),
               color = "red"
      )
    })
  
  output$instructions <- renderUI({
    HTML(paste("Brief Description of App and Data.",
               "",
               "The start date for doubling time calculations must be made where the number of cases is AT LEAST 1 for the specific healthboard.",
               "",
               "The shaded blue region on the plot describes the range of dates selected for the doubling time calculation.",
               sep ="<br/>"))
  })
  
}
  
