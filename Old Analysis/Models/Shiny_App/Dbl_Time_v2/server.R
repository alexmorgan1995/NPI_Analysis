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
  d[1, 2] <- "Ayrshire" #Trying to Force the name to be call'able
  d[13, 2] <- "Scotland" #Trying to Force the name to be call'able
  # PLOTTING EPICURVE FOR SELECTED REGION
  compute.mu<- function(dat){
    
    # format data
    dat<- dat[dat[,2] > 0,] # trim out first days if necessary (if first new number of case draw is zero) # TO CHECK: WHAT IF THE POISSON DRAW LED THE VERY FIRST non-ZERO DRAW TO BE > 1 ? 
    dat[,1]<- as.Date(dat[,1])
    dat$day_since_start<- c(0, cumsum(as.numeric(diff(dat[,1])))) # Get "day since first case" format of dates
    dat$cumIncidence<- cumsum(dat[,2]) # Get cumulative incidence # TO CHECK -> first index of cumulative incidence should always be 1 if day 0 is day of 
    t = 0
    t_di = t
    Ct = dat$cumIncidence[dat$day_since_start == t]
    while(sum(dat$cumIncidence > 2*Ct, na.rm = TRUE) > 0 ){
      t_di<- c(t_di, dat$day_since_start[which(dat$cumIncidence >= 2*Ct)[1]])
      Ct = dat$cumIncidence[dat$day_since_start == tail(t_di, 1)]
      
    }
    d_j<- diff(t_di)
    mu = length(d_j)/(sum(1/d_j)) # = This is the mean doubling time over the entire cumulative incidence sequence
  
    return(round(mu, 3))
    
  }

  sim.epi<- function(df, its){
    df0 <- df %>% filter(cumNumCases > 0) # trim data to first reported case
    df0$numNewCases<- c(df0$cumNumCases[1], diff(df0$cumNumCases)) # TO CHECK: on day 0, number of new cases should ALWAYS be 1, right?
    df0$date<- as.Date(df0$date)
    df0$day_since_start<- c(0, cumsum(as.numeric(diff(df0$date)))) # Get date in terms of "day first case"
    df0$cumIncidence<- cumsum(df0$numNewCases) # TO CHECK -> first index of cumulative incidence should always be 1 if day 0 is day of first reported case, right?
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
    
    colnames(d2)[ncol(d2)] <- (latest) # latest day does not have a colname
    
    d2<- d2 %>%
      gather('date', 'cumNumCases', 1:ncol(d2)) %>%
      mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
      arrange(date)
    
    d2.2<- d2 %>%
      filter(cumNumCases>0) %>%
      mutate(numNewCases = c(cumNumCases[1],diff(cumNumCases))) %>%
      select(date, numNewCases)
    
    # THIS IS THE DOUBLING TIME DEFINED AS THE HARMONIC MEAN OF THE CUMULATIVE INCIDENCE DOUBLING TIMES OVER THE ENTIRE SEQUENCE
    compute.mu(d2.2)
    
    # 2) CONFIDENCE INTERVAL
    d2.3<- sim.epi(d2, 1000)
    mus<- NULL
    sim.indices<- which(substr(colnames(d2.3), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets
    for(i in 1:length(sim.indices)){
      mus<- c(mus, compute.mu(d2.3[,c(1,sim.indices[i])]))
    }
    
    dtmuci<- round(CI(mus, ci = 0.95), 3)
    
    # .... harmonic mean computed from original data lies outside of CI ...
    # report mean doubling time as mean of this distribution instead?
    
    # 3) PLOT IT
    
    d2.4<- d2.3 %>%
      select(-c(numNewCases, day_since_start, cumIncidence))
    
    d2.4$sim.cumNumCases<- cumsum(d2.4[,which(substr(colnames(d2.4), 1, 1) == 'V')])
    d2.4<- d2.4 %>% select(-which(substr(colnames(d2.4), 1, 1) == 'V'))
    colnames(d2.4)[2]<- d2.4$region[1]
    d2.4<- d2.4 %>% select(-region)
    d2.4<- as.data.frame((as.matrix(d2.4)))

    
    d2.6 <- data.frame(d2.4[,1:2])
    
    ggplot(d2.6, aes(x = as.Date(date), y = as.numeric(as.character(d2.6[,2])))) + geom_line(size = 1.02, col = 'black')+
      xlab('') + ylab('Cumulative Number of Cases')+
      ggtitle(paste0('COVID-19 in ', input$healthboard))+
      theme_bw()+ scale_x_date(date_breaks = "2 days", expand = c(0,0)) + scale_y_continuous(expand = c(0,0))  +
      theme(legend.position="none",
            panel.border= element_blank(),
            axis.text.y = element_text(face="bold", colour="black", size=10),
            axis.text.x = element_text(face="bold", colour="black", size=10, angle = 45, vjust=1, hjust=1),
            axis.title.y = element_text(face="bold", colour="black", size=11),
            axis.title.x = element_text(face="bold", colour="black", size=11),
            axis.line.y = element_line(color="black", size = 0.5),
            axis.line.x = element_line(color="black", size = 0.5),
            plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))
  
  })
 

# 2) COMPUTE DOUBLING TIME  ----
  
  # GET USER INPUT DATES
    
    output$progressBox <- renderValueBox({
      
      d2<- as.data.frame(
        d %>%
          filter(Health_Board == input$healthboard) %>% # Get correct row
          select(-c(1,2, ncol(d)-1, ncol(d))))
      
      colnames(d2)[ncol(d2)] <- (latest) # latest day does not have a colname
      
      d2<- d2 %>%
        gather('date', 'cumNumCases', 1:ncol(d2)) %>%
        mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
        arrange(date)
      
      d2.2<- d2 %>%
        filter(cumNumCases>0) %>%
        mutate(numNewCases = c(cumNumCases[1],diff(cumNumCases))) %>%
        select(date, numNewCases)
      
      # THIS IS THE DOUBLING TIME DEFINED AS THE HARMONIC MEAN OF THE CUMULATIVE INCIDENCE DOUBLING TIMES OVER THE ENTIRE SEQUENCE
      Td <- round(compute.mu(d2.2), 2)
      
      # 2) CONFIDENCE INTERVAL
      d2.3<- sim.epi(d2, 1000)
      mus<- NULL
      sim.indices<- which(substr(colnames(d2.3), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets
      for(i in 1:length(sim.indices)){
        mus<- c(mus, compute.mu(d2.3[,c(1,sim.indices[i])]))
      }
      
      dtmuci<- round(CI(mus, ci = 0.95), 2)
      
      # THIS IS THE CONFIDENCE INTERVAL ON THE HARMONIC MEAN OF THE CUMULATIVE INCIDENCE DOUBLING TIMES OVER THE ENTIRE SEQUENC
      
      valueBox(
        value = tags$p(paste0("Doubling time: ", dtmuci['mean'], " Days", " (95% CI: ", dtmuci['lower'] , ' - ', dtmuci['upper'] , ")"), style = "font-size: 50%;"),
               paste0('For ', input$healthboard," between ", as.character(d2.2$date[1]), " and ", as.character(d2.2$date[nrow(d2.2)])),
               color = "red"
      )
    })
  
  output$progressBox1 <- renderValueBox({
    
    d2<- as.data.frame(
      d %>%
        filter(Health_Board == input$healthboard) %>% # Get correct row
        select(-c(1,2, ncol(d)-1, ncol(d))))
    
    colnames(d2)[ncol(d2)] <- (latest) # latest day does not have a colname
    
    d2<- d2 %>%
      gather('date', 'cumNumCases', 1:ncol(d2)) %>%
      mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
      arrange(date)
    
    d2.2<- d2 %>%
      filter(cumNumCases>0) %>%
      mutate(numNewCases = c(cumNumCases[1],diff(cumNumCases))) %>%
      select(date, numNewCases)
    
    totnumdata <- data.frame(d)
    totnum <- totnumdata[nrow(totnumdata),ncol(totnumdata)-2]
    
    # THIS IS THE CONFIDENCE INTERVAL ON THE HARMONIC MEAN OF THE CUMULATIVE INCIDENCE DOUBLING TIMES OVER THE ENTIRE SEQUENC
    
    valueBox(
      value = tags$p(paste0("Total Cases in Scotland: ", totnum), style = "font-size: 50%;"),
      paste0("Between ", as.character(d2.2$date[1]), " and ", as.character(d2.2$date[nrow(d2.2)])),
      color = "green"
    )
  })
  
  
  output$instructions <- renderUI({
    HTML(paste("<strong>This application calculates the doubling time and the resulting epidemic curve for 14 of Scotland's NHS healthboards.</strong>",
               "",
               "This epidemic doubling time is defined as the time until the number of cases in the population doubles",
               "",
               "Uncertainty for the doubling time estimates for the selected healthboard is provided by 95% confidence intervals.",
               "",
               "Doubling time estimates were calculated as described by <em>Kamalich et al, (2020; Preprint)</em>: 
               https://www.medrxiv.org/content/10.1101/2020.02.05.20020750v4",
               "",
               "Case data is provided by Health Protection Scotland and the Scottish Government.",
               "",
               "This application was developed by members of the COVID-19 modelling team at the University of Edinburgh:",
               "",
               "<em>Giles Calder-Gerver, Camille Simonet, Alex Morgan, Bram van Bunnik and with additional support from Epigroup members</em>",
               sep ="<br/>"))
  })
  
}
  
