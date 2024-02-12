library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(rdrop2); library(Rmisc); library(lubridate);library("shiny"); library("rsconnect")
library(plotly)


Server <- function(input, output) {


  
  drop_auth(rdstoken = 'tokenfile.RDS')
  
    latest<- as.character(Sys.Date()-1) 

  
  # Dropbox: gmail account: dailyscraperbox@gmail.com, pw: dsbcid19! (same pw for gmail account and dropbox)

    drop_download(paste0('/Applications/COVID_Scraper/Daily_scraper_reports/SARS-Cov-2-Scotland-', latest, '_raw.xlsx'), overwrite = TRUE)
  d <- read_excel(paste0('SARS-Cov-2-Scotland-', latest, '_raw.xlsx'))
  d[1, 2] <- "Ayrshire" 
  d[13, 2] <- "Scotland"

  
  drop_auth(rdstoken = 'tokenfile_NEW.RDS')
  
  files<- drop_dir('/applications/covid_scraper/daily_scraper_reports')
  
  #latest.death<- files[grep('death', files$name), ] %>% # We will use when plotting death as well
   # arrange(client_modified) %>%
    #tail(1) %>%
    #select(name) %>%
    #as.character()
  
  #latest.death.date<- str_sub(strsplit(latest.death, '_')[[1]][1], -10, -1)
  
  #latest.death.date<- str_sub(strsplit(latest, '_')[[1]][1], -10, -1)
  
  drop_download(paste0('/applications/covid_scraper/daily_scraper_reports/SARS-Cov-2-Scotland-', latest, '_deaths_raw.xlsx'), overwrite = TRUE)
  
  d.death<- as.data.frame(read_excel(paste0('SARS-Cov-2-Scotland-', latest, '_deaths_raw.xlsx'))[,-1]) %>%
    filter(Date == 'Deaths_New') %>%
    select(-1)
  
  colnames(d.death)[ncol(d.death)]<- (latest)
  
  d.death2<- d.death %>%
    gather('date', 'Deaths_New', 1:ncol(d.death)) %>%
    mutate(cumNumDeath = cumsum(Deaths_New)) %>%
    mutate(region = 'Scotland') %>% # for plotting
    arrange(date)
  
  
  # FUNCTION COMPUTING Td (method 1)
  # Function takes a dataframe as argument, of the form cbind(date, Number of New Cases)
  compute.mu.m1<- function(dat, user.t1, user.t2){
    user.t1<- as.Date(user.t1)
    user.t2<- as.Date(user.t2)
    dat$cumNewCases<- cumsum(dat[,2])
    dat<- dat %>% filter(cumNewCases > 0)
    nb.days<- as.numeric(as.Date(user.t2) - as.Date(user.t1))
    user.n.t1<- dat %>% filter(date == user.t1) %>% select(cumNewCases)
    user.n.t2<- dat %>% filter(date == user.t2) %>% select(cumNewCases)
    Td<- round(nb.days/(log2(user.n.t2[[1]]/user.n.t1[[1]])), 3)
    return(Td)
  }
  
  sim.epi.death<- function(df, its){
    
    
    # FORMAT DATA< GET CUMULATIVE INCIDENCE
    df0 <- df %>% filter(cumNumDeath > 0) # trim data to first reported case
    
    # MODIF TO MAKE:
    df0$numNewDeaths<- c(df0$cumNumDeath[1], diff(df0$cumNumDeath)) # TO CHECK
    #df0$numNewCases<- c(1, diff(df0$cumNumCases)) # TO CHECK
    df0$date<- as.Date(df0$date)
    df0$day_since_start<- c(0, cumsum(as.numeric(diff(df0$date)))) # Get date in terms of "day first case"
    #df0$cumIncidence<- cumsum(df0$numNewCases) # TO CHECK
    
    
    # SIMULATE DATASETS
    # Each timepoint, draw a number of new cases from a poisson distribution of mean the number of new reported cases for that day in the observed data.
    # Directly append the simulated data to dataframe for easier plotting after
    df0<- cbind(df0, as.data.frame(matrix(NA, nrow = nrow(df0), ncol = its)))
    for(j in 1:length(df0$numNewDeaths)){ 
      
      df0[j,which(substr(colnames(df0), 1, 1) == "V")]<- rpois(n = its, lambda = df0$numNewDeaths[j])
    }
    return(df0)
  } # SIMULATING FUNCTION ADAPTED FOR DEATH (just changing variable names)
  
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
    
    if(input$healthboard == "Scotland") {
      combdeathcase <- rbind(d2,
                             data.frame("date" = d.death2$date,
                                        "cumNumCases" = d.death2$cumNumDeath,
                                        "region" = "deaths"))
      
      combdeathcase$date <- as.Date(combdeathcase$date)
      
      ggplot(combdeathcase, aes(x = as.Date(date), y = cumNumCases, col = region)) +
        geom_line(size = 1.1)+ xlab('') + ylab('Cumulative Number of Cases/Deaths') +
        scale_y_continuous(expand = c(0,0)) +
        ggtitle(paste0('COVID-19 in ', input$healthboard))+
        theme_bw()+
        scale_color_manual(values=c("black", "red"), 
                           breaks=c("Scotland", "deaths"),
                           labels=c("Cases", "Deaths")) +
        theme(legend.position="bottom", legend.title = element_blank(),
              legend.text = element_text(size=12, face="bold"),
              panel.border= element_blank(),
              axis.text.y = element_text(face="bold", colour="black", size=10),
              axis.text.x = element_text(face="bold", colour="black", size=10, angle = 45, vjust=1, hjust=1),
              axis.title.y = element_text(face="bold", colour="black", size=11),
              axis.title.x = element_text(face="bold", colour="black", size=11),
              axis.line.y = element_line(color="black", size = 0.5),
              axis.line.x = element_line(color="black", size = 0.5),
              plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5)) +
        annotate("rect", xmin = as.Date(input$date-7), 
                 xmax = as.Date(input$date), ymin = 0, ymax = max(d2$cumNumCases), fill = "lightblue", alpha = .3) + 
        geom_vline(xintercept = as.Date("2020-03-13"), color = "black", linetype = "dashed", size = 1)
    } else {
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
    }
  })
  
  
  output$Plotcum <- renderPlotly({
    if(input$healthboard == "Scotland") {
    d$pop <- c(
      369670, #Ayrshire 
      371910, #Fife  
      306070, #Forth Valley
      584550, #Grampian 
      1174980, #Greater Glasgow and Clyde
      659200, #Lanarkshire 
      897770, #Lothian
      22990, #Shetland
      416080, #Tayside
      115270, #Borders 
      321800, #Highland
      148790, #Dumfries and Galloway
      5438100) #Scotland
    
    
    dprev <- d
    dprev[,-c(1:2,(ncol(dprev)-2):ncol(dprev))] <- (dprev[,-c(1:2,(ncol(dprev)-2):ncol(dprev))]/d$pop)*10000

    d.keepall<- as.data.frame(dprev) %>%
      select(-c(1, ncol(dprev)-1:2, ncol(dprev)))
    
    colnames(d.keepall)[ncol(d.keepall)]<- as.character(as.Date(colnames(d.keepall)[ncol(d.keepall)-1])+1)
    d.keepall<- d.keepall %>%
      gather('date', 'cumNumCases', 2:ncol(d.keepall)) %>%
      arrange(date)
    d.keepall.dyn <- highlight_key(d.keepall, ~Health_Board )
    

      p <- ggplot(d.keepall.dyn, aes(x = as.Date(date), y = cumNumCases, group = Health_Board))+
        geom_line(size = 0.8)+ xlab('') + ylab('Cases per 10,000 Population')+
        #ggtitle(paste0('COVID in ', user.input.region))+
        theme_bw()+
        theme(#legend.position="none",
          panel.border= element_blank(),
          axis.text.y = element_text(face="bold", colour="black", size=8),
          axis.text.x = element_text(face="bold", colour="black", size=8, angle = 45, vjust=1, hjust=1),
          axis.title.y = element_text(face="bold", colour="black", size=8),
          axis.title.x = element_text(face="bold", colour="black", size=8),
          axis.line.y = element_line(color="black", size = 0.5),
          axis.line.x = element_line(color="black", size = 0.5),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))
      gg <- ggplotly(p, tooltip = "Health_Board" )
      highlight(gg,
                on = "plotly_hover", off = "plotly_deselect",
                color = "black")
    } else {
      
      NULL
    }
    
    
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
      ci.low.basic<- round(quantile(Tds, c(0.05), method = 6), 2)
      #try(ci.low.basic <- round((2*td.obs) - quantile(Tds, probs = c(1 - (alpha/2)))[[1]], 3), silent = TRUE)
      ci.upp.basic <- NA
      ci.upp.basic<- round(quantile(Tds, c(0.95), method = 6), 2)
      #try(ci.upp.basic <- round((2*td.obs) - quantile(Tds, probs = c(alpha/2))[[1]], 3), silent = TRUE)

      valueBox(tags$p(paste0("Case Doubling time: ", round(td.obs, digits = 2), " Days (95% CI: ", round(ci.low.basic, digits = 2), ' - ', round(ci.upp.basic, digits = 2) , ")"), style = "font-size: 50%;"),
               paste0('For ', input$healthboard," between ", input$date-7, " and ", input$date),
               color = "blue"
      )
    })
  
  
  output$progressBoxsum <- renderValueBox({
    
    d$pop <- c(
      369670, #Ayrshire 
      371910, #Fife  
      306070, #Forth Valley
      584550, #Grampian 
      1174980, #Greater Glasgow and Clyde
      659200, #Lanarkshire 
      897770, #Lothian
      22990, #Shetland
      416080, #Tayside
      115270, #Borders 
      321800, #Highland
      148790, #Dumfries and Galloway
      5438100) #Scotland
    
    d$prev1000 <- (d[,ncol(d)-3]/d$pop)*10000
    
    d2<- as.data.frame(
      d %>%
        filter(Health_Board == input$healthboard) %>% # Get correct row
        select(-c(1,2, ncol(d)-1:3, ncol(d))))
    
    colnames(d2)[ncol(d2)]<- (latest) # latest day does not have a colname
    
    d2<- d2 %>%
      gather('date', 'cumNumCases', 1:ncol(d2)) %>%
      mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
      arrange(date)
    
    prev <- as.numeric(as.character((unlist(d[d$Health_Board == input$healthboard,ncol(d)]))))
    
    if(input$healthboard == "Scotland") {
      valueBox(
        value = tags$p(paste0("Total Cases in ", input$healthboard,": ", d2$cumNumCases[nrow(d2)])," Cases (", round(prev, digits = 2), " per 10,000 population)", style = "font-size: 50%;"),
        tags$strong(paste0("Total Deaths in Scotland: ", d.death2$cumNumDeath[nrow(d.death2)]), " Deaths", style = "font-size: 120%;"),
        color = "green"
      )
      
      #valueBox(
      #  value = tags$p(paste0("Total Deaths in Scotland: ", d.death2$cumNumDeath[nrow(d.death2)]," as of ", d.death2$date[nrow(d.death2)]), style = "font-size: 50%;"),
      #  paste0("The Death doubling time is: ", round(td.death.obs,2), " (", ci.low.basic.death, ' - ', ci.upp.basic.death , ")"),
      #  color = "blue"
      #)
      
    } else {
      valueBox(
        value = tags$p(paste0("Total Cases in ", input$healthboard,": ", d2$cumNumCases[nrow(d2)])," Cases as of ", as.character(d2$date[nrow(d2)]), style = "font-size: 50%;"),
        paste0(round(prev, digits = 2), " Cases per 10,000 Population"),
        color = "green"
      )
    }
    
    
    
  })
  
  
  
  output$progressBoxdeath <- renderValueBox({
    
    if(input$healthboard == "Scotland") {
      
      
      td.death.obs<- compute.mu.m1(d.death2[,c('date', 'Deaths_New')], user.t1 = input$date-7, user.t2 = input$date)
      d.death2.with.sims<- sim.epi.death(d.death2[,c('date', 'cumNumDeath')], 1000)
      
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEATH CONFIDENCE INTERVAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 1) Re-compute Td on each simulated dataset
      Tds.deaths<- NULL
      sim.indices.deaths<- which(substr(colnames(d.death2.with.sims), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets
      for(i in 1:length(sim.indices.deaths)){
        Tds.deaths<- c(Tds.deaths, compute.mu.m1(d.death2.with.sims[,c(1,sim.indices.deaths[i])], user.t1 = input$date-7, user.t2 = input$date))
      }
      
      # 2) Get CI from the bootstrap distribution
      alpha = 0.1
      ci.low.basic.death<- round(quantile(Tds.deaths, c(0.05), method = 6), 1)
      ci.upp.basic.death<- round(quantile(Tds.deaths, c(0.95), method = 6), 1)
      
      valueBox(tags$p(paste0("Death Doubling time: ", round(td.death.obs,2), " Days (95% CI: ", round(ci.low.basic.death, digits = 2), ' - ', round(ci.upp.basic.death, digits = 2) , ")"), style = "font-size: 50%;"),
               paste0('For ', input$healthboard," between ", input$date-7, " and ", input$date),
               color = "red"
      )
      
    } else {
      NULL
    }
  })
  
  
  output$instructions <- renderUI({
    HTML(paste0("<strong>This application calculates the doubling time and the resulting epidemic curve for 12 of Scotland's NHS healthboards.</strong>", br(), "Total cases and deaths are recent as of ", d.death2$date[nrow(d.death2)], ".", br(),
               "", br(),
               "This epidemic case/death doubling times are defined as the time until the number of cases/deaths in the population doubles.", br(), 
               "", br(),
               "The time window selection for the the doubling time calculations must be made where the number of cases is <strong>AT LEAST 1</strong> on the start date for the specified healthboard. NA 95% CI output occurs due to negative changes in cumulative case data for certain healthboards. Care must be taken when interpreting doubling times for these healthboards", br(),
               "", br(),
               "The <strong>shaded blue region</strong> on the plot describes the range of dates selected for the doubling time calculation. The <strong>dotted line</strong> denotes the shift in Scottish COVID-19 testing strategy on 2020-03-13 to limit testing to hospitalised cases. Healthboards may be highlighted and explored on the Cases per 10,000 population plot output", br(),
               "", br(),
               "Data is provided by Health Protection Scotland and the Scottish Government.", br(),
               "", br(),
               "This application was developed by members of the COVID-19 modelling team at the University of Edinburgh: <em>Giles Calder-Gerver, Camille Simonet, Alex Morgan, Bram van Bunnik, Mark Woolhouse and with additional support from Epigroup members.</em>", br(),
               "", br(),
               "Contact Email: alex.morgan@ed.ac.uk",
               sep ="<br/>"))
  })
  
}
  
