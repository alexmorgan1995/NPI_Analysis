rm(list=ls())
library(plotly)

library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(rdrop2); library(Rmisc); library(lubridate);library("shiny"); library("rsconnect"); library(stringr)


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
  rename(NumDeaths = Deaths_New) %>%
  arrange(date)

####################
input <- data.frame("healthboard" = "Scotland")

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
          legend.text = element_text(size=10, face="bold"),
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
