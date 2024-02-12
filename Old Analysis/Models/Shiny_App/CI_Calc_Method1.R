rm(list=ls())
library(plotly)

library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(rdrop2); library(Rmisc); library(lubridate);library("shiny"); library("rsconnect")

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



#
dprev <- d
dprev[,-c(1:2,(ncol(dprev)-2):ncol(dprev))] <- (dprev[,-c(1:2,(ncol(dprev)-2):ncol(dprev))]/d$pop)*1000


input <- data.frame("healthboard" = "Ayrshire")
user.input.region <-"Ayrshire"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  PLOT WITH ALL CURVES (--> Dynamic )  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
d.keepall<- as.data.frame(dprev) %>%
  select(-c(1, ncol(dprev)-1:2, ncol(dprev)))

colnames(d.keepall)[ncol(d.keepall)]<- as.character(as.Date(colnames(d.keepall)[ncol(d.keepall)-1])+1)
d.keepall<- d.keepall %>%
  gather('date', 'cumNumCases', 2:ncol(d.keepall)) %>%
  #filter(!Health_Board == 'Total') %>%
  arrange(date)
d.keepall.dyn <- highlight_key(d.keepall, ~Health_Board )


p <- ggplot(d.keepall.dyn, aes(x = as.Date(date), y = cumNumCases, group = Health_Board))+
  geom_line(size = 0.8)+ xlab('') + ylab('Cases per 1000 Population')+
  #ggtitle(paste0('COVID in ', user.input.region))+
  theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(face="bold", colour="black", size=10),
    axis.text.x = element_text(face="bold", colour="black", size=10, angle = 45, vjust=1, hjust=1),
    axis.title.y = element_text(face="bold", colour="black", size=11),
    axis.title.x = element_text(face="bold", colour="black", size=11),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.line.x = element_line(color="black", size = 0.5),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))
gg <- ggplotly(p, tooltip = "Health_Board" )
highlight(gg,
          on = "plotly_hover", off = "plotly_deselect",
          color = "black")

if(input$healthboard == "Scotland") {
  ggplot(d.keepall.dyn, aes(x = as.Date(date), y = cumNumCases, group = Health_Board))+
    geom_line(size = 0.8)+ xlab('') + ylab('Cases per 1000 Population')+
    #ggtitle(paste0('COVID in ', user.input.region))+
    theme_bw()+
    theme(#legend.position="none",
      panel.border= element_blank(),
      axis.text.y = element_text(face="bold", colour="black", size=10),
      axis.text.x = element_text(face="bold", colour="black", size=10, angle = 45, vjust=1, hjust=1),
      axis.title.y = element_text(face="bold", colour="black", size=11),
      axis.title.x = element_text(face="bold", colour="black", size=11),
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














#

d$prev1000 <- (d[,ncol(d)-4]/d$pop)*1000


####
input <- data.frame("healthboard" = "Ayrshire")
user.input.region <-"Ayrshire"


d2<- as.data.frame(
  d %>%
    filter(Health_Board == input$healthboard) %>% # Get correct row
    select(-c(1,2, ncol(d)-1:3, ncol(d))))

colnames(d2)[ncol(d2)]<- (latest) # latest day does not have a colname


d2<- d2 %>%
  gather('date', 'cumNumCases', 1:ncol(d2)) %>%
  mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
  arrange(date)


t1 <- "2020-03-18"

t2 <- "2020-03-25"


# APPLY FUNCTION TO COMPUTE TD OVER OBSERVED DATASET
d2.2<- d2 %>%
  mutate(numNewCases = c(cumNumCases[1], diff(cumNumCases))) %>%
  select(date, numNewCases)

td.obs<- compute.td.m1(d2.2, user.t1 = t1, user.t2 = t2)
paste0("Doubling time: ", td.obs, ' for ', user.input.region," between ", t1, " and ", t2)


# SIMULATE DATA
d2.3<- sim.epi(d2, 1000)
Tds<- NULL
sim.indices<- which(substr(colnames(d2.3), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets
for(i in 1:length(sim.indices)){
  Tds<- c(Tds, compute.mu.m1(d2.3[,c(1,sim.indices[i])], user.t1 = t1, user.t2 = t2))
}


alpha = 0.1
ci.low.basic <- NA
ci.low.basic<- round(quantile(Tds, c(0.05), method = 6), 2)
#try(ci.low.basic <- round((2*td.obs) - quantile(Tds, probs = c(1 - (alpha/2)))[[1]], 3), silent = TRUE)
ci.upp.basic <- NA
ci.upp.basic<- round(quantile(Tds, c(0.95), method = 6), 2)

ci.low.basic[ci.low.basic<0] <- 0
# TEXTBOX OUTPUT
tds.ci<- round(CI(Tds, ci = 0.95), 3)
paste0('For the region: ', user.input.region,", between ", t1, " and ", t2)
paste0("The doubling time is: ", td.obs, " (", ci.low.basic , ' - ', ci.upp.basic  , ")")


