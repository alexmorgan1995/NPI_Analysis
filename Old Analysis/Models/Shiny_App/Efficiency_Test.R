rm(list=ls())

library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(rdrop2); library(Rmisc); library(lubridate); library(plotly)
library(stringr)

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

sim.epi<- function(df, its){
  # FORMAT DATA< GET CUMULATIVE INCIDENCE
  df0 <- df %>% filter(cumNumCases > 0) # trim data to first reported case
  # MODIF TO MAKE:
  df0$numNewCases<- c(df0$cumNumCases[1], diff(df0$cumNumCases)) # to add a cumulative new cases
  #df0$numNewCases<- c(1, diff(df0$cumNumCases)) # TO CHECK
  df0$date<- as.Date(df0$date) # making sure all the character dates are actually dates
  df0$day_since_start<- c(0, cumsum(as.numeric(diff(df0$date)))) # Get date in terms of "day first case"
  df0$cumIncidence<- cumsum(df0$numNewCases) # TO CHECK
  df0<- cbind(df0, as.data.frame(matrix(NA, nrow = nrow(df0), ncol = its))) # initialise the dataframe
  for(j in 1:length(df0$numNewCases)){ # TO CHECK: this to be modified if we decide on computing mu exclding very early phase
    df0[j,which(substr(colnames(df0), 1, 1) == "V")] <- rpois(n = its, lambda = df0$numNewCases[j])
  }
  return(df0)
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
}

input <- data.frame("healthboard" = "Scotland",
                    "date" = as.Date("2020-03-27"))

#### IMPORT MOST RECENT CASES ####
drop_auth(rdstoken = 'tokenfile_NEW.RDS')

files<- drop_dir('/applications/covid_scraper/daily_scraper_reports')

latest <- files[-(grep('death', files$name)), ] %>%
  arrange(client_modified) %>%
  tail(1) %>%
  select(name) %>%
  as.character()

latest.case.date<- str_sub(strsplit(latest, '_')[[1]][1], -10, -1)

drop_download(paste0("/applications/COVID_Scraper/Daily_scraper_reports/", latest), overwrite = TRUE)
d <- read_excel(paste0(latest)); d[13, 2] <- "Scotland"

#### Data Manipulation ####


d2<- as.data.frame(
  d %>%
    filter(Health_Board == input$healthboard) %>% # Get correct row
    select(-c(1,2, ncol(d)-1, ncol(d))))

colnames(d2)[ncol(d2)]<- (latest.case.date) # latest day does not have a colname


d2<- d2 %>%
  gather('date', 'cumNumCases', 1:ncol(d2)) %>%
  mutate(region = input$healthboard) %>% # Needed for plotting with ggplot
  arrange(date)

# APPLY FUNCTION TO COMPUTE TD OVER OBSERVED DATASET
d2.2<- d2 %>%
  mutate(numNewCases = c(cumNumCases[1], diff(cumNumCases))) %>%
  select(date, numNewCases)

td.obs<- compute.mu.m1(d2.2, user.t1 = input$date-7, user.t2 = input$date)
#Doubling Time
paste0("Doubling time: ", td.obs, ' for ', input$healthboard," between ", input$date-7, " and ", input$date)

#Calculating CIs for Doubling Times 
d2.3<- sim.epi(d2, 1000)
Tds<- NULL
sim.indices<- which(substr(colnames(d2.3), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets

for(i in 1:length(sim.indices)){
  Tds <- c(Tds, compute.mu.m1(d2.3[,c(1,sim.indices[i])], user.t1 = input$date-7, user.t2 = input$date))
}

#Create CIs based on 1000 Doubling Time Calculations
alpha = 0.1
ci.low.basic<- round(quantile(Tds, c(0.05), method = 6), 2)
ci.upp.basic<- round(quantile(Tds, c(0.95), method = 6), 2)

paste0("Case Doubling time: ", round(td.obs, digits = 2), " Days (95% CI: ", round(ci.low.basic, digits = 2), ' - ', 
       round(ci.upp.basic, digits = 2) , ")")

#### IMPORT DEATHS ####

latest.death<- files[grep('death', files$name), ] %>% # We will use when plotting death as well
  arrange(client_modified) %>%
  tail(1) %>% 
  select(name) %>%
  as.character()

latest.death.date<- str_sub(strsplit(latest.death, '_')[[1]][1], -10, -1)

drop_download(paste0('/applications/covid_scraper/daily_scraper_reports/SARS-Cov-2-Scotland-', latest.death.date, '_deaths_raw.xlsx'), overwrite = TRUE)

d.death<- as.data.frame(read_excel(paste0('SARS-Cov-2-Scotland-', latest.death.date, '_deaths_raw.xlsx'))[,-1]) %>%
  filter(Date == 'Deaths_New') %>%
  select(-1)

colnames(d.death)[ncol(d.death)]<- (latest.death.date)

#### Data Manipulation - Deaths ####

d.death2<- d.death %>%
  gather('date', 'Deaths_New', 1:ncol(d.death)) %>%
  mutate(cumNumDeath = cumsum(Deaths_New)) %>%
  mutate(region = 'Scotland') %>% # for plotting
  arrange(date)

td.death.obs<- compute.mu.m1(d.death2[,c('date', 'Deaths_New')], user.t1 = input$date-7, user.t2 = input$date)
d.death2.with.sims<- sim.epi.death(d.death2[,c('date', 'cumNumDeath')], 1000)


Tds.deaths<- NULL
sim.indices.deaths<- which(substr(colnames(d.death2.with.sims), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets
for(i in 1:length(sim.indices.deaths)){
  Tds.deaths<- c(Tds.deaths, compute.mu.m1(d.death2.with.sims[,c(1,sim.indices.deaths[i])], user.t1 = input$date-7, user.t2 = input$date))
}

alpha = 0.1
ci.low.basic.death<- round(quantile(Tds.deaths, c(0.05), method = 6), 1)
ci.upp.basic.death<- round(quantile(Tds.deaths, c(0.95), method = 6), 1)

paste0("Death Doubling time: ", round(td.death.obs,2), " Days (95% CI: ", round(ci.low.basic.death, digits = 2), ' - ', round(ci.upp.basic.death, digits = 2) , ")")

