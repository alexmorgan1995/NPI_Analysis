
rm(list=ls())
library(plotly)

library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(rdrop2); library(Rmisc); library(lubridate);library("shiny"); library("rsconnect"); library(stringr)


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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SIMULATE DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEATH - GET DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

drop_auth(rdstoken = 'tokenfile_NEW.RDS')

files<- drop_dir('/applications/covid_scraper/daily_scraper_reports')

latest.death<- files[grep('death', files$name), ] %>% # We will use when plotting death as well
  arrange(client_modified) %>%
  tail(1) %>%
  select(name) %>%
  as.character()

latest.death.date<- str_sub(strsplit(latest.death, '_')[[1]][1], -10, -1)

drop_download(paste0('/applications/covid_scraper/daily_scraper_reports/', latest.death), overwrite = TRUE)

d.death<- as.data.frame(read_excel(latest.death)[,-1]) %>%
  filter(Date == 'Deaths_New') %>%
  select(-1)

colnames(d.death)[ncol(d.death)]<- (latest.death.date)

d.death2<- d.death %>%
  gather('date', 'Deaths_New', 1:ncol(d.death)) %>%
  mutate(cumNumDeath = cumsum(Deaths_New)) %>%
  mutate(region = 'Scotland') %>% # for plotting
  rename(NumDeaths = Deaths_New) %>%
  arrange(date)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OBSERVED DOUBLING TIME DEATH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t1<- '2020-03-15'
t2<- '2020-03-25'

td.death.obs<- compute.mu.m1(d.death2[,c('date', 'NumDeaths')], user.t1 = t1, user.t2 = t2)
paste0('DEATH DOUBLING TIME BETWEEN ', t1, ' AND ', t2, ' : ', td.death.obs)



d.death2.with.sims<- sim.epi.death(d.death2[,c('date', 'cumNumDeath')], 1000)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEATH CONFIDENCE INTERVAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Re-compute Td on each simulated dataset
Tds.deaths<- NULL
sim.indices.deaths<- which(substr(colnames(d.death2.with.sims), 1, 1) == 'V') # Get indices of columns corresponding to simulated datasets
for(i in 1:length(sim.indices.deaths)){
  Tds.deaths<- c(Tds.deaths, compute.mu.m1(d.death2.with.sims[,c(1,sim.indices.deaths[i])], user.t1 = t1, user.t2 = t2))
}

# 2) Get CI from the bootstrap distribution
alpha = 0.1
ci.low.basic.death<- round(quantile(Tds.deaths, c(0.05), method = 6), 1)
ci.upp.basic.death<- round(quantile(Tds.deaths, c(0.95), method = 6), 1)

paste0('Between ', t1, ' and ', t2)
paste0("The DEATH doubling time is: ", td.death.obs, " (", ci.low.basic.death, ' - ', ci.upp.basic.death , ")")
