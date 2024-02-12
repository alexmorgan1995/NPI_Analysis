library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(rdrop2); library(Rmisc); library(lubridate);library("shiny"); library("rsconnect")


setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/Shiny_App/Test")


drop_auth(rdstoken = 'tokenfile.RDS')

latest<- as.character(Sys.Date()-1) 


drop_download(paste0('/Applications/COVID_Scraper/Daily_scraper_reports/SARS-Cov-2-Scotland-', latest, '_raw.xlsx'), overwrite = TRUE)
d<- read_excel(paste0('SARS-Cov-2-Scotland-', latest, '_raw.xlsx'))


totnum <- data.frame(d)
totnum[nrow(d),ncol(d)-2]


user.input.region<- 'Lothian'


d2<- as.data.frame(
  d %>%
    filter(Health_Board == user.input.region) %>% # Get correct row
    select(-c(1,2, ncol(d)-1, ncol(d))))

colnames(d2)[ncol(d2)]<- (latest) # latest day does not have a colname


d2<- d2 %>%
  gather('date', 'cumNumCases', 1:ncol(d2)) %>%
  mutate(region = user.input.region) %>% # Needed for plotting with ggplot
  arrange(date)



ggplot(d2, aes(x = as.Date(date), y = cumNumCases, group = region))+
  geom_line(size = 1.1)+ xlab('') + ylab('Cumulative number of cases')+
  ggtitle(paste0('COVID in ', user.input.region))+
  theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(face="bold", colour="black", size=10, angle = 45, vjust=1, hjust=1),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))


t1<- '2020-03-12'
t2<- '2020-03-17'

nb.days<- as.numeric(as.Date(t2) - as.Date(t1))
n.t1<- d2 %>% filter(date == t1) %>% select(cumNumCases)
n.t2<- d2 %>% filter(date == t2) %>% select(cumNumCases)
Td<- round(nb.days/(log2(n.t2[[1]]/n.t1[[1]])), 3)


paste0("Doubling time: ", Td, ' for ', user.input.region," between ", t1, " and ", t2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)


# FUNCTION TO COMPUTE HARMONIC MEAN OF DOUBLING TIME ON ANY DATASET, of the form (date, numNewCases)
# It will trim the data itself to cut off initial time points until first observed case
compute.mu<- function(dat){
  
  # format data
  dat<- dat[dat[,2] > 0,] # trim out first days if necessary (if first new number of case draw is zero) # TO CHECK: WHAT IF THE POISSON DRAW LED THE VERY FIRST non-ZERO DRAW TO BE > 1 ? 
  dat[,1]<- as.Date(dat[,1])
  dat$day_since_start<- c(0, cumsum(as.numeric(diff(dat[,1])))) # Get "day since first case" format of dates
  dat$cumIncidence<- cumsum(dat[,2]) # Get cumulative incidence # TO CHECK -> first index of cumulative incidence should always be 1 if day 0 is day of 
  
   
  # Find the sequence of doubling timepoints
  t = 0
  t_di = t
  Ct = dat$cumIncidence[dat$day_since_start == t]
  while(sum(dat$cumIncidence > 2*Ct, na.rm = TRUE) > 0 ){
    
    t_di<- c(t_di, dat$day_since_start[which(dat$cumIncidence >= 2*Ct)[1]])
    Ct = dat$cumIncidence[dat$day_since_start == tail(t_di, 1)]
    
  }
  
  # Convert into sequence of doubling times d_j = delta(t_di)
  d_j<- diff(t_di)
  
  # Harmonic mean of doubling times of cumulative incidence
  mu = length(d_j)/(sum(1/d_j)) # = This is the mean doubling time over the entire cumulative incidence sequence
  
  return(round(mu, 3))
  
}

# FUNCTION TO SIMULATE NEW DATASET OF NUMBER OF NEW CASES
# Takes as argument: its = number of datasets to generate
# df = dataframe of the form (date, numNewCases)

sim.epi<- function(df, its){
  
  
  # FORMAT DATA< GET CUMULATIVE INCIDENCE
  df0 <- df %>% filter(cumNumCases > 0) # trim data to first reported case
  df0$numNewCases<- c(df0$cumNumCases[1], diff(df0$cumNumCases)) # TO CHECK: on day 0, number of new cases should ALWAYS be 1, right?
  df0$date<- as.Date(df0$date)
  df0$day_since_start<- c(0, cumsum(as.numeric(diff(df0$date)))) # Get date in terms of "day first case"
  df0$cumIncidence<- cumsum(df0$numNewCases) # TO CHECK -> first index of cumulative incidence should always be 1 if day 0 is day of first reported case, right?
  
  
  # SIMULATE DATASETS
  # Each timepoint, draw a number of new cases from a poisson distribution of mean the number of new reported cases for that day in the observed data.
  # Directly append the simulated data to dataframe for easier plotting after
  df0<- cbind(df0, as.data.frame(matrix(NA, nrow = nrow(df0), ncol = its)))
  for(j in 1:length(df0$numNewCases)){ # TO CHECK: this to be modified if we decide on computing mu exclding very early phase
    
    df0[j,which(substr(colnames(df0), 1, 1) == "V")]<- rpois(n = its, lambda = df0$numNewCases[j])
  }
  return(df0)
}


# 1) DOUBLING TIME
# Mu on original data --> transform cumNumCases into numNewCases, then feed this into compute.mu()
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

# THIS IS THE CONFIDENCE INTERVAL ON THE HARMONIC MEAN OF THE CUMULATIVE INCIDENCE DOUBLING TIMES OVER THE ENTIRE SEQUENCE
paste0("The doubling time CI_95 is ( ", dtmuci['lower'] , ' - ', dtmuci['upper'] , " )")
hist(mus)

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
#d2.5<- d2.4 %>% gather('set', 'cumNumCases', 2:ncol(d2.4)) # can ignore warning message


#
rearrange <- data.frame(matrix(unlist(d2.4), nrow=length(d2.4), byrow=T), stringsAsFactors=FALSE)
low1000 <- numeric(ncol(rearrange))
max1000 <- numeric(ncol(rearrange))


for(i in 1:ncol(rearrange)) {
  low1000[i] <- min(as.numeric(rearrange[2:nrow(rearrange),i]))
}

for(i in 1:ncol(rearrange)) {
  max1000[i] <- max(as.numeric(rearrange[2:nrow(rearrange),i]))
}

d2.6 <- data.frame(d2.4[,1:2])
d2.6$low <- low1000; d2.6$high <- max1000
d2.7 <- d2.6 %>% gather('set', 'cumNumCases', 2:ncol(d2.6)) # can ignore warning message

ggplot(d2.6, aes(x = as.Date(date), y = as.numeric(as.character(d2.6[,2])))) + geom_line(size = 1.02, col = 'black')+
  xlab('') + ylab('Cumulative number of cases')+
  ggtitle(paste0('COVID in ', user.input.region))+
  geom_ribbon(aes(ymin= low, ymax=high), alpha=0.2) +
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

