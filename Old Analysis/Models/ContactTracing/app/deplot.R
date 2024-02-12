library(rsconnect); library("shiny"); 
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/ContactTracing/app/ContactTracingAppGit")

runApp()

library(rsconnect)
rsconnect::deployApp(appDir = getwd(), appName = "ContactTracing")