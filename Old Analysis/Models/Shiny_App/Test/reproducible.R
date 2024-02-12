
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Models/Shiny_App/Test")

tested <- read_xlsx("SARS-Cov-2-Scotland-2020-03-20_raw.xlsx")
test1 <- read.csv("csvtest.csv")

tested[tested$Health_Board == "Ayrshire and Arran",]
tested[tested$Health_Board == "Lothian",]

test1[test1$Health_Board == "Ayrshire and Arran",]
test1[test1$Health_Board == "Lothian",]
