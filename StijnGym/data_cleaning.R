setwd("~/Desktop/StijnGym")

longdata = read.csv(file="data/longitudinal.csv", sep = ";", dec = ",", header = T)
longdata$time = 1:48
longdata$height = as.numeric(as.character(longdata$height))

save(longdata, file="data/longdata.Rdata", version = 2)
