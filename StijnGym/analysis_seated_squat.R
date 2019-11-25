setwd("~/Desktop/StijnGym")

library(ggplot2)
library(foreign)
library(nlme)

load("data/longdata.Rdata")
two_measurements = read.spss("data/results measurements.sav", to.data.frame = T)

Ids = c(12,14,21,22,23,24)
diffstart = c(27,14,28,28,23,23)
diffend = c(28, 15,29,29,26,24)
for(i in 1:length(Ids)){
  ID = Ids[i]
  start = diffstart[i]
  end = diffend[i]
  difference =  longdata$SS1RM[longdata$ID==ID & longdata$time==start] - longdata$SS1RM[longdata$ID==ID & longdata$time==end]
  longdata$SS1RM[longdata$ID==ID & longdata$time>=end] = longdata$SS1RM[longdata$ID==ID & longdata$time>=end] + difference
}

final_model = lme(fixed=SS1RM~time + Age + Seks + height + weight + waist,
                  random = ~time|ID, data = longdata, na.action = na.omit)

orig = longdata[!is.na(longdata$SS1RM),]
orig$fitted = fitted(final_model)
orig$resid = orig$SS1RM - orig$fitted

ggplot(orig)  + geom_point(aes(x=time, y=SS1RM)) +
  geom_line(aes(x=time, y=fitted)) + theme_bw() +
  theme(text=element_text(size=15)) +
  facet_wrap(~ID) + xlab("Visit number") + ylab("SS1RM")

#Now extracting random slope from the model
slope = rep(NA, 24)
slope[unique(orig$ID)] = ranef(final_model)[,2]

two_measurements$slope = slope
cor.test(two_measurements$CS30s2-two_measurements$CS30s1, 
         slope, method = "pearson", conf.level = 0.95)
cor.test(two_measurements$CS5T2-two_measurements$CS5T1, 
         slope, method = "pearson", conf.level = 0.95)
cor.test(two_measurements$meanMVC4leL-two_measurements$meanMVC1leL, 
         slope, method = "pearson", conf.level = 0.95)
cor.test(two_measurements$meanMVC4leR-two_measurements$MeanMVC1LeR, 
         slope, method = "pearson", conf.level = 0.95)


ggplot() + geom_label(aes(x=slope, label = two_measurements$ID,
                          y=two_measurements$CS30s2-two_measurements$CS30s1)) +
  theme_bw() + theme(text=element_text(size=15)) +
  xlab("Slope") + ylab("Difference")

ggplot() + geom_label(aes(x=slope, label = two_measurements$ID,
                          y=two_measurements$CS5T2-two_measurements$CS5T1)) +
  theme_bw() + theme(text=element_text(size=15)) +
  xlab("Slope") + ylab("Difference")