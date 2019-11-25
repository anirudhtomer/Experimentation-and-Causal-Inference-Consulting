setwd("~/Desktop/StijnGym")

library(ggplot2)
library(foreign)
library(nlme)

load("data/longdata.Rdata")
two_measurements = read.spss("data/results measurements.sav", to.data.frame = T)

final_model = lme(fixed=PO1RM~time + Age + Seks + height + weight + waist,
                  random = ~time|ID, data = longdata, na.action = na.omit)

orig = longdata[!is.na(longdata$PO1RM),]
orig$fitted = fitted(final_model)
orig$resid = orig$PO1RM - orig$fitted

ggplot(orig)  + geom_point(aes(x=time, y=PO1RM)) +
  geom_line(aes(x=time, y=fitted)) + theme_bw() +
  theme(text=element_text(size=15)) +
  facet_wrap(~ID) + xlab("Visit number") + ylab("PO1RM")

#please  change page number to see all patients.
#use R package ggforce
# ggplot(orig)  + geom_point(aes(x=time, y=PO1RM)) +
#   geom_line(aes(x=time, y=fitted)) + theme_bw() +
#   theme(text=element_text(size=15)) +
#   facet_wrap_paginate(~ID, ncol = 3, nrow=3, page = 1) + 
#   xlab("Visit number") + ylab("PO1RM")

# ggplot(orig)  + geom_point(aes(x=time, y=PO1RM)) +
#   geom_line(aes(x=time, y=fitted)) + theme_bw() +
#   theme(text=element_text(size=15)) +
#   facet_wrap_paginate(~ID, ncol = 3, nrow=3, page = 2) + 
#   xlab("Visit number") + ylab("PO1RM")


#Now extracting random slope from the model
slope = rep(NA, 24)
slope[unique(orig$ID)] = ranef(final_model)[,2]

two_measurements$slope = slope
cor.test(two_measurements$GS2-two_measurements$GS1, 
         slope, method = "pearson", conf.level = 0.95)

#removing outliers
outliers = c(14)
cor.test((two_measurements$GS2-two_measurements$GS1)[-outliers], 
         slope[-outliers], method = "pearson", conf.level = 0.95)


cor.test((two_measurements$meanMVC4aeR-two_measurements$MeanMVC1aeR)[-outliers], 
         slope[-outliers], method = "pearson", conf.level = 0.95)

ggplot() + geom_label(aes(x=slope, label = two_measurements$ID,
                          y=two_measurements$GS2-two_measurements$GS1)) +
  theme_bw() + theme(text=element_text(size=15)) +
  xlab("Slope") + ylab("Difference")
