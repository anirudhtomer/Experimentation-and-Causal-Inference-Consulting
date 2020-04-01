if(!require(rms)){
  install.packages("rms")
}

if(!require(foreign)){
  install.packages("foreign")
}

library(rms)
library(foreign)

df = read.spss(file.choose(), to.data.frame = T)
df = df[,2:5]
dd <- datadist(df)
options(datadist="dd")

mod1 = lrm(antibody ~ AID + Behavior + Cognition, 
           data=df, x=TRUE, y=TRUE)
mod1
summary(mod1)
plot(mod1)

plot(anova(mod1), what='proportion chisq') # relative importance
plot(Predict(mod1, fun=plogis)) # predicted values

#Validation results
# IMPORTANT: AUC = Dxy/2 + 0.5
rms::validate(mod1, method="boot", B=500) # bootstrapped validation

#Calibration plot
my.calib <- rms::calibrate(mod1, method="boot", B=500) # model calibration
plot(my.calib, las=1)

