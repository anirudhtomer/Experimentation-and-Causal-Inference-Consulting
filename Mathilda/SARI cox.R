library(foreign)
library(haven)
library(JM)
library(lattice)

#ALL DATA
data <- read.spss(file.choose(), to.data.frame = T)

#DATA WITH ONLY PEOPLE VALUE FOR CSAMUSCLE
newdata <- data[!is.na(data$CSAmuscle),]

#impute mean height for missings:
#first create new data set with only complete values for height
dataLENGTE <-subset(data,!(is.na(data["Height_cm"])))
#check mean values of height for men and women
tapply(dataLENGTE$Height_cm, dataLENGTE$Gender, mean)

#impute mean height in new data set
newdata$Height_cm <- replace(newdata$Height_cm, (is.na(newdata$Height_cm) & newdata$Gender ==2), 165.2244)
newdata$Height_cm <- replace(newdata$Height_cm, (is.na(newdata$Height_cm) & newdata$Gender ==1), 178.7520)
#create mm2perm2 value
newdata$mm2perm2 <- (newdata$CSAmuscle * 100)/((newdata$Height_cm/100) ^2)

#Make variable with time between operation and follow-up
newdata$FUdays <- (newdata$FUdate - newdata$DateofOR)
newdata$FUdays <- as.numeric(newdata$FUdays)

hist(newdata$mm2perm2)
summary(newdata$mm2perm2)
#nummer 195 HEATH klopt niet? CSA van 11. 
hist(newdata$density)

###COXMODEL BUILDING###
#first check PH assumption: double log of kaplan meiers made above --> see if lines are parallel (slide 229/259)
#Closure method
newdata$HERNIA = ifelse(newdata$HERNIA=="no", yes = 0, 1)

#Don't think Closuremeth is useful
kmfit = survfit(Surv(FUdays, HERNIA)~Closuremeth, conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)

kmfit = survfit(Surv(FUdays, HERNIA)~Smoking, conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)

kmfit = survfit(Surv(FUdays, HERNIA)~Diabetes, conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)


kmfit = survfit(Surv(FUdays, HERNIA)~I(BMI > 25), conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)


kmfit = survfit(Surv(FUdays, HERNIA)~I(Age > 64), conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)

kmfit = survfit(Surv(FUdays, HERNIA)~I(mm2perm2 > 4274), conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)


kmfit = survfit(Surv(FUdays, HERNIA)~I(CSAmuscle > 124.5), conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)

kmfit = survfit(Surv(FUdays, HERNIA)~I(density > 35), conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)


kmfit = survfit(Surv(FUdays, HERNIA)~Cardio, conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)


kmfit = survfit(Surv(FUdays, HERNIA)~I(Woundlength > 22), conf.type="log-log", data=newdata)
survminer::ggsurvplot(kmfit, risk.table = T,
                      break.time.by = 1, xlab = "Time(days)", ylim = c(0.5,1), conf.int = T)

temp = newdata[, c("Woundlength", "Cardio", "COPD", "Diabetes", "mm2perm2", "Closuremeth",
                   "Smoking", "BMI", "Age", "density", "CSAmuscle", "Gender", "FUdays", "HERNIA")]
temp = temp[complete.cases(temp), ]
temp$FUdays = temp$FUdays + 0.00001
temp$HERNIA = ifelse(temp$HERNIA=="no", 0, 1)

cox_All = coxph(Surv(FUdays, HERNIA) ~ Woundlength + Closuremeth + 
                    Cardio + COPD + Diabetes + Smoking + BMI + Age + 
                    density + Gender + mm2perm2,
                data = temp, x=T, model=T)

cv.fit = cv.glmnet(model.matrix(cox_All), Surv(temp$FUdays, temp$HERNIA), 
                   family = "cox", alpha=1, nfolds = 15, lambda = exp(seq(from = 20, to = -20, by = -0.01)),
                   standardize = T)
plot(cv.fit)
log(cv.fit$lambda.min)

fit = glmnet(model.matrix(cox_All), Surv(temp$FUdays, temp$HERNIA), 
             family = "cox", alpha=1, standardize = T)

coef(fit, s = cv.fit$lambda.min)


############################################################
KMClosuremeth <- survfit(Surv(FUdays, HERNIA) ~ Closuremeth, data = newdata)
KMClosuremeth

par(mfrow = c(1,2))

plot(KMClosuremeth, mark.time = FALSE, xlab = "Days", ylab = "Survival", col = 1:2, main = "Closuremethod")
plot(KMClosuremeth, mark.time = FALSE, fun = function (s) -log(-log(s)), xlab = "Days", ylab = "-log(-log(Survival))", col = 1:2)
#PH violated...

#for smoking
KMsmoke <- survfit(Surv(FUdays, HERNIA) ~ Smoking, data=newdata)
KMsmoke

plot(KMsmoke, mark.time = FALSE, xlab = "Days", ylab = "Survival", col = 1:2, main = "Smoking")
plot(KMsmoke, mark.time = FALSE, fun = function (s) -log(-log(s)), xlab = "Days", ylab = "-log(-log(Survival))", col = 1:2)
#PH violated

#for Diabetes
KMdm <- survfit(Surv(FUdays, HERNIA) ~ Diabetes, data=newdata)
KMdm

plot(KMdm, mark.time = FALSE, xlab = "Days", ylab = "Survival", col = 1:2, main = "Diabetes")
plot(KMdm, mark.time = FALSE, fun = function (s) -log(-log(s)), xlab = "Days", ylab = "-log(-log(Survival))", col = 1:2)
#PH violated

#For checking PH in continuous variables: Schoenfeld residuals (with a time-dependant coefficient (sldie 240))
#For mm2perm2
PHmuscle <- coxph(Surv(FUdays, HERNIA) ~ mm2perm2, data = newdata)
zph1 <- cox.zph(PHmuscle, transform = "identity")
zph2 <- cox.zph(PHmuscle, transform = "log")
zph3 <- cox.zph(PHmuscle, transform = "km")

zph1
zph2
zph3

par(mfrow = c(1,1))

plot(zph1)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph2)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph3)
abline(h=0, lty=2, lwd = 2, col = "red")

#For density
PHdensity <- coxph(Surv(FUdays, HERNIA) ~ density, data=newdata)
PHdensity

zph1 <- cox.zph(PHdensity, transform = "identity")
zph2 <- cox.zph(PHdensity, transform = "log")
zph3 <- cox.zph(PHdensity, transform = "km")

zph1
zph2
zph3

plot(zph1)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph2)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph3)
abline(h=0, lty=2, lwd = 2, col = "red")

#for age
PHAge <- coxph(Surv(FUdays, HERNIA) ~ Age, data=newdata)
PHAge

zph1 <- cox.zph(PHAge, transform = "identity")
zph2 <- cox.zph(PHAge, transform = "log")
zph3 <- cox.zph(PHAge, transform = "km")

zph1
zph2
zph3

plot(zph1)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph2)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph3)
abline(h=0, lty=2, lwd = 2, col = "red")

#for BMI
PHbmi <- coxph(Surv(FUdays, HERNIA) ~ BMI, data=newdata)
PHbmi

zph1 <- cox.zph(PHbmi, transform = "identity")
zph2 <- cox.zph(PHbmi, transform = "log")
zph3 <- cox.zph(PHbmi, transform = "km")

zph1
zph2
zph3

plot(zph1)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph2)
abline(h=0, lty=2, lwd = 2, col = "red")
plot(zph3)
abline(h=0, lty=2, lwd = 2, col = "red")

#In most cases PH assumption is violated... Can I use a Coxmodel, or should I add a time-dependent variable?

#BUILD MODEL, if necessary relax assumptions of linearity and additivity

#Check linearity assumption for continuous variables
#linearity assumption mm2perm2
lin.mpm <- coxph(Surv(FUdays, HERNIA) ~ mm2perm2, data = newdata)
nlin.mpm <- coxph(Surv(FUdays, HERNIA) ~ ns(mm2perm2, 3), data = newdata)

ND.mpm <- with(newdata, data.frame(mm2perm2 = seq(min(mm2perm2), max(mm2perm2), length.out = 100)))

prs.l.mpm <- predict(lin.mpm, ND.mpm, se.fit = TRUE, type = "lp")
ND.mpm$pred.l.mpm <- prs.l.mpm[[1]]
ND.mpm$se.l.mpm <- prs.l.mpm[[2]]
ND.mpm$lo.l.mpm <- ND.mpm$pred.l.mpm - 1.96 * ND.mpm$se.l.mpm
ND.mpm$up.l.mpm <- ND.mpm$pred.l.mpm + 1.96 * ND.mpm$se.l.mpm

prs.nl.mpm <- predict(nlin.mpm, ND.mpm, se.fit = TRUE, type = "lp")
ND.mpm$pred.nl.mpm <- prs.nl.mpm[[1]]
ND.mpm$se.nl.mpm <- prs.nl.mpm[[2]]
ND.mpm$lo.nl.mpm <- ND.mpm$pred.nl.mpm - 1.96 * ND.mpm$se.nl.mpm
ND.mpm$up.nl.mpm <- ND.mpm$pred.nl.mpm + 1.96 * ND.mpm$se.nl.mpm

plot(range(ND.mpm$mm2perm2), range(ND.mpm$lo.l.mpm, ND.mpm$up.l.mpm, ND.mpm$lo.nl.mpm, ND.mpm$up.nl.mpm),
     type = "n", xlab = "mm2perm2", ylab = "log Hazard Ratio")
abline(h = 0, lty = 2.5, lwd = 2)
lines(ND.mpm$mm2perm2, ND.mpm$pred.nl.mpm, col = "blue", lwd = 2)
lines(ND.mpm$mm2perm2, ND.mpm$lo.nl.mpm, lty = 2, col = "blue", lwd = 2)
lines(ND.mpm$mm2perm2, ND.mpm$up.nl.mpm, lty = 2, col = "blue", lwd = 2)
lines(ND.mpm$mm2perm2, ND.mpm$pred.l.mpm, col = "red", lwd = 2)
lines(ND.mpm$mm2perm2, ND.mpm$lo.l.mpm, lty = 2, col = "red", lwd = 2)
lines(ND.mpm$mm2perm2, ND.mpm$up.l.mpm, lty = 2, col = "red", lwd = 2)
rug(newdata$mm2perm2)
legend("topright", c("linear", "spline"), lty = 1, lwd = 2,
       col = c("red", "blue"), cex = 1.5, bty = "n")

#Conclusion: mm2perm2 is linearly related to the Log Hazard Ratio

#linearity assumption age
lin.age <- coxph(Surv(FUdays, HERNIA) ~ Age, data = newdata)
nlin.age <- coxph(Surv(FUdays, HERNIA) ~ ns(Age, 3), data = newdata)

ND.age <- with(newdata, data.frame(Age = seq(min(Age), max(Age), length.out = 100)))

prs.l.age <- predict(lin.age, ND.age, se.fit = TRUE, type = "lp")
ND.age$pred.l.age <- prs.l.age[[1]]
ND.age$se.l.age <- prs.l.age[[2]]
ND.age$lo.l.age <- ND.age$pred.l.age - 1.96 * ND.age$se.l.age
ND.age$up.l.age <- ND.age$pred.l.age + 1.96 * ND.age$se.l.age

prs.nl.age <- predict(nlin.age, ND.age, se.fit = TRUE, type = "lp")
ND.age$pred.nl.age <- prs.nl.age[[1]]
ND.age$se.nl.age <- prs.nl.age[[2]]
ND.age$lo.nl.age <- ND.age$pred.nl.age - 1.96 * ND.age$se.nl.age
ND.age$up.nl.age <- ND.age$pred.nl.age + 1.96 * ND.age$se.nl.age

plot(range(ND.age$Age), range(ND.age$lo.l.age, ND.age$up.l.age, ND.age$lo.nl.age, ND.age$up.nl.age),
     type = "n", xlab = "Age", ylab = "log Hazard Ratio")
abline(h = 0, lty = 2.5, lwd = 2)
lines(ND.age$Age, ND.age$pred.nl.age, col = "blue", lwd = 2)
lines(ND.age$Age, ND.age$lo.nl.age, lty = 2, col = "blue", lwd = 2)
lines(ND.age$Age, ND.age$up.nl.age, lty = 2, col = "blue", lwd = 2)
lines(ND.age$Age, ND.age$pred.l.age, col = "red", lwd = 2)
lines(ND.age$Age, ND.age$lo.l.age, lty = 2, col = "red", lwd = 2)
lines(ND.age$Age, ND.age$up.l.age, lty = 2, col = "red", lwd = 2)
rug(newdata$Age)
legend("topright", c("linear", "spline"), lty = 1, lwd = 2,
       col = c("red", "blue"), cex = 1.5, bty = "n")


#linearity assumption bmi
lin.bmi <- coxph(Surv(FUdays, HERNIA) ~ BMI, data = newdata)
nlin.bmi <- coxph(Surv(FUdays, HERNIA) ~ ns(BMI, 3), data = newdata)

ND.bmi <- with(newdata, data.frame(BMI = seq(min(BMI), max(BMI), length.out = 100)))

prs.l.bmi <- predict(lin.bmi, ND.bmi, se.fit = TRUE, type = "lp")
ND.bmi$pred.l.bmi <- prs.l.bmi[[1]]
ND.bmi$se.l.bmi <- prs.l.bmi[[2]]
ND.bmi$lo.l.bmi <- ND.bmi$pred.l.bmi - 1.96 * ND.bmi$se.l.bmi
ND.bmi$up.l.bmi <- ND.bmi$pred.l.bmi + 1.96 * ND.bmi$se.l.bmi

prs.nl.bmi <- predict(nlin.bmi, ND.bmi, se.fit = TRUE, type = "lp")
ND.bmi$pred.nl.bmi <- prs.nl.bmi[[1]]
ND.bmi$se.nl.bmi <- prs.nl.bmi[[2]]
ND.bmi$lo.nl.bmi <- ND.bmi$pred.nl.bmi - 1.96 * ND.bmi$se.nl.bmi
ND.bmi$up.nl.bmi <- ND.bmi$pred.nl.bmi + 1.96 * ND.bmi$se.nl.bmi

plot(range(ND.bmi$BMI), range(ND.bmi$lo.l.bmi, ND.bmi$up.l.bmi, ND.bmi$lo.nl.bmi, ND.bmi$up.nl.bmi),
     type = "n", xlab = "BMI", ylab = "log Hazard Ratio")
abline(h = 0, lty = 2.5, lwd = 2)
lines(ND.bmi$BMI, ND.bmi$pred.nl.bmi, col = "blue", lwd = 2)
lines(ND.bmi$BMI, ND.bmi$lo.nl.bmi, lty = 2, col = "blue", lwd = 2)
lines(ND.bmi$BMI, ND.bmi$up.nl.bmi, lty = 2, col = "blue", lwd = 2)
lines(ND.bmi$BMI, ND.bmi$pred.l.bmi, col = "red", lwd = 2)
lines(ND.bmi$BMI, ND.bmi$lo.l.bmi, lty = 2, col = "red", lwd = 2)
lines(ND.bmi$BMI, ND.bmi$up.l.bmi, lty = 2, col = "red", lwd = 2)
rug(newdata$BMI)
legend("topright", c("linear", "spline"), lty = 1, lwd = 2,
       col = c("red", "blue"), cex = 1.5, bty = "n")

anova(lin.bmi, nlin.bmi)

#linearity assumption density
lin.dens <- coxph(Surv(FUdays, HERNIA) ~ density, data = newdata)
nlin.dens <- coxph(Surv(FUdays, HERNIA) ~ ns(density, 3), data = newdata)

ND.dens <- with(newdata, data.frame(density = seq(min(density), max(density), length.out = 100)))

prs.l.dens <- predict(lin.dens, ND.dens, se.fit = TRUE, type = "lp")
ND.dens$pred.l.dens <- prs.l.dens[[1]]
ND.dens$se.l.dens <- prs.l.dens[[2]]
ND.dens$lo.l.dens <- ND.dens$pred.l.dens - 1.96 * ND.dens$se.l.dens
ND.dens$up.l.dens <- ND.dens$pred.l.dens + 1.96 * ND.dens$se.l.dens

prs.nl.dens <- predict(nlin.dens, ND.dens, se.fit = TRUE, type = "lp")
ND.dens$pred.nl.dens <- prs.nl.dens[[1]]
ND.dens$se.nl.dens <- prs.nl.dens[[2]]
ND.dens$lo.nl.dens <- ND.dens$pred.nl.dens - 1.96 * ND.dens$se.nl.dens
ND.dens$up.nl.dens <- ND.dens$pred.nl.dens + 1.96 * ND.dens$se.nl.dens

plot(range(ND.dens$density), range(ND.dens$lo.l.dens, ND.dens$up.l.dens, ND.dens$lo.nl.dens, ND.dens$up.nl.dens),
     type = "n", xlab = "density", ylab = "log Hazard Ratio")
abline(h = 0, lty = 2.5, lwd = 2)
lines(ND.dens$density, ND.dens$pred.nl.dens, col = "blue", lwd = 2)
lines(ND.dens$density, ND.dens$lo.nl.dens, lty = 2, col = "blue", lwd = 2)
lines(ND.dens$density, ND.dens$up.nl.dens, lty = 2, col = "blue", lwd = 2)
lines(ND.dens$density, ND.dens$pred.l.dens, col = "red", lwd = 2)
lines(ND.dens$density, ND.dens$lo.l.dens, lty = 2, col = "red", lwd = 2)
lines(ND.dens$density, ND.dens$up.l.dens, lty = 2, col = "red", lwd = 2)
rug(newdata$density)
legend("topright", c("linear", "spline"), lty = 1, lwd = 2,
       col = c("red", "blue"), cex = 1.5, bty = "n")

anova(lin.dens, nlin.dens)


#Check additivity: interaction terms. 
fitfull <- coxph(Surv(FUdays, HERNIA) ~ mm2perm2 + Closuremeth + Age + BMI + Smoking + Diabetes + mm2perm2:Closuremeth + mm2perm2:Age + Closuremeth:Smoking + Closuremeth:Diabetes + Age:BMI + Age:Smoking + Age:Diabetes + BMI:Smoking + BMI: Diabetes + Smoking:Diabetes, data=newdata)
#problem with convergence?
fitsimple <- coxph(Surv(FUdays, HERNIA) ~ mm2perm2 + Closuremeth + Age + BMI + Smoking + Diabetes, data=newdata)
anova(fitfull, fitsimple)
#models seem comparable, use fitsimple?

summary(fitsimple)
#mm2perm2 no significant predictor


fitfull2 <- coxph(Surv(FUdays, HERNIA) ~ density + Closuremeth + Age + BMI + Smoking + Diabetes + density:Closuremeth + density:Age + Closuremeth:Smoking + Closuremeth:Diabetes + Age:BMI + Age:Smoking + Age:Diabetes + BMI:Smoking + BMI: Diabetes + Smoking:Diabetes, data = newdata)
summary(fitfull2)
fitsimple2 <- coxph(Surv(FUdays, HERNIA) ~ density + Closuremeth + Age + BMI + Smoking + Diabetes, data = newdata)
anova(fitfull2, fitsimple2)
#Models seem comparable, use fitsimple2?

summary(fitsimple2)
#density no significant predictor

#Useful to make effect plot?
#Transform to survival?
#Check hypothesis of interest: mm2perm2 and density significant predictors?
#PH violated --> does it matter? 
#Competing risks?
