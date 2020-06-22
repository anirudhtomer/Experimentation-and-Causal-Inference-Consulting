#First we load packages and install missing ones
if(!require(foreign)){
  install.packages("foreign")
}
if(!require(nlme)){
  install.packages("nlme")
}
if(!require(ggplot2)){
  install.packages("ggplot2")
}
if(!require(lme4)){
  install.packages("lme4")
}
if(!require(GGally)){
  install.packages("GGally")
}
if(!require(reshape2)){
  install.packages("reshape2")
}
if(!require(compiler)){
  install.packages("compiler")
}
if(!require(parallel)){
  install.packages("parallel")
}
if(!require(boot)){
  install.packages("boot")
}
if(!require(lattice)){
  install.packages("lattice")
}

library(splines)

#Now we load the data. A pop will appear. You must load data in SPSS .sav file format only
patient_data = read.spss(file = file.choose('/Users/jacolinevandendriest/Documents/Huisartsgeneeskunde/Data/syntax/Datasetlongformat_ITT.sav'), to.data.frame = T)

#Data summary
summary(patient_data)

#Treatment has 3 levels in the data, one of them is unnecessary
patient_data$Treatment = droplevels(patient_data$Treatment)
#make usual care as the reference category, so that results are printed for duloxetine
patient_data = within(patient_data, Treatment <- relevel(Treatment, ref = "usual care"))

#Remove all missing WOMAC_pain, missing BMI, missing HADS_depression, missing HADS_anxiety
patient_data = patient_data[!is.na(patient_data$WOMAC_pain) & !is.na(patient_data$BMI) &
                              !is.na(patient_data$HADS_depression) & !is.na(patient_data$HADS_anxiety) &
                              !is.na(patient_data$WOMAC_Stiffness)& !is.na(patient_data$WOMAC_function)&
                              !is.na(patient_data$Com_min2)& !is.na(patient_data$Score_Paindetect) &
                              !is.na(patient_data$Paindetect_categorie)& !is.na(patient_data$FU13herstel_dich)&
                              !is.na(patient_data$FU52herstel_dich)& !is.na(patient_data$FU13_imp)&
                              !is.na(patient_data$FU52_imp)& !is.na(patient_data$FU13AV11_FU13AV11)&
                              !is.na(patient_data$FU26AV11_FU26AV11)& !is.na(patient_data$FU39AV11_FU39AV11)&
                              !is.na(patient_data$EQ_index)& !is.na(patient_data$NVAS)&
                              !is.na(patient_data$FU52AV11_FU52AV11),]

#we make a spline for time effect with knot at median follow up of 26 weeks
print(median(patient_data$Time))

#WOMAC pain.

#The unadjusted model. you can replace ns(Time, knots=26) with just Time to remove splines
unadjusted_model = lme(data=patient_data,
                       fixed = WOMAC_pain~Treatment* ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                       random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                       na.action = na.omit,
                       control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(unadjusted_model)
#95% confidence interval for fixed effectss
intervals(unadjusted_model, which = "fixed")
#p-values
anova(unadjusted_model)
#Complete summary
summary(unadjusted_model)

#The full model now
full_model = lme(data=patient_data,
                 fixed = WOMAC_pain~(age + sex + Treatment + Com_min2 + Score_Paindetect +
                                       HADS_depression) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                 random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                 na.action = na.omit,
                 control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(full_model)
#95% confidence interval for fixed effectss
intervals(full_model, which = "fixed")
#p-values
anova(full_model)
#Complete summary
summary(full_model)

#full model, assumption of homoscedasticity of residuals is satisfied
plot(full_model)

#qqplot of residuals made by two lines is also fine
qqnorm(residuals(full_model))
qqline(residuals(full_model))

#Now check model fit for 9 randomly selected patients. 
#please always run the set.seed(2020) before the ggplot code line
patient_data$modelfit = fitted(full_model)
set.seed(2020)
ggplot(data=patient_data[patient_data$Patientnr %in% sample(levels(patient_data$Patientnr), size = 9, replace = F),]) +
  geom_point(aes(x=Time, y=WOMAC_pain)) +
  geom_line(aes(x=Time, y=modelfit)) + facet_wrap(~factor(Patientnr)) +
  theme_bw() + 
  xlab("Time (weeks)") + ylab("WOMAC Pain")


##########
# Difference in womac pain duloxetine - usual care at 13,26,39,52 weeks
# patient characteristics are same under both duloxetin and usual care 
##########
test_times = c(6, 13, 26, 39, 52)
fixed_effects = fixef(full_model)
for(test_time in test_times){
  contrast = rep(0, length(fixed_effects))
  names(contrast) = names(fixed_effects)
  #number 6 is fixed effect of treatment
  contrast[6] = 1
  #number 19 and 20 are treatment time interactions
  contrast[c(19,20)] = ns(test_time, knots=26, Boundary.knots = c(0,52))
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

##########
# Difference in womac pain duloxetine at 13,26,39,52 weeks and baseline
# The chosen patient's characteristics are set below 
##########
test_times = c(6, 13, 26, 39, 52)
fixed_effects = fixef(full_model)
temp_patient = patient_data[1,]
temp_patient$HADS_depression = median(temp_patient$HADS_depression)
temp_patient$Score_Paindetect = median(temp_patient$Score_Paindetect)
temp_patient$age = median(temp_patient$age)
temp_patient$Com_min2 = factor("nee", levels=c("nee", "ja"))
temp_patient$sex = factor("male", levels = c("female", "male"))
temp_patient$joint = factor("heup", levels = c("heup", "knie"))
temp_patient$Treatment = factor("duloxetine", levels = c("usual care", "duloxetine"))
#confirm patient data before proceeding
View(temp_patient)

for(test_time in test_times){
  temp_patient$Time = test_time
  contrast = c(model.matrix(~(age + sex + Treatment + Com_min2 + Score_Paindetect +
                                HADS_depression) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                            temp_patient))
  contrast[1:8]=0
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

#Final summary plots
#Everything below this point for plotting data
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

#We will make plots for patient with median age, median BMI, median HADS scores, 
#for both joints and both genders, both treatments over time
pred_patient_data = data.frame(expand.grid(Time=seq(0,52,length.out = 10),
                                           Treatment=levels(patient_data$Treatment),
                                           sex=levels(patient_data$sex),
                                           joint=levels(patient_data$joint),
                                           Com_min2=levels(patient_data$Com_min2),
                                           HADS_depression=median(patient_data$HADS_depression),
                                           Score_Paindetect=median(patient_data$Score_Paindetect),
                                           age=median(patient_data$age)))
                                           

plotData = effectPlotData(full_model, newdata = pred_patient_data, patient_data)
#This is the patient data for whom we make a plot
View(plotData)

summary(plotData)

ggplot(data=plotData) + geom_line(aes(x=Time, y=pred, 
                                      color=Treatment, 
                                      group=Treatment)) + 
  geom_ribbon(aes(ymin=low, ymax=upp, x=Time, fill = Treatment), alpha = 0.3)+
  facet_grid(joint~sex)+
  theme_bw() + 
  xlab("Time (weeks)") + ylab("WOMAC Pain")

# WOMAC Stiffness.

#The unadjusted model. you can replace ns(Time, knots=26) with just Time to remove splines
unadjusted_model = lme(data=patient_data,
                       fixed = WOMAC_stiffness~Treatment* ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                       random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                       na.action = na.omit,
                       control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(unadjusted_model)
#95% confidence interval for fixed effectss
intervals(unadjusted_model, which = "fixed")
#p-values
anova(unadjusted_model)
#Complete summary
summary(unadjusted_model)

#The full model now
full_model = lme(data=patient_data,
                 fixed = WOMAC_Stiffness~(BMI + age + sex + joint + Treatment + 
                                       HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                 random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                 na.action = na.omit,
                 control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(full_model)
#95% confidence interval for fixed effectss
intervals(full_model, which = "fixed")
#p-values
anova(full_model)
#Complete summary
summary(full_model)

#full model, assumption of homoscedasticity of residuals is satisfied
plot(full_model)

#qqplot of residuals made by two lines is also fine
qqnorm(residuals(full_model))
qqline(residuals(full_model))

#Now check model fit for 9 randomly selected patients. 
#please always run the set.seed(2020) before the ggplot code line
patient_data$modelfit = fitted(full_model)
set.seed(2023)
ggplot(data=patient_data[patient_data$Patientnr %in% sample(levels(patient_data$Patientnr), size = 9, replace = F),]) +
  geom_point(aes(x=Time, y=WOMAC_Stiffness)) +
  geom_line(aes(x=Time, y=modelfit)) + facet_wrap(~factor(Patientnr)) +
  theme_bw() + 
  xlab("Time (weeks)") + ylab("WOMAC Stiffness")


##########
# Difference in womac pain duloxetine - usual care at 13,26,39,52 weeks
# patient characteristics are same under both duloxetin and usual care 
##########
test_times = c(6, 13, 26, 39, 52)
fixed_effects = fixef(full_model)
for(test_time in test_times){
  contrast = rep(0, length(fixed_effects))
  names(contrast) = names(fixed_effects)
  #number 6 is fixed effect of treatment
  contrast[6] = 1
  #number 19 and 20 are treatment time interactions
  contrast[c(19,20)] = ns(test_time, knots=26, Boundary.knots = c(0,52))
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

##########
# Difference in womac pain duloxetine at 13,26,39,52 weeks and baseline
# The chosen patient's characteristics are set below 
##########
test_times = c(13, 26, 39, 52)
fixed_effects = fixef(full_model)
temp_patient = patient_data[1,]
temp_patient$HADS_depression = median(temp_patient$HADS_depression)
temp_patient$HADS_anxiety = median(temp_patient$HADS_anxiety)
temp_patient$age = median(temp_patient$age)
temp_patient$BMI = median(temp_patient$BMI)
temp_patient$sex = factor("male", levels = c("female", "male"))
temp_patient$joint = factor("heup", levels = c("heup", "knie"))
temp_patient$Treatment = factor("duloxetine", levels = c("usual care", "duloxetine"))
#confirm patient data before proceeding
View(temp_patient)

for(test_time in test_times){
  temp_patient$Time = test_time
  contrast = c(model.matrix(~(BMI + age + sex + joint + Treatment + 
                                HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                            temp_patient))
  contrast[1:8]=0
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

#Final summary plots
#Everything below this point for plotting data
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

#We will make plots for patient with median age, median BMI, median HADS scores, 
#for both joints and both genders, both treatments over time
pred_patient_data = data.frame(expand.grid(Time=seq(0,52,length.out = 10),
                                           Treatment=levels(patient_data$Treatment),
                                           sex=levels(patient_data$sex),
                                           joint=levels(patient_data$joint),
                                           HADS_depression=median(patient_data$HADS_depression),
                                           HADS_anxiety=median(patient_data$HADS_anxiety),
                                           age=median(patient_data$age),
                                           BMI=median(patient_data$BMI)))

plotData = effectPlotData(full_model, newdata = pred_patient_data, patient_data)
#This is the patient data for whom we make a plot
View(plotData)

summary(plotData)

ggplot(data=plotData) + geom_line(aes(x=Time, y=pred, 
                                      color=Treatment, 
                                      group=Treatment)) + 
  geom_ribbon(aes(ymin=low, ymax=upp, x=Time, fill = Treatment), alpha = 0.3)+
  facet_grid(joint~sex)+
  theme_bw() + 
  xlab("Time (weeks)") + ylab("WOMAC Stiffness")

#WOMAC Function.

#The unadjusted model. you can replace ns(Time, knots=26) with just Time to remove splines
unadjusted_model = lme(data=patient_data,
                       fixed = WOMAC_function~Treatment* ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                       random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                       na.action = na.omit,
                       control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(unadjusted_model)
#95% confidence interval for fixed effectss
intervals(unadjusted_model, which = "fixed")
#p-values
anova(unadjusted_model)
#Complete summary
summary(unadjusted_model)

#The full model now
full_model = lme(data=patient_data,
                 fixed = WOMAC_function~(BMI + age + sex + joint + Treatment + 
                                       HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                 random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                 na.action = na.omit,
                 control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(full_model)
#95% confidence interval for fixed effectss
intervals(full_model, which = "fixed")
#p-values
anova(full_model)
#Complete summary
summary(full_model)

#full model, assumption of homoscedasticity of residuals is satisfied
plot(full_model)

#qqplot of residuals made by two lines is also fine
qqnorm(residuals(full_model))
qqline(residuals(full_model))

#Now check model fit for 9 randomly selected patients. 
#please always run the set.seed(2020) before the ggplot code line
patient_data$modelfit = fitted(full_model)
set.seed(2020)
ggplot(data=patient_data[patient_data$Patientnr %in% sample(levels(patient_data$Patientnr), size = 9, replace = F),]) +
  geom_point(aes(x=Time, y=WOMAC_pain)) +
  geom_line(aes(x=Time, y=modelfit)) + facet_wrap(~factor(Patientnr)) +
  theme_bw() + 
  xlab("Time (weeks)") + ylab("WOMAC function")


##########
# Difference in womac pain duloxetine - usual care at 13,26,39,52 weeks
# patient characteristics are same under both duloxetin and usual care 
##########
test_times = c(6, 13, 26, 39, 52)
fixed_effects = fixef(full_model)
for(test_time in test_times){
  contrast = rep(0, length(fixed_effects))
  names(contrast) = names(fixed_effects)
  #number 6 is fixed effect of treatment
  contrast[6] = 1
  #number 19 and 20 are treatment time interactions
  contrast[c(19,20)] = ns(test_time, knots=26, Boundary.knots = c(0,52))
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

##########
# Difference in womac pain duloxetine at 13,26,39,52 weeks and baseline
# The chosen patient's characteristics are set below 
##########
test_times = c(13, 26, 39, 52)
fixed_effects = fixef(full_model)
temp_patient = patient_data[1,]
temp_patient$HADS_depression = median(temp_patient$HADS_depression)
temp_patient$HADS_anxiety = median(temp_patient$HADS_anxiety)
temp_patient$age = median(temp_patient$age)
temp_patient$BMI = median(temp_patient$BMI)
temp_patient$sex = factor("male", levels = c("female", "male"))
temp_patient$joint = factor("heup", levels = c("heup", "knie"))
temp_patient$Treatment = factor("duloxetine", levels = c("usual care", "duloxetine"))
#confirm patient data before proceeding
View(temp_patient)

for(test_time in test_times){
  temp_patient$Time = test_time
  contrast = c(model.matrix(~(BMI + age + sex + joint + Treatment + 
                                HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                            temp_patient))
  contrast[1:8]=0
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

#Final summary plots
#Everything below this point for plotting data
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

#We will make plots for patient with median age, median BMI, median HADS scores, 
#for both joints and both genders, both treatments over time
pred_patient_data = data.frame(expand.grid(Time=seq(0,52,length.out = 10),
                                           Treatment=levels(patient_data$Treatment),
                                           sex=levels(patient_data$sex),
                                           joint=levels(patient_data$joint),
                                           HADS_depression=median(patient_data$HADS_depression),
                                           HADS_anxiety=median(patient_data$HADS_anxiety),
                                           age=median(patient_data$age),
                                           BMI=median(patient_data$BMI)))

plotData = effectPlotData(full_model, newdata = pred_patient_data, patient_data)
#This is the patient data for whom we make a plot
View(plotData)

summary(plotData)

ggplot(data=plotData) + geom_line(aes(x=Time, y=pred, 
                                      color=Treatment, 
                                      group=Treatment)) + 
  geom_ribbon(aes(ymin=low, ymax=upp, x=Time, fill = Treatment), alpha = 0.3)+
  facet_grid(joint~sex)+
  theme_bw() + 
  xlab("Time (weeks)") + ylab("WOMAC function")

# Nominated VAS.

#The unadjusted model. you can replace ns(Time, knots=26) with just Time to remove splines
unadjusted_model = lme(data=patient_data,
                       fixed = NVAS~Treatment* ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                       random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                       na.action = na.omit,
                       control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(unadjusted_model)
#95% confidence interval for fixed effectss
intervals(unadjusted_model, which = "fixed")
#p-values
anova(unadjusted_model)
#Complete summary
summary(unadjusted_model)

#The full model now
full_model = lme(data=patient_data,
                 fixed = NVAS~(BMI + age + sex + joint + Treatment + 
                                            HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                 random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                 na.action = na.omit,
                 control = lmeControl(maxIter = 1500, opt="optim", optimMethod = "L-BFGS-B"))

#Print results
#1. average effect size for each covariate in the model
print(full_model)
#95% confidence interval for fixed effectss
intervals(full_model, which = "fixed")
#p-values
anova(full_model)
#Complete summary
summary(full_model)

#full model, assumption of homoscedasticity of residuals is satisfied
plot(full_model)

#qqplot of residuals made by two lines is also fine
qqnorm(residuals(full_model))
qqline(residuals(full_model))

#Now check model fit for 9 randomly selected patients. 
#please always run the set.seed(2020) before the ggplot code line
patient_data$modelfit = fitted(full_model)
set.seed(2027)
ggplot(data=patient_data[patient_data$Patientnr %in% sample(levels(patient_data$Patientnr), size = 9, replace = F),]) +
  geom_point(aes(x=Time, y=WOMAC_Stiffness)) +
  geom_line(aes(x=Time, y=modelfit)) + facet_wrap(~factor(Patientnr)) +
  theme_bw() + 
  xlab("Time (weeks)") + ylab("NVAS")


##########
# Difference in womac pain duloxetine - usual care at 13,26,39,52 weeks
# patient characteristics are same under both duloxetin and usual care 
##########
test_times = c(13, 52)
fixed_effects = fixef(full_model)
for(test_time in test_times){
  contrast = rep(0, length(fixed_effects))
  names(contrast) = names(fixed_effects)
  #number 6 is fixed effect of treatment
  contrast[6] = 1
  #number 19 and 20 are treatment time interactions
  contrast[c(19,20)] = ns(test_time, knots=26, Boundary.knots = c(0,52))
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

##########
# Difference in womac pain duloxetine at 13,26,39,52 weeks and baseline
# The chosen patient's characteristics are set below 
##########
test_times = c(13, 52)
fixed_effects = fixef(full_model)
temp_patient = patient_data[1,]
temp_patient$HADS_depression = median(temp_patient$HADS_depression)
temp_patient$HADS_anxiety = median(temp_patient$HADS_anxiety)
temp_patient$age = median(temp_patient$age)
temp_patient$BMI = median(temp_patient$BMI)
temp_patient$sex = factor("male", levels = c("female", "male"))
temp_patient$joint = factor("heup", levels = c("heup", "knie"))
temp_patient$Treatment = factor("duloxetine", levels = c("usual care", "duloxetine"))
#confirm patient data before proceeding
View(temp_patient)

for(test_time in test_times){
  temp_patient$Time = test_time
  contrast = c(model.matrix(~(BMI + age + sex + joint + Treatment + 
                                HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                            temp_patient))
  contrast[1:8]=0
  meanEff = contrast %*% fixed_effects
  stdErr = sqrt(c(t(contrast) %*% vcov(full_model) %*% contrast))
  
  print("*************")
  print("*************")
  print(paste0("Time = ", test_time, " weeks"))
  print(paste0("Avg difference = ", meanEff))
  print(paste0("95% confidence interval difference = ",qnorm(p=c(0.025, 0.975), mean=meanEff, sd=stdErr)))
}

#Final summary plots
#Everything below this point for plotting data
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

#We will make plots for patient with median age, median BMI, median HADS scores, 
#for both joints and both genders, both treatments over time
pred_patient_data = data.frame(expand.grid(Time=seq(0,52,length.out = 10),
                                           Treatment=levels(patient_data$Treatment),
                                           sex=levels(patient_data$sex),
                                           joint=levels(patient_data$joint),
                                           HADS_depression=median(patient_data$HADS_depression),
                                           HADS_anxiety=median(patient_data$HADS_anxiety),
                                           age=median(patient_data$age),
                                           BMI=median(patient_data$BMI)))

plotData = effectPlotData(full_model, newdata = pred_patient_data, patient_data)
#This is the patient data for whom we make a plot
View(plotData)

summary(plotData)

ggplot(data=plotData) + geom_line(aes(x=Time, y=pred, 
                                      color=Treatment, 
                                      group=Treatment)) + 
  geom_ribbon(aes(ymin=low, ymax=upp, x=Time, fill = Treatment), alpha = 0.3)+
  facet_grid(joint~sex)+
  theme_bw() + 
  xlab("Time (weeks)") + ylab("WOMAC Stiffness")

#OARSI-OMERACT
#unadjusted.
m<-glmer(FU13_imp ~ Treatment + (1| Cluster), data=patient_data, family=binomial,
control=glmerControl(optimizer="bobyqa"), nAGQ=10)

print (m, corr=FALSE)

se <-sqrt(diag(vcov(m)))
(tab<-cbind(Est=fixef(m), LL=fixef(m)-1.96*se, UL=fixef(m)+1.96*se))

exp (tab)

#adjusted
ma<-glmer(FU13_imp ~ Treatment + age + sex + Score_Paindetect+ HADS_depression + Com_min2 + (1| Cluster), data=patient_data, family=binomial,
         control=glmerControl(optimizer="bobyqa"), nAGQ=10)

print (ma, corr=FALSE)

se1 <-sqrt(diag(vcov(ma)))
(tab<-cbind(Est=fixef(ma), LL=fixef(ma)-1.96*se, UL=fixef(ma)+1.96*se))

exp (tab)
