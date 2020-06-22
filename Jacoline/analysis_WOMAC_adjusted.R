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
library(splines)

if(!require(JointAI)){
  install.packages("JointAI")
}
library(JointAI)

#Now we load the data. A pop will appear. You must load data in SPSS .sav file format only
patient_data = read.spss(file = "Datasetlongformat_WOMAC.sav", to.data.frame = T)

#Data summary
summary(patient_data)

#Treatment has 3 levels in the data, one of them is unnecessary
patient_data$Treatment = droplevels(patient_data$Treatment)
#make usual care as the reference category, so that results are printed for duloxetine
patient_data = within(patient_data, Treatment <- relevel(Treatment, ref = "usual care"))

patient_data = patient_data[complete.cases(patient_data[, c("sex", "joint", "age", "BMI", "Com_min2", "Treatment",
                                                            "HADS_anxiety", "Score_Paindetect", "HADS_depression")]),]

#we make a spline for time effect with knot at median follow up of 26 weeks
print(median(patient_data$Time))

#full_model_imputation = lme_imp(data = patient_data,
#                                fixed = WOMAC_pain~(BMI + age + sex + joint + Treatment + Com_min2 + Score_Paindetect+
#                                                     HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)),
#                                random = ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr,
#                                n.iter=1000)

#The full model now
full_model = lme(data=patient_data,
                 fixed = WOMAC_pain~(Treatment + BMI + age + sex + joint +  Com_min2 + Score_Paindetect + HADS_depression + HADS_anxiety)*ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                 random = list(~1|Cluster, ~1 + ns(Time, knots = 26, Boundary.knots = c(0,52)) | Patientnr), 
                 na.action = na.omit,
                 control = lmeControl(maxIter = 3000, msMaxIter=200, opt="optim", optimMethod = "L-BFGS-B", niterEM = 1000))


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
test_times = c(13, 26, 39, 52)
fixed_effects = fixef(full_model)
for(test_time in test_times){
  contrast = rep(0, length(fixed_effects))
  names(contrast) = names(fixed_effects)
  #number 2 is fixed effect of treatment
  contrast[2] = 1
  #number 13 and 14 are treatment time interactions
  contrast[c(13,14)] = ns(test_time, knots=26, Boundary.knots = c(0,52))
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
temp_patient$Com_min2 = factor("ja", levels = c("ja", "nee"))
temp_patient$Score_Paindetect = median(temp_patient$Score_Paindetect)
#confirm patient data before proceeding
View(temp_patient)

for(test_time in test_times){
  temp_patient$Time = test_time
  contrast = c(model.matrix(~(Treatment + BMI + age + sex + joint + Com_min2 + Score_Paindetect+
                                HADS_depression + HADS_anxiety) * ns(Time, knots = 26, Boundary.knots = c(0,52)), 
                            temp_patient))
  #set all fixed effects to 0,  except for effect of time and all the interactions with time
  contrast[1:10]=0
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
                                           Com_min2 = "ja",
                                           Score_Paindetect = median(patient_data$Score_Paindetect),
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
  xlab("Time (weeks)") + ylab("WOMAC Pain (Com_min2 = 'ja')")
