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

#Now we load the data. A pop will appear. You must load data in SPSS .sav file format only
patient_data = read.spss(file = "Datasetlongformat_NVAS.sav", to.data.frame = T)

#Data summary
summary(patient_data)

#Treatment has 3 levels in the data, one of them is unnecessary
patient_data$Treatment = droplevels(patient_data$Treatment)
#make usual care as the reference category, so that results are printed for duloxetine
patient_data = within(patient_data, Treatment <- relevel(Treatment, ref = "usual care"))

patient_data = patient_data[complete.cases(patient_data[, c("NVAS", "sex", "joint", "age", "BMI", "Com_min2", "Treatment",
                                                             "Score_Paindetect")]),]

full_model = lme(data=patient_data,
                 fixed = NVAS~(BMI + age + sex + joint + Treatment + Com_min2 + Score_Paindetect) * Time, 
                 random = list(~1|Cluster, ~1 + Time | Patientnr), 
                 na.action = na.omit,
                 control = lmeControl(maxIter = 5000, msMaxIter=500, opt="optim", optimMethod = "BFGS"))

anova(full_model)
summary(full_model)


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
  xlab("Time (weeks)") + ylab("NVAS (Com_min2 = 'ja')")
