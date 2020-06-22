#First we load packages and install missing ones
if(!require(foreign)){
  install.packages("foreign")
}
if(!require(GLMMadaptive)){
  install.packages("GLMMadaptive")
}
if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(splines)

#Now we load the data. A pop will appear. You must load data in SPSS .sav file format only
patient_data = read.spss(file = "Datasetlongformat_improvement.sav", to.data.frame = T)

#Data summary
summary(patient_data)

#Treatment has 3 levels in the data, one of them is unnecessary
patient_data$Treatment = droplevels(patient_data$Treatment)
#make usual care as the reference category, so that results are printed for duloxetine
patient_data = within(patient_data, Treatment <- relevel(Treatment, ref = "usual care"))

patient_data = patient_data[complete.cases(patient_data[, c("Improvement", "sex", "joint", "age", "BMI", "Com_min2", "Treatment",
                                                            "Score_Paindetect")]),]

patient_data$Improvement = ifelse(patient_data$Improvement=="geen klachten", yes = 1, no = 0)

full_model = mixed_model(data=patient_data,
                         fixed = Improvement~(BMI + age + sex + joint + Treatment + Com_min2 + Score_Paindetect) * Time, 
                         random = ~1 + Time | Patientnr,
                         family = 'binomial', control=list(iter_EM=200))

summary(full_model)

pred_patient_data = data.frame(expand.grid(Time=seq(0,52,length.out = 10),
                                           Treatment=levels(patient_data$Treatment),
                                           sex=levels(patient_data$sex),
                                           joint=levels(patient_data$joint),
                                           Com_min2 = "ja",
                                           Score_Paindetect = median(patient_data$Score_Paindetect),
                                           age=median(patient_data$age),
                                           BMI=median(patient_data$BMI)))

plotData = effectPlotData(full_model, newdata = pred_patient_data)
#This is the patient data for whom we make a plot
View(plotData)

summary(plotData)

ggplot(data=plotData) + geom_line(aes(x=Time, y=plogis(pred), 
                                      color=Treatment, 
                                      group=Treatment)) + 
  geom_ribbon(aes(ymin=plogis(low), ymax=plogis(upp), x=Time, fill = Treatment), alpha = 0.3)+
  facet_grid(joint~sex)+
  theme_bw() + 
  xlab("Time (weeks)") + ylab("Patient improvement log odds (Com_min2 = 'ja')")

ggplot(data=plotData) + geom_line(aes(x=Time, y=pred, 
                                      color=Treatment, 
                                      group=Treatment)) + 
  geom_ribbon(aes(ymin=low, ymax=upp, x=Time, fill = Treatment), alpha = 0.3)+
  facet_grid(joint~sex)+
  theme_bw() + 
  xlab("Time (weeks)") + ylab("Patient improvement probability (Com_min2 = 'ja')")
