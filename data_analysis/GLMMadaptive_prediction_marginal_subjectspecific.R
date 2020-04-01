dataOA_ingrid$scaled_vit_D = scale(dataOA_ingrid$vit.D, center = T, scale = T)

#Please do not use GEE. It works under MCAR which is never the case.

#Four GLMM models are of interest
#1. model with additive terms and linear vitamin D
#2. model with additive terms and spline vitamin D
#3. model with all first order pairwise interactions and spline vitamin D
#4. model with all first order pairwise interactions and linear vitamin D


model1 = glmer(value ~ (sex + group + kneeOA_incidence + scaled_vit_D + sc.age +sc.BMI + sc.months) + (1|id), data=dataOA_ingrid, family=binomial, nAGQ=10, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
model4 = glmer(value ~ (sex + group + kneeOA_incidence + scaled_vit_D + sc.age +sc.BMI + sc.months)^2 + (1|id), data=dataOA_ingrid, family=binomial, nAGQ=10, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model2 = glmer(value ~ (sex + group + kneeOA_incidence + ns(scaled_vit_D, df=3) + sc.age +sc.BMI + sc.months) + (1|id), data=dataOA_ingrid, family=binomial, nAGQ=10, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=5e5)))
model3 = glmer(value ~ (sex + group + kneeOA_incidence + ns(scaled_vit_D, df=3) + sc.age +sc.BMI + sc.months)^2 + (1|id), data=dataOA_ingrid, family=binomial, nAGQ=10, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=5e5)))

AIC(model1, model2, model3, model4)

#Model 4 has the least AIC, however it is only 10 units smaller. 
#Model 1, 2, 3 are all same. Better use model 1 because it is simpler

#We compare model 1 and model4
library(ggplot2)
library(ggpubr)
fit_model1 = ggplot() + geom_point(aes(x=fitted(model1),y=dataOA_ingrid$value)) + 
  geom_smooth(aes(x=fitted(model1),y=dataOA_ingrid$value)) + geom_abline(slope = 1, intercept = 0) + 
  xlab("Fitted Probability") + ylab("Observed Value") + scale_y_continuous(breaks = c(0,1)) + 
  ggtitle("Model-1")

fit_model4 = ggplot() + geom_point(aes(x=fitted(model4),y=dataOA_ingrid$value)) + 
  geom_smooth(aes(x=fitted(model4),y=dataOA_ingrid$value)) + geom_abline(slope = 1, intercept = 0) + 
  xlab("Fitted Probability") + ylab("Observed Value") + scale_y_continuous(breaks = c(0,1)) + 
  ggtitle("Model-4")

ggpubr::ggarrange(fit_model1, fit_model4, ncol=2, nrow=1)

#Both are terrible. I would still choose model1 for its simplicity

#To make subject specific predictions, I use Dimitris' package
#https://drizopoulos.github.io/GLMMadaptive/index.html
library(GLMMadaptive)

model1_dimitris = mixed_model(fixed = value ~  sex + group + kneeOA_incidence + scaled_vit_D + sc.age +sc.BMI + sc.months, 
                              random = ~ 1 | id, data = dataOA_ingrid, 
                              family = binomial(), iter_EM = 0, 
                              optimizer = "nlminb", iter_qN_incr = 5, 
                              initial_values = list(betas=summary(model1)$coefficients[,1], D=c(summary(model1)$varcor[[1]])))

##################################
# Subject specific predictions 
##################################
#Choose a patient's training data
training_data_patient_12 = dataOA_ingrid[dataOA_ingrid$id %in% c(12),]
#Choose some vitamin D levels for which you want prediction. 
#I generate a sequence of them based on the quantiles in the original dataset
vit.D_levels = seq(quantile(dataOA_ingrid$vit.D, probs=0.05), quantile(dataOA_ingrid$vit.D, probs=0.95), length.out = 100)
#Scale them as per the scale of training data
scaled_vitD_levels = (vit.D_levels - mean(dataOA_ingrid$vit.D))/(sd(dataOA_ingrid$vit.D))

#Creata a prediction dataset. I choose patient 12 and his right knee for illustration
prediction_data_patient_12 = dataOA_ingrid[dataOA_ingrid$id==12 & dataOA_ingrid$kneeOA_incidence %in% c("kneeOA_R_incidence"),]
prediction_data_patient_12 = prediction_data_patient_12[rep(1, length(scaled_vitD_levels)),]
#add vitamin D levels to this patient
prediction_data_patient_12$vit.D = vit.D_levels
prediction_data_patient_12$scaled_vit_D = scaled_vitD_levels

#predict the probability of outcome
prediction_result = predict(model1_dimitris, training_data_patient_12, prediction_data_patient_12,
                            type_pred = "response", 
                            type = c("subject_specific"), 
                            se.fit = T,level = 0.95, return_newdata = T)

#You can change the ylim to zoom into the graph
ggplot(data=prediction_result$newdata2) + geom_line(aes(x=vit.D, y=pred)) + 
  geom_ribbon(aes(x=vit.D, ymin=low, ymax=upp), fill="grey", alpha=0.5) + 
  ylim(0,1) + xlab("Vitamin D") + ylab("Probability of Outcome = 1")

#Risk ratio
reference_vit_D_level = median(unique(dataOA_ingrid$vit.D))
reference_vit_D_level_scaled = (reference_vit_D_level - mean(dataOA_ingrid$vit.D))/(sd(dataOA_ingrid$vit.D))
reference_data_patient_12 = dataOA_ingrid[dataOA_ingrid$id==12 & dataOA_ingrid$kneeOA_incidence %in% c("kneeOA_R_incidence"),]
#add vitamin D levels to this patient
reference_data_patient_12$vit.D = reference_vit_D_level
reference_data_patient_12$scaled_vit_D = reference_vit_D_level_scaled
reference_prediction = predict(model1_dimitris, training_data_patient_12, reference_data_patient_12,
        type_pred = "response", 
        type = c("subject_specific"), 
        se.fit = T,level = 0.95, return_newdata = T)$newdata2$pred


#Risk ratio plot
ggplot() + geom_line(aes(x=prediction_result$newdata2$vit.D, y=prediction_result$newdata2$pred/reference_prediction)) + 
  xlab("Vitamin D") + ylab("Risk ratio (Outcome = 1)") + ggtitle("Reference Vit D = 55.605")


###########################################
# Marginal predictions
###########################################
#No training data is required, only prediction data is needed
#Choose some vitamin D levels for which you want prediction. 
#I generate a sequence of them based on the quantiles in the original dataset
vit.D_levels = seq(quantile(dataOA_ingrid$vit.D, probs=0.05), quantile(dataOA_ingrid$vit.D, probs=0.95), length.out = 100)
#Scale them as per the scale of training data
scaled_vitD_levels = (vit.D_levels - mean(dataOA_ingrid$vit.D))/(sd(dataOA_ingrid$vit.D))

#Creata a prediction dataset. I choose patient 12 and his right knee for illustration
prediction_data_patient_12 = dataOA_ingrid[dataOA_ingrid$id==12 & dataOA_ingrid$kneeOA_incidence %in% c("kneeOA_R_incidence"),]
prediction_data_patient_12 = prediction_data_patient_12[rep(1, length(scaled_vitD_levels)),]
#add vitamin D levels to this patient
prediction_data_patient_12$vit.D = vit.D_levels
prediction_data_patient_12$scaled_vit_D = scaled_vitD_levels

#predict the probability of outcome
prediction_result = predict(model1_dimitris,prediction_data_patient_12,
                            type_pred = "response", 
                            type = c("marginal"), 
                            se.fit = T,level = 0.95, return_newdata = F)

prediction_result_df = prediction_data_patient_12
prediction_result_df$pred = plogis(prediction_result$pred)
prediction_result_df$upp = plogis(qnorm(0.975, prediction_result$pred, prediction_result$se.fit))
prediction_result_df$low = plogis(qnorm(0.025, prediction_result$pred, prediction_result$se.fit))

#You can change the ylim to zoom into the graph
ggplot(data=prediction_result_df) + geom_line(aes(x=vit.D, y=pred)) + 
  geom_ribbon(aes(x=vit.D, ymin=low, ymax=upp), fill="grey", alpha=0.5) + 
  ylim(0,1) + xlab("Vitamin D") + ylab("Probability of Outcome = 1")

#Risk ratio
reference_vit_D_level = median(unique(dataOA_ingrid$vit.D))
reference_vit_D_level_scaled = (reference_vit_D_level - mean(dataOA_ingrid$vit.D))/(sd(dataOA_ingrid$vit.D))
reference_data_patient_12 = dataOA_ingrid[dataOA_ingrid$id==12 & dataOA_ingrid$kneeOA_incidence %in% c("kneeOA_R_incidence"),]
#add vitamin D levels to this patient
reference_data_patient_12$vit.D = reference_vit_D_level
reference_data_patient_12$scaled_vit_D = reference_vit_D_level_scaled
reference_prediction = plogis(predict(model1_dimitris, reference_data_patient_12,
                               type_pred = "response", 
                               type = c("marginal"), 
                               se.fit = T,level = 0.95, return_newdata = F)$pred)


#Risk ratio plot
ggplot() + geom_line(aes(x=prediction_result_df$vit.D, y=prediction_result_df$pred/reference_prediction)) + 
  xlab("Vitamin D") + ylab("Risk ratio (Outcome = 1)") + ggtitle("Reference Vit D = 55.605")
