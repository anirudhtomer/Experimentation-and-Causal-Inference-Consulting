library(nlme)
library(ggplot2)

#After you load the Rdata object run this
DF = DF[!is.na(DF$methylation),]
DF$time = DF$time / 12

ggplot(data=DF) + geom_line(aes(x=time,y=methylation, group=id)) + 
  facet_grid(temperature~material)


#Random intercept
model1 = lme(data=DF, 
             fixed = methylation ~ material + time * temperature + 
               material * time + material * temperature, 
             random = ~1|id,
             control = lmeControl(opt = "optim", optimMethod = "BFGS"))

#Random intercept + slope
model2 = lme(data=DF, 
             fixed = methylation ~ material + time * temperature + 
               material * time + material * temperature, 
             random = ~time|id,
             control = lmeControl(opt = "optim", optimMethod = "BFGS"))

#Random intercept + slope but nested inside each patient.
#So basically each temperature setting has its own random intercept and slope but nested per patient
model3 = lme(data=DF, 
             fixed = methylation ~ material + time * temperature + 
               material * time + material * temperature, 
             random = ~time|id/temperature,
             control = lmeControl(opt = "optim", optimMethod = "BFGS"))

#See some parameters of the covariance matrix are too small
summary(model3)

anova(model1, model2, model3)
AIC(model1, model2, model3)
#I will prefer model2 coz later after 24 months I really do expect random slopes

#F test for each of the fixed effect parameters in model2 shows that all parameters are necessary
# Effect of material is close to none, but lets still keep
anova(model2)

#Remove interaction of material and temperature
model_final = lme(data=DF, 
                  fixed = methylation ~ material + time * temperature + 
                    material * time, 
                  random = ~time|id,
                  control = lmeControl(opt = "optim", optimMethod = "BFGS"))

##########
#Final Model diagnostics based on subject specific fits and residuals
#########

DF$fitted = model_final$fitted[,2]
DF$residuals = model_final$residuals[,2]

#QQplot is fine
qqnorm(DF$residuals)
qqline(DF$residuals)

#Now fitted profiles for "+4C"
ggplot(data=DF[DF$temperature=="+4C",]) + geom_line(aes(x=time,y=fitted, color=material)) + 
  geom_point(aes(x=time,y=methylation, shape=material,color=material), size=3) + theme_bw() +
  theme(text=element_text(size=15)) + xlab("Time (years)") + ylab("Methylation") + 
  facet_wrap(~id)

#Now fitted profiles for "-20C"
ggplot(data=DF[DF$temperature=="-20C",]) + geom_line(aes(x=time,y=fitted,  color=material)) + 
  geom_point(aes(x=time,y=methylation, shape=material,color=material), size=3) + theme_bw() +
  theme(text=element_text(size=15)) + xlab("Time (years)") + ylab("Methylation") + 
  facet_wrap(~id)

#Now fitted profiles for "-80C"
ggplot(data=DF[DF$temperature=="-80C",]) + geom_line(aes(x=time,y=fitted, color=material)) + 
  geom_point(aes(x=time,y=methylation, shape=material,color=material), size=3) + theme_bw() +
  theme(text=element_text(size=15)) + xlab("Time (years)") + ylab("Methylation") + 
  facet_wrap(~id)

##############
#Question: Are the profiles stable over time. 
#To answer this question. Remove the effect of time completely
##############
model_test_with_time = update(model_final, method = "ML")

model_test_without_time = lme(data=DF, 
                              fixed = methylation ~ material + temperature,
                              random = ~1|id,
                              control = lmeControl(opt = "optim", optimMethod = "BFGS"),
                              method = "ML")

AIC(model_test_without_time, model_test_with_time)
anova(model_test_without_time, model_test_with_time)
#So answer to your question: Yes effect of time is necessary in our model

################
# Predictions with 95% CI
###############
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

predDf = data.frame(expand.grid(id=1,material=c("DNA", "Blood"), 
                                temperature=c("+4C", "-20C", "-80C"), 
                                time=seq(0,1, 0.25)))
predDf = effectPlotData(model_final, predDf, DF)

#Predictions with 95% CI
ggplot(data=predDf) + geom_line(aes(x=time, y=pred)) + 
  geom_ribbon(aes(x=time, ymin=low, ymax=upp), fill='grey', alpha=0.5) +
  facet_grid(material~temperature) + theme_bw() +
  theme(text=element_text(size=15)) + xlab("Time (years)") + ylab("Methylation")
