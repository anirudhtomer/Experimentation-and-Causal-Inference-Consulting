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

#put list of covariates here that are in the model
newDF <- with(ds, expand.grid(fup.duration = seq(0, 10, by = 0.25),
                                previos.operation.yn = factor(c(1)),
                                pulmonary.allograft = factor(c(0)),
                              rvot.augmentation.yn = factor(0),
                              abo.mismatch  =factor(1),
                                donor.age = 38,
                                dcell = factor(c(0,1))))
#Put model object name here
temp = effectPlotData(model_3_usefulInt, newDF, ds)

ggplot(data=temp) + geom_line(aes(y=pred, x=fup.duration, color=dcell, group=dcell)) + 
    geom_ribbon(aes(x=fup.duration, ymin=low, ymax=upp, group=dcell), fill = "gray", alpha=0.4) + 
    xlab("Follow-up time (Years)") + scale_y_continuous(breaks=seq(0, 200, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 20, by=2)) + 
    ylab("peak.rvot.grad") + theme(text = element_text(size=13))
