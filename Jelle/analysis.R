#put list of covariates here that are in the model
newDF <- with(ds, expand.grid(age = 68,
                              sex = "male",
                              ASA.dich = "1 - 2", 
                              Time=seq(1, 7, 0.25)))
#Put model object name here
temp = effectPlotData(model, newDF, ds)

ggplot(data=temp) + geom_line(aes(y=pred, x=Time)) + 
  geom_ribbon(aes(x=Time, ymin=low, ymax=upp), fill = "gray", alpha=0.4) + 
  xlab("Follow-up time (Years)") +
  scale_x_continuous(breaks = seq(0, 20, by=2)) + 
  ylab("peak.rvot.grad") + theme(text = element_text(size=13))
ds = read.spss(file = file.choose(), to.data.frame = T)
ggplot(data=ds, aes(x = Time, y=CRP)) + geom_line(aes(group=PID)) + stat_summary(fun.y=mean, colour="red", geom="line", size = 3)

model_1 = lme(data=ds[!is.na(ds$CRP),], fixed=CRP ~ ns(Time, df=3) + age + sex +ASA.dich,
random = ~ns(Time, df=3)|PID,
control = lmeControl(opt = "optim"), method="ML")

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
