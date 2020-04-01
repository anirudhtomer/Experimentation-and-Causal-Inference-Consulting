#First we automatically check and install packages that are needed for analysis
if(!require(xlsx)){
  install.packages("xlsx")
}
if(!require(nlme)){
  install.packages("nlme")
}

if(!require(ggplot2)){
  install.packages("ggplot2")
}

library(xlsx)
library(nlme)
library(ggplot2)

#Read the data from the excel file format.
#a windows file explorer dialog will open, you can choose a file it.
df = read.xlsx(file=file.choose(), sheetIndex = 1, as.data.frame = T)
#the data set that will load is called "df"
#Set the reference category for Animal Type, side of brain, and age group
df$Animal.Type = relevel(x = df$Animal.Type, ref = "WTw")
df$Side.of.Brain = relevel(x = df$Side.of.Brain, ref = "A")
df$Age.Group = relevel(x = df$Age.Group, ref = "Young")

#we create a model with all possible interactions
full_model1 = lme(Lipid..2..Sphingosine. ~ Side.of.Brain * Age.Group * Animal.Type, 
                 data = df, random = ~ 1|Animal.ID, 
                 weights = varIdent(form=~1))

#a model where variance of error depends on categorical subgroups
full_model2 = lme(Lipid..2..Sphingosine. ~ Side.of.Brain * Age.Group * Animal.Type, 
            data = df, random = ~ 1|Animal.ID, 
            weights = varIdent(form=~1|Animal.Type*Age.Group*Side.of.Brain),
            control = lmeControl(maxIter = 500, opt="optim", optimMethod = "L-BFGS-B"))

anova(full_model1, full_model2)

#pairwise t-test p-values for comparing 15 combinations 
#with the reference category combination that we set earlier after loading data
summary(full_model)

#F-test to test each of the main effects and interactions overall
anova(full_model)

#F test to get a p-value for side of brain
anova(full_model, Terms = c("Side.of.Brain", "Side.of.Brain:Age.Group", "Side.of.Brain:Animal.Type", "Side.of.Brain:Age.Group:Animal.Type"))

#F test to get a p-value for age group
anova(full_model, Terms = c("Age.Group", "Side.of.Brain:Age.Group", "Age.Group:Animal.Type", "Side.of.Brain:Age.Group:Animal.Type"))

#F test to get a p-value for animal type
anova(full_model, Terms = c("Animal.Type", "Side.of.Brain:Animal.Type", "Age.Group:Animal.Type", "Side.of.Brain:Age.Group:Animal.Type"))

#check the model assumptions for homoscedasticity of error variance
plot(full_model)
#qqplot for normality assumption
qqnorm(residuals(full_model))
qqline(residuals(full_model))

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

pred_df = data.frame(expand.grid(Animal.Type=levels(df$Animal.Type),
                               Side.of.Brain=levels(df$Side.of.Brain),
                                 Age.Group=levels(df$Age.Group)))

plotData = effectPlotData(full_model2, newdata = pred_df, df)
plotData = effectPlotData(ttmodel, newdata = pred_df, tt)

ggplot(data=plotData) + geom_point(aes(x=Side.of.Brain, y=pred, 
                                       color=Age.Group, group=Age.Group)) + 
  geom_errorbar(aes(x=Side.of.Brain, ymin=low, ymax=upp, 
                    color=Age.Group, group=Age.Group), width=0.1) + 
  facet_grid(.~Animal.Type) + theme_bw() + 
  xlab("Side of brain") + ylab("Outcome") +
  scale_y_continuous(breaks = 1:25, limits = c(1,25))
