library(foreign)
library(ggplot2)
dataset = read.spss(file.choose(), to.data.frame = T)

ggplot(data=dataset) + 
  geom_boxplot(aes(x=Sex, y=DP_amplitude), outlier.shape = NA) +
  facet_grid(Ear~Frequency)+
  theme_bw()

dataset$Frequency = factor(dataset$Frequency)

#Nested random effects
library(nlme)
model_nested = lme(fixed=DP_amplitude~Ear * Sex, random = list(~1|id, ~1|Frequency), 
            data = dataset, control = lmeControl(optimMethod = "L-BFGS-B"))

summary(model_nested)
anova(model_nested)

#Crossed random effects
library(lme4)
model_crossed = lmer(DP_amplitude~Ear * Sex + (1 | id) + (1 |Frequency), data = dataset)
summary(model_crossed)
anova(model_crossed)
