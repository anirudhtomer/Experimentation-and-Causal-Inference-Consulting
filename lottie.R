library(ggplot2)

data = read.csv(file = file.choose(), header = T)
data$gender = ifelse(data$gender==1, "Male", "Female")

colnames(data)[1] = "PatientCode"

#No effect of gender
ggplot(data=data) + geom_point(aes(x=gender, y=Kernscore))

#Some effect of ga_days
ggplot(data=data) + geom_point(aes(x=ga_days, y=Kernscore)) + geom_smooth(aes(x=ga_days, y=Kernscore))

model = lm(data = data, formula = Kernscore~gender + I(ga_days/7))
summary(model)


#Outcome number 2

#No effect of gender
ggplot(data=data) + geom_point(aes(x=gender, y=ExpresIndex))

#Some effect of ga_days
ggplot(data=data) + geom_point(aes(x=ga_days, y=ExpresIndex)) + geom_smooth(aes(x=ga_days, y=ExpresIndex))

model_linear = lm(data = data, formula = ExpresIndex~gender + I(ga_days/7))
model_nonlinear = lm(data = data, formula = ExpresIndex~gender + ns(ga_days, df = 2))

AIC(model_linear)
AIC(model_nonlinear)

ggplot() + geom_line(aes(x=168:223, y=predict(model_nonlinear, newdata = data.frame(gender="Female", ga_days=168:223)), color="Female"))+
  geom_line(aes(x=168:223, y=predict(model_nonlinear, newdata = data.frame(gender="Male", ga_days=168:223)), color="Male")) +
  geom_vline(xintercept = 207) +
  xlab("ga_days") + ylab("ExpresIndex")

#Outcome 3
model_linear = lm(data = data, formula = RecIndex~gender + I(ga_days/7))
model_nonlinear = lm(data = data, formula = RecIndex~gender + ns(ga_days, df = 2))

AIC(model_linear)
AIC(model_nonlinear)


################### Round 2: Correlation between the outcomes #######
data = read.csv(file=file.choose())
