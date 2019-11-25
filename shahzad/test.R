#install.packages("coxed")

#library(coxed)
set.seed(2019)

N=30000
alzeihmer = sample(c(0,1), size = N, replace = T, prob=c(0.9, 0.1))
biomarker = sapply(alzeihmer, function(x){
  rnorm(n = 1, ifelse(x, 1, 5), ifelse(x, 0.5, 2))
})

df = data.frame(alzeihmer, biomarker)

simdata = sim.survdata(N=N, T=1000, X = df,
                       num.data.frames=1,
                       hazard.fun = function(t){
                         (4/3) * (t/3)^(4-1)
                       })
simdata$betas
exp(simdata$betas)

#fit <- survfit(Surv(y, failed) ~ alzeihmer, data = simdata$data)
#survminer::ggsurvplot(fit)

summary(survreg(Surv(y, failed) ~ alzeihmer + biomarker, 
                data=simdata$data, dist='weibull'))

summary(coxph(Surv(y, failed) ~ alzeihmer + biomarker, 
              data=simdata$data))

alzeihmerRows = (1:N)[simdata$data$alzeihmer==T]
nonAlzeihmerRows = sample((1:N)[simdata$data$alzeihmer==F],
                          size = length(alzeihmerRows), replace = F)

df2 = simdata$data[c(alzeihmerRows, nonAlzeihmerRows),]
summary(coxph(Surv(y, failed) ~ alzeihmer + biomarker, 
              data=df2))
