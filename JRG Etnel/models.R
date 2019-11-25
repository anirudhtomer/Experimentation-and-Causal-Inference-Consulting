model_all = lme(data=ds, fixed=peak.rvot.grad ~ male + 
                ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) *age + 
                diagnosis + etiology.new + nyha.binom + sah + 
                congenital.anomaly + smoking + 
                previos.operation.yn + associated.pathologies.yn + 
                urgent + concomitant.procedures.yn + 
                ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) *pulmonary.allograft +
                ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) *rvot.augmentation.yn +
                ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) *allograft.diameter + donor.age +
                donor.male + ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) *abo.mismatch + 
                bmi + 
                ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) * dcell,
              random = ~ns(fup.duration, df=3, Boundary.knots = c(0, 6.7))|id,
              control = lmeControl(opt = "optim"), method="ML")

model_3_usefulInt = lme(data=ds, fixed=peak.rvot.grad ~ male + 
                    ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) *age + 
                    diagnosis + etiology.new + nyha.binom + sah + 
                    congenital.anomaly + smoking + 
                    previos.operation.yn + associated.pathologies.yn + 
                    urgent + concomitant.procedures.yn + 
                    pulmonary.allograft +
                     ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) * rvot.augmentation.yn +
                    allograft.diameter + donor.age +
                    donor.male + 
                    ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) * abo.mismatch + bmi + 
                    ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) * dcell,
                  random = ~ns(fup.duration, df=3, Boundary.knots = c(0, 6.7))|id,
                  control = lmeControl(opt = "optim"), method="ML")

#A simpler final model if the interest in only in interaction with dcell
model_final = lme(data=ds, fixed=peak.rvot.grad ~ previos.operation.yn + 
                pulmonary.allograft + donor.age +
                ns(fup.duration, df=3, Boundary.knots = c(0, 6.7)) * dcell,
              random = ~ns(fup.duration, df=3, Boundary.knots = c(0, 6.7))|id,
              control = lmeControl(opt = "optim"), method="ML")

#Linear evolution like the Brazilian team
model_linear = lme(data=ds, fixed=peak.rvot.grad ~ 
                  fup.duration * dcell,
                random = ~fup.duration|id,
                control = lmeControl(opt = "optim"), method="REML")

anova(model_1, model_2, model_3, model_4)

anova(model_4, Terms=c("ns(fup.duration, df = 3, Boundary.knots = c(0, 15.58)):dcell"))
anova(model_3_int, Terms=c("ns(fup.duration, df = 3, Boundary.knots = c(0, 15.58)):allograft.diameter"))
anova(model_4, model_5)

