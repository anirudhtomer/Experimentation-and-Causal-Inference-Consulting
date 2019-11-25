library(foreign)
library(ggplot2)
library(nlme)

clair = read.spss(file.choose(), to.data.frame=TRUE)
clair_alod = clair[,1:7]
clair_alod$agedecade=clair_alod$age/10
clair_alod = clair_alod[complete.cases(clair_alod),]

clair_alod = cbind(clair_alod, mcage = scale(clair_alod$agedecade,center = T, scale = F), 
                   mccaod=scale(clair_alod$CA_OD, center = T, scale = F))

ggplot(data=clair_alod[1:200,], aes(y=AL_OD, x=age)) + geom_line(aes(group=IDC))

model_0 = lme(fixed = AL_OD~age + CA_OD + ETHNI_EUR + GENDER,
              random=~I(age/10)|IDC, data=clair_alod, method="ML")

model_1 = lme(fixed = AL_OD~I(age/10) + CA_OD + ETHNI_EUR + GENDER +
            I(age/10)*ETHNI_EUR + I(age/10)*GENDER + CA_OD*ETHNI_EUR +
              CA_OD *GENDER + I(age/10)*CA_OD,
    random=~I(age/10)|IDC, data=clair_alod, method="ML")

model_2 = lme(fixed = AL_OD~I(age/10) + CA_OD + ETHNI_EUR + GENDER +
                I(age/10)*ETHNI_EUR + I(age/10)*CA_OD,
              random=~I(age/10)|IDC, data=clair_alod, method="ML")


model_3 = lme(fixed = AL_OD~mcage + mccaod + ETHNI_EUR + GENDER +
                mcage*ETHNI_EUR + mcage*mccaod+ ETHNI_EUR * GENDER,
              random=~mcage|IDC, data=clair_alod, method="ML")

#Model 4 is grand (AIC is almost same)
model_4 = lme(fixed = AL_OD~mcage + mccaod + ETHNI_EUR + GENDER +
                mcage*ETHNI_EUR + mcage*mccaod,
              random=~mcage|IDC, data=clair_alod, method="ML")

anova(model_4, type="marginal")

#For the other response variable
clair_alodcr = clair[,c(1:6,8)]
clair_alodcr$agedecade=clair_alodcr$age/10
clair_alodcr = clair_alodcr[complete.cases(clair_alodcr),]
clair_alodcr = cbind(clair_alodcr, mcage = scale(clair_alodcr$agedecade,center = T, scale = F), 
                   mccaod=scale(clair_alodcr$CA_OD, center = T, scale = F))

ggplot(data = clair_alodcr, aes(x=mccaod, y=ALCR_OD)) + geom_line(aes(group=IDC)) + stat_smooth()

model_cr_1 = lme(fixed = ALCR_OD~mcage * mccaod + mcage*ETHNI_EUR + GENDER,
    random=~mcage|IDC, data=clair_alodcr)

summary(model_cr_1)
anova(model_cr_1, type="marginal")

###############################
#@Clair: Everything that you need begins here.
################################
clair$IDC = factor(clair$IDC)
clair$IDC=droplevels(clair$IDC)

clair$slope_alod = unlist(by(clair$IDC, data = clair, FUN = function(ds){
    if(nrow(ds)==2){
      slope = (ds$AL_OD[2]-ds$AL_OD[1])/(ds$age[2]-ds$age[1])
      return(rep(slope,2))
    }else{
      return(NA)
    }
  }))

clair$slope_alcr_od = unlist(by(clair$IDC, data = clair, FUN = function(ds){
  if(nrow(ds)==2){
    slope = (ds$ALCR_OD[2]-ds$ALCR_OD[1])/(ds$age[2]-ds$age[1])
    return(rep(slope,2))
  }else{
    return(NA)
  }
}))

clair_first_visit_slopealod = clair[complete.cases(clair[,-8]) & clair$Visit==1, ]

ggplot(data = clair_first_visit_slopealod, aes(y=slope_alod, x=AL_OD)) + 
  geom_point(aes(color=CA_OD>1.25, alpha=1/10)) + 
  stat_smooth(aes(group=CA_OD>1.25, color=CA_OD>1.25)) + 
  facet_grid(GENDER~ETHNI_EUR) + xlim(19, 25) + ylab("Slope(AL)")

ggplot(data = clair_first_visit_slopealod, aes(y=slope_alod, x=CA_OD)) + 
  geom_point(aes(color=AL_OD, alpha=1/10)) + 
  stat_smooth() + 
  facet_grid(GENDER~ETHNI_EUR) + ylab("Slope(AL)")

clair_first_visit_slopealod_alcr = clair[complete.cases(clair[,-7]) & clair$Visit==1, ]

ggplot(data = clair_first_visit_slopealod_alcr, aes(y=slope_alod, x=ALCR_OD)) + 
  geom_point(aes(color=CA_OD>1.25, alpha=1/10)) + 
  stat_smooth(aes(group=CA_OD>1.25, color=CA_OD>1.25)) + 
  facet_grid(GENDER~ETHNI_EUR) + xlim(2.5, 3.125) + ylab("Slope(AL)")

ggplot(data = clair_first_visit_slopealod_alcr, aes(y=slope_alod, x=CA_OD)) + 
  geom_point(aes(color=ALCR_OD, alpha=1/10)) + 
  stat_smooth() + 
  facet_grid(GENDER~ETHNI_EUR) + ylab("Slope(AL)")

model_slope_1 = lm(slope_alod~AL_OD + CA_OD + ETHNI_EUR + GENDER, 
                   data = clair_first_visit_slopealod[clair_first_visit_slopealod$AL_OD>=19 & clair_first_visit_slopealod$AL_OD<=25,])
summary(model_slope_1)


model_slope_2 = lm(slope_alod~ALCR_OD + CA_OD + ETHNI_EUR + GENDER, 
                   data = clair_first_visit_slopealod_alcr[clair_first_visit_slopealod_alcr$ALCR_OD>=2.5 & clair_first_visit_slopealod_alcr$ALCR_OD<=3.125,])
summary(model_slope_2)
