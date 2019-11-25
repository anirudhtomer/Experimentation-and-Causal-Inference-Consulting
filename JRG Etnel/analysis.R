library(xlsx)
library(nlme)
library(splines)
library(ggplot2)

ds = read.xlsx(file.choose(), sheetIndex = 1)

ds$dcell = as.factor(ds$dcell)
ds$male = as.factor(ds$male)
ds$diagnosis = as.factor(ds$diagnosis)
ds$etiology.new  = as.factor(ds$etiology.new)
ds$nyha.binom  = as.factor(ds$nyha.binom)
ds$sah  = as.factor(ds$sah)
ds$congenital.anomaly  = as.factor(ds$congenital.anomaly)
ds$smoking  = as.factor(ds$smoking)
ds$previos.operation.yn  = as.factor(ds$previos.operation.yn)
ds$associated.pathologies.yn  = as.factor(ds$associated.pathologies.yn)
ds$copd  = as.factor(ds$copd)
ds$previous.cv.events.yn  = as.factor(ds$previous.cv.events.yn)
ds$urgent  = as.factor(ds$urgent)
ds$concomitant.procedures.yn  = as.factor(ds$concomitant.procedures.yn)
ds$rvot.augmentation.yn  = as.factor(ds$rvot.augmentation.yn)
ds$pulmonary.allograft  = as.factor(ds$pulmonary.allograft)
ds$abo.mismatch  = as.factor(ds$abo.mismatch)
ds$donor.male  = as.factor(ds$donor.male)
ds$bmi = as.numeric(ds$bmi)

#Up to here you run the code to clean the data set

ds.id = ds[!duplicated(ds$id),]

#Assign pair-id
ds.id$pairId= NA
pairIdCounter = 1
for(i in 1:nrow(ds.id)){
  #check if paired patient has a pairId, else give a new one
  matchedId = ds.id$id.matched.pair[i]
  
  if(matchedId %in% ds.id$id){
    matchedPairId = ds.id$pairId[ds.id$id == matchedId]
  }else{
    matchedPairId = NA
  }
  
  if(is.na(matchedPairId)){
    matchedPairId = pairIdCounter
    pairIdCounter = pairIdCounter + 1
  }
  
  ds.id$pairId[i] = matchedPairId
}

ds$pairId = rep(ds.id$pairId, c(by(ds, ds$id, nrow)))


#Overall pattern
ggplot(data = ds, aes(x = fup.duration, y = (peak.rvot.grad))) + 
    geom_line(aes(group=id)) + geom_smooth() + 
    scale_x_continuous(breaks = seq(0, 30, by = 1))

#Pattern for individuals
ggplot(data = ds[ds$id==sample(unique(ds$id),size = 1),], aes(x = fup.duration, y = peak.rvot.grad)) + 
    geom_line() + geom_point()  + scale_x_continuous(breaks = seq(0, 30, by = 1)) + 
    geom_vline(xintercept = 1)

#Overall pattern by dcell
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=dcell, color=factor(dcell))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by male
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=male, color=factor(male))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1)) 

#Overall pattern by diagnosis
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=diagnosis, color=factor(diagnosis))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by etiology
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=etiology.new, color=factor(etiology.new))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by nyha.binom
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=nyha.binom, color=factor(nyha.binom))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by sah
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=sah, color=factor(sah))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by congenital.anomaly
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=congenital.anomaly, color=factor(congenital.anomaly))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by smoking
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=smoking, color=factor(smoking))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by previos.operation.yn
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=previos.operation.yn, color=factor(previos.operation.yn))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by associated.pathologies.yn
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=associated.pathologies.yn, color=factor(associated.pathologies.yn))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by previous.cv.events.yn
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=previous.cv.events.yn, color=factor(previous.cv.events.yn))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by copd
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=copd, color=factor(copd))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by urgent
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=urgent, color=factor(urgent))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by concomitant.procedures.yn
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=concomitant.procedures.yn, color=factor(concomitant.procedures.yn))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by pulmonary.allograft
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=pulmonary.allograft, color=factor(pulmonary.allograft))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by rvot.augmentation.yn
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=rvot.augmentation.yn, color=factor(rvot.augmentation.yn))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by donor.male
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=donor.male, color=factor(donor.male))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

#Overall pattern by abo.mismatch
ggplot(data = ds, aes(x = fup.duration, y = peak.rvot.grad)) + 
  geom_line(aes(group=id)) + geom_smooth(aes(group=abo.mismatch, color=factor(abo.mismatch))) + 
  scale_x_continuous(breaks = seq(0, 30, by = 1))

