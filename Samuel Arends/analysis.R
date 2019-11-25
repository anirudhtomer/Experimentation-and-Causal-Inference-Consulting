library(xlsx)

ds = read.xlsx(file="data.xlsx", sheetIndex = 1)
ds = ds[order(ds$ID, ds$Hand),]
ds$rownum = 1:nrow(ds)

#Excluding those with 2 measurements: and one measurement has high noise, or is zero
selectedRows = c(by(data=ds, INDICES = ds$ID, function(x){
  if(nrow(x)==1){
    return(x$rownum)
  }else{
    if(sum(x$No_Sedative)==0 | (x$No_Sedative[1]>0 & x$No_Sedative[2]>0)){
      return(x$rownum[which.min(x$Noise)])
    }else{
      return(x$rownum[which.max(x$No_Sedative)])
    }
  }
}))

ds_filtered = ds[ds$rownum %in% selectedRows,]

setdiff(ds$rownum, ds_filtered$rownum)
#Check if manually you were doing the same
manualExcludeRows = c(2,3,5,8,9,12,17,21,22,25,27,29,32,37)

#Match densities
percentageZero = sum(ds_filtered$No_Sedative==0)/length(ds_filtered$No_Sedative)
realNoSed = ds_filtered$No_Sedative[ds_filtered$No_Sedative > 0]
simNoSed = rnorm(100, mean=mean(log(realNoSed)), sd=sd(log(realNoSed)))
plotDf = data.frame(noSed = c(log(realNoSed), simNoSed), 
                    type=rep(c("real", "sim"), c(length(realNoSed),length(simNoSed))))

ggplot(data=plotDf) + geom_density(aes(x=(noSed), color=type))

getPValue = function(noSed, sed){
  wilcox.test(noSed, sed, alternative = "greater", paired=T, exact = F)$p.value
}

generateSample = function(sampleSize, percDecrease=0.1){
  
  #percentageZero = 0.6
  zeroCount = sum(rbinom(sampleSize,1, percentageZero))
  nonZeroVals = exp(rnorm(sampleSize - zeroCount, mean=mean(log(realNoSed)), sd=sd(log(realNoSed))))
  
  noSed = c(rep(0, zeroCount), nonZeroVals)
  sed = noSed * (1-percDecrease)
  
  return(list(noSed=noSed, sed=sed))
}

#Test if generate sample works
plotDf = data.frame(noSed = c(ds_filtered$No_Sedative, generateSample(100)$noSed), 
                    type=rep(c("real", "sim"), c(length(ds_filtered$No_Sedative),100)))

ggplot(data=plotDf) + geom_density(aes(x=log(noSed+1), color=type))

#Sample size calculations
k = 0
nList = seq(5, 100, by=2)
rejectCount = rep(0, length(nList))
nSimDs = 1000
for(n in nList){
  k = k + 1
  for(i in 1:nSimDs){
    sample = generateSample(n, percDecrease = 0.1)
    if(getPValue(sample$noSed, sample$sed)<0.01){
      rejectCount[k] = rejectCount[k] + 1
    }
  }
}
names(rejectCount) = paste(nList, "subjects", sep = " ")
power = rejectCount/nSimDs

print(power)

