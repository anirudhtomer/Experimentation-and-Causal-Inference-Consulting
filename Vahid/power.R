install.packages("ggplot2")
library(ggplot2)

getPValue = function(sample){
    t.test(x=sample[,2],y=sample[,1], paired=T, alternative = "greater")$p.value
}

generateSample = function(sampleSize, correlation=0.9){
  return(MASS::mvrnorm(n=sampleSize, c(1.63,2.52), 
                       Sigma=matrix(c(0.44^2,correlation*0.44*1.05,correlation*0.44*1.05,1.05^2), nrow = 2)))
}

k = 0
nList = seq(5, 100, by=2)
rejectCount = rep(0, length(nList))
for(n in nList){
    k = k + 1
    for(i in 1:1000){
        if(getPValue(generateSample(n))<0.01){  #alpha, type 1 error rate
            rejectCount[k] = rejectCount[k] + 1
        }
    }
}

#Power =0.8, you can change this number and try sample size for other power
sampleSize = nList[which(rejectCount/1000 > 0.99)][1]

qplot(x=c(generateSample(1500)), geom="density", color=rep(c("start", "end"), each=1500))
