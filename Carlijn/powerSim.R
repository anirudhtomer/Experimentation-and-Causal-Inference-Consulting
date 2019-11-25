install.packages("gtools")
library(gtools)

getPValue = function(sample){
  wilcox.test(sample[,1], sample[,2], paired = T)$p.value
}

nCategories = 6
probabilities = rdirichlet(1, rep(1, nCategories*nCategories))

generateSample = function(sampleSize){
  data = c(rmultinom(1, sampleSize, prob = probabilities))
  
  datamatrix = matrix(nrow = 0, ncol = 2)

  for(i in 1:nCategories){
    for(j in 1:nCategories){
      count = data[(i-1)*nCategories + j]
      if(count > 0){
        datamatrix = rbind(datamatrix, t(replicate(count, c(i,j))))
      }
    }
  }
  
  return(datamatrix)
}

k = 0
nList = seq(10, 400, by =50)
rejectCount = rep(0, length(nList))
for(n in nList){
  k = k + 1
  for(i in 1:1000){
    if(getPValue(generateSample(n))<0.05){
      rejectCount[k] = rejectCount[k] + 1
    }
  }
}

rejectCount/1000

