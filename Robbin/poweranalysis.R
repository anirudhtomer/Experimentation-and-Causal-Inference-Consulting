getPValue = function(qol_rs, qol_cp){
    wilcox.test(qol_rs, qol_cp, alternative = "g", exact = T)$p.value
}

generateSample = function(sampleSize){
    qol_rs = rnorm(n = sampleSize, mean = 40.3, sd = 16.1)
    qol_cp = rnorm(n = sampleSize, mean = 31.2, sd = 10.4)
    
    qol_rs[qol_rs<18] = 18
    qol_cp[qol_cp<18] = 18
    
    qol_rs[qol_rs>126] = 126
    qol_cp[qol_cp>126] = 126
    
    return(list(qol_rs = qol_rs, qol_cp = qol_cp))
}


k = 0
nList = seq(20, 200, by =10)
rejectCount = rep(0, length(nList))
nSimDs = 1000
for(n in nList){
    k = k + 1
    for(i in 1:nSimDs){
        sample = generateSample(n)
        if(getPValue(sample$qol_rs, sample$qol_cp)<0.05){
            rejectCount[k] = rejectCount[k] + 1
        }
    }
}
names(rejectCount) = paste(nList, "subjects", sep = " ")
power = rejectCount/nSimDs

print(power)
