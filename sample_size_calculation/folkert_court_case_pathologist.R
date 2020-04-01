makeTable = function(h0prob=0.2, haprob=0.25, type1error=0.01, 
                     type2error=c(0.01, 0.05, 0.1, 0.2)){
    
    getPValue = function(sample){
        binom.test(x=sample[1], n=sample[2], p=h0prob, alternative="greater")$p.value
    }
    
    generateSample = function(sampleSize){
        return(c(sum(rbinom(n=sampleSize, 1, haprob)), sampleSize))
    }
    
    k = 0
    nList = seq(2, 1000, by=5)
    rejectCount = rep(0, length(nList))
    for(n in nList){
        k = k + 1
        for(i in 1:1000){
            pval = getPValue(generateSample(n))
            #print(pval)
            if(pval<type1error){
                rejectCount[k] = rejectCount[k] + 1
            }
        }
    }
    tt=round(rejectCount/1000,2)
    
    resMatrix = matrix(nrow = 1, ncol = length(type2error))
    resMatrix[1,] = sapply(type2error, FUN = function(err){nList[which(tt>=(1-err))][1]})
    colnames(resMatrix) = paste(type2error*100, "%", sep="")
    rownames(resMatrix) = c("Number of pathologists: ")
    
    cat(paste("Indictor claims that the chance of misdiagnosis was ", h0prob*100, "%. \n", sep=""))
    cat(paste("The real of chance of misdiagnosis is however ", haprob*100, "%. \n \n", sep=""))
    
    cat("The following shows the number of pathologists required, such that a type 1 error of ", 
        type1error*100,"%, and type 2 error as given in the columns of the table is not exceeded. \n \n", sep="")
    print(resMatrix)
    
    cat("\n \n Please note the following: \n")
    cat("Type 1 error is incorrectly rejecting the indictor's claim. ", 
        "That is, incorrectly letting the pathologist go away even though he did a wrong diagnosis. \n",
        sep="")
    cat("Type 2 error is incorrectly accepting indictor's claim. ",
        "That is, incorrectly indicting the pathologist for wrongdoing.", sep="")
}

