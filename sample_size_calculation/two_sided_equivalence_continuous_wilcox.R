#BIC: Bone to Implant Conact

#Actual difference in mean BIC contact between groups
actual_mean_difference = 0

#Actual standard deviation (SD) of BIC in either group, assuming SD is equal in both groups it is 0.157mm
actual_sd_lowspeed = 0.157
actual_sd_highspeed = 0.157

#Equivalance margin is 0.2mm
equi_margin = 0.2

#two-sided test with 5% alpha
alpha = 0.05

power = c()
nsim = 10000
for(sample_size in seq(15, 30, 1)){
  res = sapply(1:nsim, FUN = function(i){
    lowspeed_bic = rnorm(sample_size, mean = 0, sd = actual_sd_lowspeed)
    highspeed_bic = rnorm(sample_size, mean = 0, sd = actual_sd_highspeed)
    
    testres = wilcox.test(x = lowspeed_bic, y = highspeed_bic, alternative = "two.sided",
                conf.level = 1-alpha, conf.int = T, exact = T, paired = F)
    
    return(all(abs(testres$conf.int) < equi_margin))
  })
  
  print(paste("Power for", round((sample_size*2/0.9)), "patients including 10% dropout:", (table(res)/nsim)["TRUE"]*100, "%"))
}


#Project number 2 on volume data
actual_mean_difference = 0

#Actual standard deviation (SD) of volume in either group, assuming SD is equal in both groups
actual_sd_groupA = 2
actual_sd_groupB = 2

#equivalence margin
equi_margin = 3.375

#two-sided test with 5% alpha
alpha = 0.05

power = c()
nsim = 10000
for(sample_size in seq(10, 30, 1)){
  res = sapply(1:nsim, FUN = function(i){
    groupA = rnorm(sample_size, mean = 0, sd = actual_sd_groupA)
    groupB = rnorm(sample_size, mean = 0, sd = actual_sd_groupB)
    
    testres = wilcox.test(x = groupA, y = groupB, alternative = "two.sided",
                          conf.level = 1-alpha, conf.int = T, exact = T, paired = F)
    
    return(all(abs(testres$conf.int) < equi_margin))
  })
  
  print(paste("Power for", round((sample_size*2/0.9)), "patients including 10% dropout:", (table(res)/nsim)["TRUE"]*100, "%"))
}
