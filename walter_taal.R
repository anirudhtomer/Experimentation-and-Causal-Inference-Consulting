nsim = 10000

pa = 0.3
p0 = 0.1

power = c()
for(sample_size in seq(10, 50, 5)){
  res = sapply(1:nsim, FUN = function(i){
    success_count = sum(rbinom(n = sample_size, size = 1, p=pa))
    t=binom.test(x = success_count, n = sample_size, p = p0, 
                 alternative = "greater", conf.level=0.975)$conf.int[1] > p0
  })

  print(paste("Power for", sample_size, "patients:", (table(res)/nsim)["TRUE"]*100, "%"))
}

