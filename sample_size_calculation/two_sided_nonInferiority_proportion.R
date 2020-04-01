nsim = 1000

p_big = 0.16
p_small = 0.11
non_inf = 0.20

p2 = p_small + non_inf

power = c()
for(sample_size in seq(5, 1000, 50)){
  res = sapply(1:nsim, FUN = function(i){
    success_big = sum(rbinom(n = sample_size, size = 1, p=p_big))
    success_2 = sum(rbinom(n = sample_size, size = 1, p=p2))
    
    t = prop.test(x = c(success_2, success_big), n = rep(sample_size,2),
              alternative = "two.sided", conf.level = 0.95, correct = F)$conf.int[1] > 0
  })

  print(paste("Power for", sample_size*2, "patients:", (table(res)/nsim)["TRUE"]*100, "%"))
}

