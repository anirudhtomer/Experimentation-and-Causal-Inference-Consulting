prevalence = 0.5531

null_hypothesis_sens = 0.9
null_hypothesis_spec = 0.9

alt_hypothesis_sens = 0.8
alt_hypothesis_spec = 0.8

alpha = 0.05
boneferroni_alpha = 0.025

set.seed(2019)
n_simulations = 5000
for(sample_size in seq(50, 500, by = 50)){
  res = rep(0, n_simulations)
  
  for(i in 1:n_simulations){
    true_status = rbinom(n = sample_size, size = 1, prevalence)
    
    predicted_status = true_status
    predicted_status[true_status == 1] = rbinom(n=length(predicted_status[true_status == 1]),size = 1, prob=alt_hypothesis_sens)
    predicted_status[true_status == 0] = rbinom(n=length(predicted_status[true_status == 0]),size = 1, prob=1-alt_hypothesis_spec)
    
    #first test for sensitivity
    sens_success  = sum(predicted_status[true_status==1]==1)
    sens_test = binom.test(x = sens_success, 
                           n = sum(true_status==1), alternative = "less",
                           p=null_hypothesis_sens, conf.level = 1-boneferroni_alpha)
    
    #next test for specificity
    spec_success = sum(predicted_status[true_status==0]==0)
    spec_test = binom.test(x = spec_success, 
                           n = sum(true_status==0), alternative = "less",
                           p=null_hypothesis_spec, conf.level = 1-boneferroni_alpha)
    
    
    if(sens_test$conf.int[2] < null_hypothesis_sens & 
       spec_test$conf.int[2] < null_hypothesis_spec){
      res[i] = 1
    }
  }
  
  print(paste("Power for sample size", sample_size, "is", 100*sum(res)/n_simulations, "%"))
}