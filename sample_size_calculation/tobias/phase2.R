prevalence = 0.2

gp_sens = 0.4
ai_sens = 0.95

gp_spec = 0.8
ai_spec = 0.8

alpha = 0.05

set.seed(2019)
n_simulations = 5000
for(sample_size in seq(50, 500, by = 50)){
  powermeter = 0
  
  for(i in 1:n_simulations){
    true_status = rbinom(n = sample_size, size = 1, prevalence)
    
    predicted_status_gp = true_status
    predicted_status_ai = true_status
    
    predicted_status_gp[true_status == 1] = rbinom(n=length(predicted_status_gp[true_status == 1]),size = 1, prob=gp_sens)
    predicted_status_gp[true_status == 0] = rbinom(n=length(predicted_status_gp[true_status == 0]),size = 1, prob=1-gp_spec)
    
    predicted_status_ai[true_status == 1] = rbinom(n=length(predicted_status_ai[true_status == 1]),size = 1, prob=ai_sens)
    predicted_status_ai[true_status == 0] = rbinom(n=length(predicted_status_ai[true_status == 0]),size = 1, prob=1-ai_spec)
    
    #first test for sensitivity
    sens_success_gp  = sum(predicted_status_gp[true_status==1]==1)
    sens_success_ai  = sum(predicted_status_ai[true_status==1]==1)
    
    sens_test = prop.test(x = c(sens_success_gp, sens_success_ai), 
                          n = c(sum(true_status), sum(true_status)),  
                          alternative = "less", conf.level = 1-alpha/2)
    
    if(sens_test$conf.int[2] < 0){
      powermeter = powermeter + 1
    }
  }
  
  print(paste("Power for sample size", sample_size, "is", 100*powermeter/n_simulations, "%"))
}