##################
# Power calculation for Thijs. 
# Primary outcome: proportion of time in one year in which average pain is below a certain threshold 
##################

library(nlme)
library(MASS)
library(splines)
library(ggplot2)

set.seed(2020)

#Total number of patients (keep it 1000)
sample_size = 50

#One sided test with 1% alpha is what we are trying
alpha = 0.05
#dont change powermeter
powermeter = 0

#Higher the better  but also means more time
total_sims_to_run = 1000
#Confidence interval for proportion of time pain is below the threshold is not easy to calculate
#Hence I am using bootstrap. More bootstrap means more computation time, but also means more accurate results
n_bootstrap = 100

total_remaining_sims = total_sims_to_run
while(total_remaining_sims>0){
  startTime = Sys.time()
  print(paste0("Total simulations remaining: ", total_remaining_sims))
  
  patient_df_week = generatePatientData(sample_size)
  
  propSurgery = rep(NA, n_bootstrap)
  propStandard = rep(NA, n_bootstrap)
  
  bs = 1
  while(bs <= n_bootstrap){
    old_ids = sample(1:sample_size, size = sample_size, replace = T)
    new_ids = 1:sample_size
    patient_df_week_bs = do.call('rbind', lapply(1:length(old_ids), function(x){
      gg = patient_df_week[patient_df_week$id == old_ids[x],]
      gg$id = new_ids[x]
      return(gg)
    }))
    
    fitted_model = try(fitModel(patient_df_week_bs), silent = T)
    
    if(!inherits(fitted_model, "try-error")){
      pdata = fitted_model$pdata
      pdata$below_threshold = pdata$pred <= pain_threshold
      propStandard[bs] = mean(pdata$below_threshold[pdata$trt=='Standard'])
      propSurgery[bs] = mean(pdata$below_threshold[pdata$trt=='Surgery'])
      
      bs=bs+1
    }
  }
  endTime=Sys.time()
  print("Simulations round took time:")
  print(endTime - startTime)
  
  #NULL hypothesis is that proportion of time is higher in surgery compared to both botox and standard treatment
  if(quantile(propSurgery - propStandard, probs = c(alpha/2)) > 0){
    powermeter = powermeter + 1
  }
  
  total_remaining_sims = total_remaining_sims - 1
  
  print(paste0("Current running estimate of the power is: ", 100 * powermeter / (total_sims_to_run-total_remaining_sims), "%"))
}

print(paste0("Final estimate of the power is: ", 100 * powermeter / total_sims_to_run, "%"))