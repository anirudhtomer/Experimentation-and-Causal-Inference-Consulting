library(ggplot2)
library(nlme)

#lbm is Lean body mass

#The mean lbm for all patients at baseline
beta_0 = 40.2

#The SD of baseline lbm between patients
sigma_b0 = 7.1

#The effect of placebo (kg/year)
beta_time = 0.01

#The effect of treatment excluding effect of placebo (2.26 - 0.01 = 2.25 kg/year)
beta_trt_time = 2.26 - beta_time

min_relevant_beta_trt_time = 0.5 - beta_time

betas = c(beta_0, beta_time, beta_trt_time)

#The measurement error SD
sigma_error = 0.432

getX2 = function(time, firstYear, Trt_first){
  return(time - 1 - firstYear * (time - 1) + 
           Trt_first * (2 + 2 * firstYear * (time - 1) - time))
}

#Model of lbm over time
lme_fixed = obs_lbm ~ 1 + time + x2
fixedFormula = ~ 1 + time + x2
randomFormula = ~ 1

sample_sizes = c(10:50)

two_sided_alpha = 0.05

#The seed allows you to recreate exact results when you run multiple times
set.seed(2019)
nSim = 1000
for(n in sample_sizes){
  
  res = sapply(1:nSim, function(sim_nr){
    patient_id = 1:n
    
    patient_data_list = lapply(patient_id, function(i){
      b0_i = rnorm(n = 1, mean = 0, sd = sigma_b0)  
      Trt_first = sample(c(0,1), size = 1, prob = c(0.5, 0.5))
      
      visitTimes = c(0, 6, runif(1, 12, 13), 18, 24) / 12
      firstYear = ifelse(visitTimes<=1, 1, 0)
      
      df_i = data.frame(id = i, firstYear = firstYear, 
                        time = visitTimes, Trt_first)
      
      df_i$x2 = getX2(df_i$time, df_i$firstYear, df_i$Trt_first)
      
      X = model.matrix(fixedFormula, data=df_i)
      Z = model.matrix(randomFormula, data=df_i)
      
      df_i$true_lbm = X %*% betas + Z %*% b0_i
      
      measurement_error = rnorm(n=length(visitTimes), mean = 0, sd=sigma_error)
      df_i$obs_lbm = df_i$true_lbm + measurement_error
      
      return(df_i)
    })
    
    patient_df = do.call('rbind', patient_data_list)
    
    fitted_model = lme(fixed = lme_fixed, random = ~1|id, data=patient_df)
    confIntFixed = intervals(fitted_model, level = (1-two_sided_alpha/2), which="fixed")
    lowerConfLimit_X2 = confIntFixed$fixed["x2","lower"]
    
    return(lowerConfLimit_X2 > min_relevant_beta_trt_time)
  })
  
  print(paste("Power for", n, "patients:", (table(res)/nSim)["TRUE"]*100, "%"))
}