library(nlme)
library(MASS)

set.seed(2020)

#Total number of patients (keep it 1000)
sample_size = 50

#One sided test with 1% alpha is what we are trying
alpha = 0.05
#dont change powermeter
powermeter = 0

#Higher the better  but also means more time
total_sims_to_run = 1000
total_remaining_sims = total_sims_to_run
while(total_remaining_sims>0){
  print(paste0("Total simulations remaining: ", total_remaining_sims))
  
  patient_df_week = generatePatientData(sample_size)
  fitted_model = try(fitModel(patient_df_week), silent = T)
  
  if(!inherits(fitted_model, "try-error")){
    fitted_model = fitted_model$fitted_model
    fixed_effects = fixef(fitted_model)
    contrast = rep(0, length(fixed_effects))
    
    contrast[4:5] = bs(1, knots = 1/12, Boundary.knots = c(0,1), degree=1)
    meanEff = contrast %*% fixed_effects
    stdErr = sqrt(c(t(contrast) %*% vcov(fitted_model) %*% contrast))
    diff_surgery_standard_upper = qnorm(p=1-alpha/2, mean=meanEff, sd=stdErr)
    
    if(diff_surgery_standard_upper < 0){
      powermeter = powermeter + 1
    }
    
    total_remaining_sims = total_remaining_sims - 1
    
    print(paste0("Current running estimate of the power is: ", 100 * powermeter / (total_sims_to_run-total_remaining_sims), "%"))
  }
}

print(paste0("Power is: ", 100 * powermeter / total_sims_to_run, "%"))