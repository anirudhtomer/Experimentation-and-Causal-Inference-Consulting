library(nlme)
library(MASS)

#Time is in years
max_followup_time = 1
dropout_prop = 0.2
dropout_start_time = 0.5

treatments = c("Standard", "Botox", "Surgery")

# Per day patient may have headaches many times
# We will do analysis of average pain score per week. 
# This is to make sure of any autocorrelation (patients repeating same data for multiple days) 
# in reported results by patients is taken care of.Although in the sample size calculation 
# I will not generate patient data with autocorrelation

#First assumption: pain score of a patient per day and per week is normally distributed
#Following parameters are for the pain score per day.

# Standard deviation of pain score fluctuation when a patient reports it, assuming pain score fluctuates as per normal distribution
sd_ps_fluc = 0.5

#Average pain score just before any treatment, after standard treatment, after botox and after surgery 
avg_ps_t0 = 6
avg_ps_t1_std = 5.5
avg_ps_t1_botox = 4
avg_ps_t1_surg = 2.5

beta_0 = avg_ps_t0
beta_time_std = (avg_ps_t1_std - avg_ps_t0)/max_followup_time
beta_time_botox = (avg_ps_t1_botox - avg_ps_t0)/max_followup_time - beta_time_std
beta_time_surg = (avg_ps_t1_surg - avg_ps_t0)/max_followup_time - beta_time_std

betas = c(beta_0, beta_time_std, beta_time_botox, beta_time_surg)

#How much does the pain score vary (give a standard deviation) between patients before any treatment?
sd_ps_t0 = 1
#How much does the pain score vary (give a standard deviation) between patients after standard treatment?
sd_ps_t1_std = 1
#How much does the pain score vary (give a standard deviation) between patients after botox treatment?
sd_ps_t1_botox = 1
#How much does the pain score vary (give a standard deviation) between patients after surgery?
sd_ps_t1_surg = 1
#How similar will be the pain score before and after the study period if nothing was done
cor_psa_t0_t1 = 0.95

set.seed(2019)

alpha = 0.01
#dont change powermeter
powermeter = 0

#Total number of patients
n = 50
#Higher the better
total_sims_to_run = 1000
total_remaining_sims = total_sims_to_run
while(total_remaining_sims>0){
  trt = factor(sample(treatments, size = n, replace = T))
  trt = relevel(trt, ref = "Standard")
  max_follow_up = sample(c(1,-1), size = n, prob=c(1-dropout_prop,dropout_prop),replace = T)
  max_follow_up[max_follow_up==-1] = runif(length(max_follow_up[max_follow_up==-1]), dropout_start_time ,1)
  
  patient_df.id = data.frame(id=1:n, trt=trt, max_follow_up=max_follow_up)
  patient_df.id$real_ps_t0 = patient_df.id$real_ps_t1 = NA
  patient_df.id[patient_df.id$trt == "Standard", c("real_ps_t0", "real_ps_t1")] = mvrnorm(n = sum(trt=="Standard"), 
                                                                                          mu = c(avg_ps_t0, avg_ps_t1_std),
                                                                                          Sigma = matrix(c(sd_ps_t0^2, 
                                                                                                           cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_std,
                                                                                                           cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_std,
                                                                                                           sd_ps_t1_std^2), nrow = 2, ncol=2))
  patient_df.id[patient_df.id$trt == "Botox", c("real_ps_t0", "real_ps_t1")] = mvrnorm(n = sum(trt=="Botox"), 
                                                                                       mu = c(avg_ps_t0, avg_ps_t1_botox),
                                                                                       Sigma = matrix(c(sd_ps_t0^2, 
                                                                                                        cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_botox,
                                                                                                        cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_botox,
                                                                                                        sd_ps_t1_botox^2), nrow = 2, ncol=2))
  patient_df.id[patient_df.id$trt == "Surgery", c("real_ps_t0", "real_ps_t1")] = mvrnorm(n = sum(trt=="Surgery"), 
                                                                                         mu = c(avg_ps_t0, avg_ps_t1_surg),
                                                                                         Sigma = matrix(c(sd_ps_t0^2, 
                                                                                                          cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_surg,
                                                                                                          cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_surg,
                                                                                                          sd_ps_t1_surg^2), nrow = 2, ncol=2))
  
  patient_df.id$slope = (patient_df.id$real_ps_t1 - patient_df.id$real_ps_t0)/max_followup_time
  
  patient_df = patient_df.id[rep(1:n, each=365),]
  patient_df$measurement_time = rep(0:364, n) / 365
  patient_df$measurement_time_day = rep(0:364, n)
  patient_df$real_ps = patient_df$real_ps_t0 + patient_df$measurement_time * patient_df$slope
  patient_df$obs_ps = patient_df$real_ps + rnorm(nrow(patient_df), mean = 0, sd = sd_ps_fluc)
  
  patient_df_week = patient_df[(patient_df$measurement_time_day %% 7)==0,]
  
  for(i in 1:nrow(patient_df_week)){
    day = patient_df_week$measurement_time_day[i]
    pid = patient_df_week$id[i]
    if(day > 0){
      day_range = (pid-1)*365 + (((day-7 + 1):day) + 1)
      patient_df_week$obs_ps[i] = mean(patient_df$obs_ps[day_range])
    }
  }
  
  patient_df_week = patient_df_week[patient_df_week$measurement_time <= patient_df_week$max_follow_up,]
  
  patient_df_week = droplevels(patient_df_week[patient_df_week$trt %in% c("Standard", "Surgery"),])
  
  fitted_model = try(expr = lme(fixed = obs_ps~1 + measurement_time + measurement_time:trt, 
                                random = ~1 + measurement_time + measurement_time:trt | id,
                                control = lmeControl(maxIter=100, msMaxIter = 100, niterEM = 50, opt="optim", optimMethod = "L-BFGS-B"),
                                data = patient_df_week, method = "REML"), silent = T)
  if(!inherits(fitted_model, "try-error")){
    if(intervals(fitted_model, which = 'fixed', level = 1-alpha/2)$fixed[3, 'upper'] < 0){
      powermeter = powermeter + 1
    }
    
    total_remaining_sims = total_remaining_sims - 1
  }
}

print(paste0("Power is: ", 100 * powermeter / total_sims_to_run, "%"))