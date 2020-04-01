#In this code we do sample size calculation of contrasts. that is, 
#we check if linear sum of covariates is equal to zero. This is done when we have multiple treatments
#we calculate stderr of the avg estimated contrast to find its confidence interval

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

#Average pain score (ps) just before any treatment, after standard treatment, after botox and after surgery 
avg_ps_t0 = 6
avg_ps_t1_std = 5.5
avg_ps_t1_botox = 4
avg_ps_t1_surg = 3.5

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
powermeter = 0
total_sims_to_run = 10
total_remaining_sims = total_sims_to_run
#Total number of patients
n = 20
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
  patient_df$visit_time = rep(0:364, n) / 365
  patient_df$real_ps = patient_df$real_ps_t0 + patient_df$visit_time * patient_df$slope
  patient_df$obs_ps = patient_df$real_ps + rnorm(nrow(patient_df), mean = 0, sd = sd_ps_fluc)
  
  # patient_df_weekly_ps = patient_df
  # patient_df_weekly_ps$week = ceiling(patient_df$visit_time * 365 / 7)
  # patient_df_weekly_ps$week_alpha = "0"
  # patient_df_weekly_ps$week_alpha[patient_df_weekly_ps$week != 0] = sort(c(letters, LETTERS))[patient_df_weekly_ps$week[patient_df_weekly_ps$week!=0]]
  # patient_df_weekly_ps$week_per_patient = paste(patient_df_weekly_ps$id, patient_df_weekly_ps$week_alpha, sep = "-")
  # 
  # patient_df_weekly = patient_df_weekly_ps[!duplicated(patient_df_weekly_ps$week_per_patient),]
  # patient_df_weekly$avg_ps = sapply(split(patient_df_weekly_ps$obs_ps, 
  #                                         f = patient_df_weekly_ps$week_per_patient), mean)
  # 
  # 
  # patient_df_weekly = patient_df_weekly[patient_df_weekly$week <= ceiling(patient_df_weekly$max_follow_up*365/7), c("id", "trt", "real_ps_t0", "slope", "week", "avg_ps")]
  
  #Converting week to time in years
  #patient_df_weekly$week_in_years = patient_df_weekly$week / 52
  
  fitted_model = try(expr = lme(fixed = obs_ps~1 + visit_time + visit_time:trt, 
                                random = ~1 + visit_time + visit_time:trt | id,
                                control = lmeControl(maxIter=100, msMaxIter = 100, niterEM = 50,opt="optim", optimMethod = "L-BFGS-B"),
                                data = patient_df, method = "REML"), silent = T)
  if(inherits(fitted_model, "try-error")){
    next
  }
  
  #contrast part. we have 3 treatments. so 1 intercept, 1 effect of time, 2 interactions of time and treatment. 
  #no effect of treatment at time 0 and hence no main effect of treatment
  cont1 = c(0.5, 1, 0, 1)
  meanEff = cont1 %*% fixef(fitted_model)
  sdErr = sqrt(c(t(cont1) %*% vcov(fitted_model) %*% cont1))
  #tStat = meanEff / sdErr
  if(qnorm(p=1-alpha, mean = meanEff, sd=sdErr) < 0){
    powermeter = powermeter + 1
  }
  
  total_remaining_sims = total_remaining_sims - 1
}

print(paste0("Power is: ", 100 * powermeter / total_sims_to_run, "%"))