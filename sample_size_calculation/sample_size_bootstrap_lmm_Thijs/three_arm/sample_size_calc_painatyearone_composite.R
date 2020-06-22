##################
# Power calculation for Thijs. 
# Primary outcome: proportion of time in one year in which average pain is below a certain threshold 
##################

library(nlme)
library(MASS)
library(splines)
library(ggplot2)

#Total number of patients (keep it 1000)
n = 1000

#Maximum follow-up time
max_followup_time = 1

#How much proportion of patients will drop out over follow-up
dropout_prop = 0.2

#On average how many years into the study will patients start dropping out.
dropout_start_time = 0.5

#What are the three treatments
treatments = c("Standard", "Botox", "Surgery")

#Leave this function as it is
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

# Per day patient may have headaches many times
# We will do analysis of average pain score per week. 
# This is to make sure of any autocorrelation (patients repeating same data for multiple days) 
# in reported results by patients is taken care of.

#First assumption: pain score of a patient per day and per week is normally distributed
#Following parameters are for the pain score per day.

# Standard deviation of pain score fluctuation when a patient reports it, assuming pain score fluctuates as per normal distribution
sd_ps_fluc = 0.5

#Average pain score just before any treatment, and after one year in: standard treatment, after botox and after surgery 
avg_ps_t0 = 6
avg_ps_t1_std = 5.5
avg_ps_t1_botox = 6
avg_ps_t6weeks_surg = 2.5
avg_ps_t1_surg = 2.0

#This parameter is the time (in years) for which botox effect remains: 3 months = 3/12 years
timeBotoxEffectLasts = 3/12

#This parameter is the average of the maximum pain reduction possible with botox.
maxBotoxPainReduction = 0.5

#Do not change this function. 
botoxCurveOverTime = function(time){
  0.5 * (maxBotoxPainReduction * cos(time * 2*pi/timeBotoxEffectLasts) - maxBotoxPainReduction)
}

#What is the pain threshold below which you want to check the proportion
pain_threshold = 4

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
sd_ps_t6weeks_surg = 0.5
#How much does the pain score vary (give a standard deviation) between patients after surgery?
sd_ps_t1_surg = 1
#How similar will be the pain score before and after the study period if nothing was done
cor_psa_t0_t1 = 0.95

set.seed(2019)

#One sided test with 1% alpha is what we are trying
alpha = 0.05
#dont change powermeter
powermeter = 0

#Total number of patients
n = 50
#Higher the better
total_sims_to_run = 1000
total_remaining_sims = total_sims_to_run
while(total_remaining_sims>0){
  startTime = Sys.time()
  print(paste0("Total simulations remaining: ", total_remaining_sims))
  
  trt = factor(sample(treatments, size = n, replace = T))
  trt = relevel(trt, ref = "Standard")
  max_follow_up = sample(c(1,-1), size = n, prob=c(1-dropout_prop,dropout_prop),replace = T)
  max_follow_up[max_follow_up==-1] = runif(length(max_follow_up[max_follow_up==-1]), dropout_start_time ,1)
  
  patient_df.id = data.frame(id=1:n, trt=trt, max_follow_up=max_follow_up)
  patient_df.id$real_ps_t0 = NA
  patient_df.id$real_ps_t6weeks = NA
  patient_df.id$real_ps_t1 = NA
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
  patient_df.id[patient_df.id$trt == "Surgery", c("real_ps_t0", "real_ps_t6weeks", "real_ps_t1")] = mvrnorm(n = sum(trt=="Surgery"), 
                                                                                                            mu = c(avg_ps_t0, avg_ps_t6weeks_surg, avg_ps_t1_surg),
                                                                                                            Sigma = matrix(c(sd_ps_t0^2, 
                                                                                                                             cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t6weeks_surg,
                                                                                                                             cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_surg,
                                                                                                                             cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t6weeks_surg,
                                                                                                                             sd_ps_t6weeks_surg^2,
                                                                                                                             cor_psa_t0_t1 * sd_ps_t1_surg * sd_ps_t6weeks_surg,
                                                                                                                             cor_psa_t0_t1 * sd_ps_t0 * sd_ps_t1_surg,
                                                                                                                             cor_psa_t0_t1 * sd_ps_t6weeks_surg * sd_ps_t1_surg,
                                                                                                                             sd_ps_t1_surg^2), nrow = 3, ncol=3))
  
  patient_df.id$slope = (patient_df.id$real_ps_t1 - patient_df.id$real_ps_t0)/max_followup_time
  
  patient_df = patient_df.id[rep(1:n, each=365),]
  patient_df$measurement_time = rep(0:364, n) / 365
  patient_df$measurement_time_day = rep(0:364, n)
  patient_df$real_ps = patient_df$real_ps_t0 + patient_df$measurement_time * patient_df$slope
  patient_df$real_ps[patient_df$trt=="Botox"] = patient_df$real_ps[patient_df$trt=="Botox"] + botoxCurveOverTime(patient_df$measurement_time[patient_df$trt=="Botox"])
  
  filter1 = patient_df$trt=="Surgery" & patient_df$measurement_time<=0.125
  patient_df$real_ps[filter1] = patient_df$real_ps_t0[filter1] + patient_df$measurement_time[filter1] * (patient_df$real_ps_t6weeks[filter1] - patient_df$real_ps_t0[filter1])/0.125
  filter2 =  patient_df$trt=="Surgery" & patient_df$measurement_time > 0.125
  patient_df$real_ps[filter2] = patient_df$real_ps_t6weeks[filter2] + patient_df$measurement_time[filter2] * (patient_df$real_ps_t1[filter2] - patient_df$real_ps_t6weeks[filter2])/0.875
  
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
  
  fitted_model = lme(fixed = obs_ps~1 + 
                       ns(measurement_time, knots = c(1.5/12, 3/12, 4.5/12, 6/12, 7.5/12, 9/12, 10.5/12), Boundary.knots = c(0,1)) + 
                       ns(measurement_time, knots = c(1.5/12, 3/12, 4.5/12, 6/12, 7.5/12, 9/12, 10.5/12), Boundary.knots = c(0,1)):trt, 
                     random = list(~1 + measurement_time | id),
                     control = lmeControl(maxIter=100, msMaxIter = 100, niterEM = 50,opt="optim", optimMethod = "L-BFGS-B"),
                     data = patient_df_week, method = "REML")
  
  corrected_alpha = alpha / 4
  if(!inherits(fitted_model, "try-error")){
    
    fixed_effects = fixef(fitted_model)
    contrast = rep(0, length(fixed_effects))
    
    contrast[18:25] = ns(1, knots = c(1.5/12, 3/12, 4.5/12, 6/12, 7.5/12, 9/12, 10.5/12), Boundary.knots = c(0,1))
    meanEff = contrast %*% fixed_effects
    stdErr = sqrt(c(t(contrast) %*% vcov(fitted_model) %*% contrast))
    diff_surgery_standard_upper = qnorm(p=1-corrected_alpha, mean=meanEff, sd=stdErr)
    
    contrast[10:17] = -ns(1, knots = c(1.5/12, 3/12, 4.5/12, 6/12, 7.5/12, 9/12, 10.5/12), Boundary.knots = c(0,1))
    
    meanEff = contrast %*% fixed_effects
    stdErr = sqrt(c(t(contrast) %*% vcov(fitted_model) %*% contrast))
    diff_surgery_botox_upper = qnorm(p=1-corrected_alpha, mean=meanEff, sd=stdErr)
    
    if(diff_surgery_standard_upper < 0  & diff_surgery_botox_upper < 0){
      powermeter = powermeter + 1
    }
    
    total_remaining_sims = total_remaining_sims - 1
    
    print(paste0("Current running estimate of the power is: ", 100 * powermeter / (total_sims_to_run-total_remaining_sims), "%"))
  }
}
print(paste0("Power is: ", 100 * powermeter / total_sims_to_run, "%"))