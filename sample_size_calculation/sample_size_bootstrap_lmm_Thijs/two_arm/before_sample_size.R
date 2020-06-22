##################
# Power calculation for Thijs. 
# Primary outcome: proportion of time in one year in which average pain is below a certain threshold 
##################

library(nlme)
library(MASS)
library(splines)
library(ggplot2)
library(ggpubr)

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

#Total number of patients (keep it 1000)
sample_size = 1000

#study period total follow-up time
study_period = 1

#How much proportion of patients will drop out over follow-up
dropout_prop = 0.2

#On average how many years into the study will patients start dropping out.
dropout_start_time = 0.25

#What are the three treatments
treatments = c("Standard", "Surgery")

#What is the pain threshold below which you want to check the proportion
pain_threshold = 4

# Standard deviation of pain score fluctuation when a patient reports it, assuming pain score fluctuates as per normal distribution
sd_ps_fluc = 1

#Average pain score just before any treatment, and after one year in: standard treatment, after botox and after surgery 
avg_ps_t0 = 9
avg_ps_t1_std = 6

avg_ps_t4weeks_surg = 5
avg_ps_t1_surg = 4

#How much does the pain score vary (give a standard deviation) between patients before any treatment?
sd_ps_t0 = 1
#How much does the pain score vary (give a standard deviation) between patients after a full study period of standard treatment?
sd_ps_t1_std = 1
#How similar will be the pain score in standard treatment before and after the study period if nothing was done
cor_ps_t0_t1_std = 0.9

#How much does the pain score vary at 4 weeks (give a standard deviation) between patients after surgery?
sd_ps_t4weeks_surg = 1
#How much does the pain score vary at full study period (give a standard deviation) between patients after surgery?
sd_ps_t1_surg = 1

#How similar will be the pain score in surgery treatment at 0 and 4 weeks
cor_ps_t0_t4weeks_surg = 0.9
#How similar will be the pain score in surgery treatment at 4 weeks and 1 year
cor_ps_t4weeks_t1_surg = 0.7
#How similar will be the pain score in surgery treatment at 0 and at 1 year
cor_ps_t0_t1_surg = 0.4


########################################
# DO NOT change anything below this
#########################################

generatePatientData = function(n){
  trt = factor(sample(treatments, size = n, replace = T))
  max_follow_up = sample(c(1,-1), size = n, prob=c(1-dropout_prop,dropout_prop),replace = T)
  max_follow_up[max_follow_up==-1] = runif(length(max_follow_up[max_follow_up==-1]), dropout_start_time , study_period)
  
  patient_df.id = data.frame(id=1:n, trt=trt, max_follow_up=max_follow_up)
  patient_df.id$real_ps_t0 = NA
  patient_df.id$real_ps_t4weeks = NA
  patient_df.id$real_ps_t1 = NA
  patient_df.id[patient_df.id$trt == "Standard", c("real_ps_t0", "real_ps_t1")] = mvrnorm(n = sum(trt=="Standard"), 
                                                                                          mu = c(avg_ps_t0, avg_ps_t1_std),
                                                                                          Sigma = matrix(c(sd_ps_t0^2, 
                                                                                                           cor_ps_t0_t1_std * sd_ps_t0 * sd_ps_t1_std,
                                                                                                           cor_ps_t0_t1_std * sd_ps_t0 * sd_ps_t1_std,
                                                                                                           sd_ps_t1_std^2), nrow = 2, ncol=2))
  
  patient_df.id[patient_df.id$trt == "Surgery", c("real_ps_t0", "real_ps_t4weeks", "real_ps_t1")] = mvrnorm(n = sum(trt=="Surgery"), 
                                                                                                            mu = c(avg_ps_t0, avg_ps_t4weeks_surg, avg_ps_t1_surg),
                                                                                                            Sigma = matrix(c(sd_ps_t0^2, 
                                                                                                                             cor_ps_t0_t4weeks_surg * sd_ps_t0 * sd_ps_t4weeks_surg,
                                                                                                                             cor_ps_t0_t1_surg * sd_ps_t0 * sd_ps_t1_surg,
                                                                                                                             cor_ps_t0_t4weeks_surg * sd_ps_t0 * sd_ps_t4weeks_surg,
                                                                                                                             sd_ps_t4weeks_surg^2,
                                                                                                                             cor_ps_t4weeks_t1_surg * sd_ps_t1_surg * sd_ps_t4weeks_surg,
                                                                                                                             cor_ps_t0_t1_surg * sd_ps_t0 * sd_ps_t1_surg,
                                                                                                                             cor_ps_t4weeks_t1_surg * sd_ps_t4weeks_surg * sd_ps_t1_surg,
                                                                                                                             sd_ps_t1_surg^2), nrow = 3, ncol=3, byrow = T))
  
  patient_df = patient_df.id[rep(1:n, each=365),]
  patient_df$measurement_time = rep(0:364, n) / 365
  patient_df$measurement_time_day = rep(0:364, n)
  
  patient_df$real_ps = 0
  patient_df$real_ps[patient_df$trt=="Standard"] = patient_df$real_ps_t0[patient_df$trt=="Standard"] + patient_df$measurement_time[patient_df$trt=="Standard"] * (patient_df$real_ps_t1[patient_df$trt=="Standard"]-patient_df$real_ps_t0[patient_df$trt=="Standard"])/study_period
  
  four_weeks = 1/12
  filter1 = patient_df$trt=="Surgery" & patient_df$measurement_time<=four_weeks
  patient_df$real_ps[filter1] = patient_df$real_ps_t0[filter1] + patient_df$measurement_time[filter1] * (patient_df$real_ps_t4weeks[filter1] - patient_df$real_ps_t0[filter1])/four_weeks
  filter2 =  patient_df$trt=="Surgery" & patient_df$measurement_time > four_weeks
  patient_df$real_ps[filter2] = patient_df$real_ps_t4weeks[filter2] + patient_df$measurement_time[filter2] * (patient_df$real_ps_t1[filter2] - patient_df$real_ps_t4weeks[filter2])/(study_period-four_weeks)
  
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
  
  return(patient_df_week)
}

fitModel = function(pat_data){
  fitted_model = lme(fixed = obs_ps~1 + 
                       bs(measurement_time, knots = c(1/12), Boundary.knots = c(0,1), degree = 1) + 
                       bs(measurement_time, knots = c(1/12), Boundary.knots = c(0,1), degree = 1):trt, 
                     random = ~1 + bs(measurement_time, knots = c(1/12), Boundary.knots = c(0,1), degree=1)|id,
                     control = lmeControl(maxIter=100, msMaxIter = 100, niterEM = 50,opt="optim", optimMethod = "L-BFGS-B"),
                     data = pat_data, method = "REML")
  
  total_samples = 1000
  test_times = seq(0, 1, length.out = total_samples)
  pdata = patient_df_week[rep(1, length(test_times)),]
  pdata$measurement_time = test_times
  
  pdata = pdata[rep(1:nrow(pdata), 2),]
  pdata$trt = factor(rep(c("Standard", "Surgery"), each=total_samples))
  
  pdata = effectPlotData(fitted_model, pdata, orig_data = patient_df_week)
  return(list(fitted_model = fitted_model, pdata=pdata))
}

set.seed(2020)
patient_df_week = generatePatientData(n=sample_size)
fitted_model = fitModel(patient_df_week)
pdata = fitted_model$pdata
pdata$below_threshold = pdata$pred <= pain_threshold

#This is how fitted profiles look like
ggplot(pdata) + geom_line(aes(x=measurement_time*12, y=pred, color=trt)) + 
  geom_ribbon(aes(x=measurement_time*12, ymin=low, ymax=upp, group=trt, fill=trt), alpha=0.2)+
  geom_hline(aes(yintercept = pain_threshold, color='Pain Threshold'), 
             size=1, linetype='dashed') +
  theme_bw()+
  scale_y_continuous(breaks = seq(0,10, 1), limits = c(0,10)) + ylab("Pain score") +
  scale_x_continuous(breaks = 0:12, limits = c(0,12))+
  xlab("Time (months)") + theme(legend.position = 'bottom', legend.direction = 'horizontal')


#This is how observed data looks like after taking weekly averages
sampled_patients = sample(1:sample_size,size = 100,replace = F)
ggplot(patient_df_week[patient_df_week$id %in% sampled_patients,]) + 
  geom_line(aes(x=measurement_time*12, y=obs_ps, group=id)) + 
  geom_smooth(aes(x=measurement_time*12, y=obs_ps)) + 
  facet_wrap(~trt)+
  geom_hline(aes(yintercept = pain_threshold, color='Pain Threshold'), 
             size=1, linetype='dashed') +
  theme_bw()+
  scale_y_continuous(breaks = seq(0,10, 1), limits = c(0,10)) + ylab("Obs Pain score") +
  scale_x_continuous(breaks = 0:12, limits = c(0,12))+
  xlab("Time (months)") + theme(legend.position = 'bottom', legend.direction = 'horizontal')


ggplot(patient_df_week[patient_df_week$id %in% sampled_patients,]) + 
  geom_line(aes(x=measurement_time*12, y=real_ps, group=id)) + 
  facet_wrap(~trt)+
  geom_hline(aes(yintercept = pain_threshold, color='Pain Threshold'), 
             size=1, linetype='dashed') +
  theme_bw()+
  scale_y_continuous(breaks = seq(0,10, 1), limits = c(0,10)) + ylab("Real Pain score") +
  scale_x_continuous(breaks = 0:12, limits = c(0,12))+
  xlab("Time (months)") + theme(legend.position = 'bottom', legend.direction = 'horizontal')

print("Proportion of time below thresholds")
by(data=pdata$below_threshold, INDICES = pdata$trt, function(x){paste0(mean(x)*100, "%")})
