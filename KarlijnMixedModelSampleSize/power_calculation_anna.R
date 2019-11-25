library(ggplot2)
library(nlme)
library(MASS)
library(nlme)

#"Placebo" is coded as 1, "GH" is coded as 0
sd_measurement_error = 0.05

#What is the scale of the measurement error for LBM?
sd_lbm = 0.1

#What is the average LBM of all patients at time t=0?
t0lbm_mean = 40
#What is the between patient standard deviation of LBM at t=0?
t0lbm_sd = 8

#What is the (rate of change of LBM) change in LBM per year, when going through placebo first?
PLfirstchangelbm_mean = 0.01
#What is the between patient standard deviation of this rate of change?
PLfirstchangelbm_sd = 0.1

#What is the (rate of change of LBM) change in LBM per year, when going through GH first?
GHfirstchangelbm_mean = 2.26
#What is the between patient standard deviation of this rate of change?
GHfirstchangelbm_sd = 0.755

#What is the (rate of change of LBM) change in LBM per year, when going through placebo second?
PLsecondchangelbm_mean = -0.5
#What is the between patient standard deviation of this rate of change?
PLsecondchangelbm_sd = 0.1

#What is the (rate of change of LBM) change in LBM per year, when going through GH second?
GHsecondchangelbm_mean = GHfirstchangelbm_mean
#What is the between patient standard deviation of this rate of change?
GHsecondchangelbm_sd = GHfirstchangelbm_sd

#for the following questions remember
#positive correlation will mean that the two values in question increase together
#zero correlation means the two values in question have no effect on each other
#negative correlation will mean that the two values in question decrease together

#What is the correlation between a patient's LBM at t=0 and the change in his LBM after going through placebo in year 1?
cor_t0lbm_PLfirstchangelbm = 0
#What is the correlation between a patient's LBM at t=0 and the change in his LBM after going through GH in year 1?
cor_t0lbm_GHfirstchangelbm = 0
#What is the correlation between a patient's LBM at t=0 and the change in his LBM after going through placebo in year 2?
cor_t0lbm_PLsecondchangelbm = 0
#What is the correlation between a patient's change in LBM at t=1 after going thru Placebo in year 1, and the change in his LBM after going through GH in year 1?
#that is when a patient gains weight quickly in year due to GH, is he also likely to have a high weight increase in general even if there was no GH
cor_PLfirstchangelbm_GHfirstchangelbm = 0
#What is the correlation correlation between a patient's change in LBM at t=1 after going thru Placebo in year 1, and the change in his LBM after going through PL in year 2?
cor_PLfirstchangelbm_PLsecondchangelbm = 0
#What is the correlation correlation between a patient's change in LBM at t=1 after going thru GH in year 1, and the change in his LBM after going through PL in year 2?
cor_GHfirstchangelbm_PLsecondchangelbm = 0

D_like_matrix = diag(4)
D_like_matrix[1,1] = t0lbm_sd^2
D_like_matrix[2,2] = PLfirstchangelbm_sd^2
D_like_matrix[3,3] = GHfirstchangelbm_sd^2
D_like_matrix[4,4] = PLsecondchangelbm_sd^2

D_like_matrix[1,2] = cor_t0lbm_PLfirstchangelbm * t0lbm_sd * PLfirstchangelbm_sd
D_like_matrix[2,1] = D_like_matrix[1,2]

D_like_matrix[1,3] = cor_t0lbm_GHfirstchangelbm * t0lbm_sd * GHfirstchangelbm_sd
D_like_matrix[3,1] = D_like_matrix[1,3]

D_like_matrix[1,4] = cor_t0lbm_PLsecondchangelbm * t0lbm_sd * PLsecondchangelbm_sd
D_like_matrix[4,1] = D_like_matrix[1,4]

D_like_matrix[2,3] = cor_PLfirstchangelbm_GHfirstchangelbm * PLfirstchangelbm_sd * GHfirstchangelbm_sd
D_like_matrix[3,2] = D_like_matrix[2,3]

D_like_matrix[2,4] = cor_PLfirstchangelbm_PLsecondchangelbm * PLfirstchangelbm_sd * PLsecondchangelbm_sd
D_like_matrix[4,2] = D_like_matrix[2,4]

D_like_matrix[3,4] = cor_GHfirstchangelbm_PLsecondchangelbm * GHfirstchangelbm_sd * PLsecondchangelbm_sd
D_like_matrix[4,3] = D_like_matrix[3,4]

set.seed(2019)
n = 10
alpha = 0.05
n_sims = 100
total_sims_done = 0

power_calc = rep(0, n_sims)

while(total_sims_done<n_sims){
  data.id = data.frame(pid = 1:n, PLfirst = sample(x = c(0,1), size = n, replace = T))
  data.id$GHfirst = 1-data.id$PLfirst
  
  temp = mvrnorm(n=n, mu = c(t0lbm_mean, PLfirstchangelbm_mean,
                             GHfirstchangelbm_mean, PLsecondchangelbm_mean), Sigma = D_like_matrix)
  
  data.id$t0lbm = temp[, 1]
  data.id$PLfirstchangelbm = temp[, 2]
  data.id$GHfirstchangelbm = temp[, 3]
  data.id$PLsecondchangelbm = temp[, 4]
  data.id$GHsecondchangelbm = data.id$GHfirstchangelbm
  
  data_long = data.id[rep(1:n, each=6),]
  data_long$time = c(0,6,12,15,21,27)/12
  
  #This is not trt given, but the effect of treatment. so for example, effect of GH remains until t=1.25 even if it is stopped at t=0
  data_long$trt = NA
  data_long$trt[data_long$PLfirst==0 & data_long$time<=1.25] = 0
  data_long$trt[data_long$PLfirst==0 & data_long$time>1.25] = 1
  data_long$trt[data_long$PLfirst==1 & data_long$time<=1.25] = 1
  data_long$trt[data_long$PLfirst==1 & data_long$time>1.25] = 0
  
  data_long$lbm = NA
  
  filter1 = data_long$PLfirst == 0 & data_long$time<=1.25
  data_long$lbm[filter1] = data_long$t0lbm[filter1] + data_long$time[filter1] * data_long$GHfirstchangelbm[filter1]
  
  filter2 = data_long$PLfirst == 0 & data_long$time>1.25
  data_long$lbm[filter2] = data_long$t0lbm[filter2] + 1.25 * data_long$GHfirstchangelbm[filter2]
  data_long$lbm[filter2] = data_long$lbm[filter2] + (data_long$time[filter2] - 1.25) * data_long$PLsecondchangelbm[filter2]
  
  filter3 = data_long$PLfirst == 1 & data_long$time<=1.25
  data_long$lbm[filter3] = data_long$t0lbm[filter3] + data_long$time[filter3] * data_long$PLfirstchangelbm[filter3]
  
  filter4 = data_long$PLfirst == 1 & data_long$time>1.25
  data_long$lbm[filter4] = data_long$t0lbm[filter4] + 1.25 * data_long$PLfirstchangelbm[filter4]
  data_long$lbm[filter4] = data_long$lbm[filter4] + (data_long$time[filter4] - 1.25) * data_long$GHsecondchangelbm[filter4]
  
  data_long$lbm = data_long$lbm + rnorm(nrow(data_long), 0, sd_measurement_error)
  
  data_long$tminus1pt25 = (data_long$time - 1.25) * (data_long$time > 1.25)
  
  fitted_model = try(lme(fixed = lbm~time + time:PLfirst + tminus1pt25 + tminus1pt25:PLfirst,
                         random = ~1| pid,
                         data = data_long,
                         control = lmeControl(opt="optim", optimMethod = "L-BFGS-B")), silent = T)
  
  if(!inherits(fitted_model, "try-error")){
    total_sims_done = total_sims_done + 1
    
    confIntFixed = intervals(fitted_model, level = 1-alpha/2, which="fixed")
    lowerConfLimit = confIntFixed$fixed["time:PLfirst","lower"]
    upperConfLimit = confIntFixed$fixed["time:PLfirst","upper"]
    
    first_condition = (lowerConfLimit > 0 & upperConfLimit>0) | (lowerConfLimit<0 & upperConfLimit<0)
    
    lowerConfLimit = confIntFixed$fixed["PLfirst:tminus1pt25","lower"]
    upperConfLimit = confIntFixed$fixed["PLfirst:tminus1pt25","upper"]
    
    second_condition = (lowerConfLimit > 0 & upperConfLimit>0) | (lowerConfLimit<0 & upperConfLimit<0)
    
    if(first_condition==T & second_condition==T){
      power_calc[total_sims_done] = 1
    }
  }
}

print(paste0("Power is ", (sum(power_calc)/n_sims)*100, " %"))
