library(geepack)
#Time is in years

#The total patients for which you want to calculate power
total_patients = 300

#Total time for which you intend to run the trial
total_study_time = 1

#Total clusters
total_clusters = 13

#Lets say the minimum time for which you intend to provide intervention or treatment in any cluster is 1 month
min_time_intervention_or_control = 1/12

#probability of appropriate care
prob_care_control = 0.5
prob_care_intervention = 0.7

#lets say appropriate care probability differs
sd_logodds_care_between_clusters = 0.25
print("Control variation 95%")
print(paste(plogis(log(prob_care_control/(1-prob_care_control)) + c(-2*sd_logodds_care_between_clusters, 2*sd_logodds_care_between_clusters))))

print("Intervention variation 95%")
print(paste(plogis(log(prob_care_intervention/(1-prob_care_intervention)) + c(-2*sd_logodds_care_between_clusters, 2*sd_logodds_care_between_clusters))))

intervention_start_time = seq(min_time_intervention_or_control,
                              total_study_time - min_time_intervention_or_control, length.out = total_clusters)

cluster_log_odds = rnorm(n=total_clusters, mean = 0, sd=sd_logodds_care_between_clusters)
logodds_intervention = log(prob_care_intervention/(1-prob_care_intervention)) + cluster_log_odds
logodds_control = log(prob_care_control/(1-prob_care_control)) + cluster_log_odds
prob_care_control = plogis(logodds_control)
prob_care_intervention = plogis(logodds_intervention)

set.seed(2020)

alpha = 0.05

powermeter = 0
total_sims_to_run = 200
for(i in 1:total_sims_to_run){
  #print(paste("Running simulation number: ", i))
  patient_df = data.frame(cluster = sample(x = 1:total_clusters, size = total_patients, replace = T),
                          arrival_time = runif(n = total_patients, min = 0, max = total_study_time))
  
  patient_df = patient_df[order(patient_df$cluster, patient_df$arrival_time),]
  
  patient_df$cluster_intervention_start_time = intervention_start_time[patient_df$cluster]
  patient_df$treatment = factor(ifelse(test = patient_df$arrival_time>=patient_df$cluster_intervention_start_time,
                                       yes = "Intervention", no="Control"))
  patient_df$treatment = relevel(patient_df$treatment, ref = "Control")
  
  patient_df$prob_care = 0
  patient_df$prob_care[patient_df$treatment=="Intervention"] = prob_care_intervention[patient_df$cluster[patient_df$treatment=="Intervention"]]
  patient_df$prob_care[patient_df$treatment=="Control"] = prob_care_control[patient_df$cluster[patient_df$treatment=="Control"]]
  
  patient_df$appropriate_care = rbinom(n=total_patients, size=1, prob=patient_df$prob_care)
  
  fit_model = geeglm(formula = appropriate_care~treatment, id=cluster, 
                     data=patient_df, family=binomial(link = 'logit'), corstr = 'exchangeable')
  intervention_coef = summary(fit_model)$coefficients[2, c(1,2)]
  
  if(intervention_coef[1] + qnorm(p=alpha/2) *intervention_coef[2] >0){
    powermeter = powermeter + 1
  }
}
print("*************************")
print(paste("Power for", total_patients, "patients is: ", 100*powermeter/total_sims_to_run, "%"))
