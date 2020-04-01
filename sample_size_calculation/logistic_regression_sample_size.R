if (!require(contrast)) install.packages('contrast')
library(contrast)

REGION = c("ROT", "AMS", "GRO", "NIJ")

#Time is in months
TIME_UNTIL_FULL_ECMO_EXPERIENCE = 3
ROTTERDAM_ECMO_START_TIME = 3
MAX_STUDY_TIME = 36

HELI_ECMO_START_TIMES = c("ROT"=3, "AMS"=6, "GRO"=9, "NIJ"=12)
HELI_ECMO_FULL_EXPERIENCE_TIMES = HELI_ECMO_START_TIMES + TIME_UNTIL_FULL_ECMO_EXPERIENCE

#Probabilities of having an alive discharge (primary outcome) 
#in various groups. We have 4 groups
#1. There is no ecmo on heli, and no ecmo given ever to patient (15% survival)
#2. There is no ecmo on heli, but ecmo is given at hospital after 65 mins (20% survival)
#3. There is ecmo on heli, but inexperienced team is present. ecmo after 45 mins (30% survival)
#4. There is ecmo on heli, and experienced team is present. ecmo after 35 mins (35% survival)
PR_NEVER_ECMO = 0.15
PR_ECMO_AT_HOSPITAL = 0.20
PR_ECMO_HELI_INEXPERIENCED = 0.3
PR_ECMO_HELI_EXPERIENCED = 0.35

#NULL hypothesis that we want to reject:
#Using ECMO on heli is not useful because the 
#probability of being discharged alive is never more than 2.5% above the scenario with no ECMO

MAX_ECMO_PROB_DIFF = 0.025
NULL_PROB = MAX_ECMO_PROB_DIFF + PR_NO_ECMO
NULL_LOG_ODDS = log(NULL_PROB / (1-NULL_PROB))

#Type I error occurs when NULL hypothesis is true, but it is rejected
TYPE_I_ERROR = 0.025

#Type II error occurs when NULL hypothesis is not rejected, despite it being false

#NR_SIMS decides how accurate our sample size calculations are. 
#Increase until you see no difference due to any further increase
NR_SIMS = 1000
set.seed(2019)
for(NR_PATIENTS_STUDY_PERIOD in c(100, 200, 300, 400, 500)){
  success = sapply(1:NR_SIMS, function(sim_nr){
    #First we generate covariates
    patient_df = data.frame(pid=1:NR_PATIENTS_STUDY_PERIOD)
    
    patient_df$visit_time = sort(runif(n=NR_PATIENTS_STUDY_PERIOD, min = 0, max = MAX_STUDY_TIME), decreasing = F)
    patient_df$region = sample(x=REGION, size = NR_PATIENTS_STUDY_PERIOD, replace = T)
    
    patient_df$ecmo_on_heli = ifelse(patient_df$visit_time >= HELI_ECMO_START_TIMES[patient_df$region], 1, 0)
    patient_df$ecmo_full_experience = ifelse(patient_df$visit_time >= HELI_ECMO_FULL_EXPERIENCE_TIMES[patient_df$heli],1,0)
    #First give heli access to all
    patient_df$access_to_heli = 1
    #remove randomly heli access from control patients
    patient_df$access_to_heli[patient_df$ecmo_on_heli == 0] = sample(c(0,1), size=sum(patient_df$ecmo_on_heli == 0), replace=T)
    
    patient_df$alive_discharge_prob = PR_NEVER_ECMO
    patient_df$alive_discharge_prob[patient_df$access_to_heli==1] = PR_ECMO_AT_HOSPITAL
    
    #Then we generate outcomes
    patient_df$alive_discharge = rbinom(n = NR_PATIENTS_STUDY_PERIOD, size = 1, patient_df$alive_discharge_prob)
    
    #Now we fit a model and do a statistical test
    model = glm(formula = alive_discharge~1 + ecmo, data = patient_df, 
                family = binomial(link = "logit"))
    fitted_log_odds_ecmo = contrast(model, list(ecmo=1), conf.int = 1-TYPE_I_ERROR)
    return(fitted_log_odds_ecmo$Lower > NULL_LOG_ODDS)
  })
  
  print(paste0("Total patients: ", NR_PATIENTS_STUDY_PERIOD, ", and power is ", 100 * mean(success), " %"))
}
