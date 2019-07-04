
# This code generates realistic example data for the Transparent Psi Project
library(progress) # to create progress bar

######################################################################
#                                                                    #
#               Set parameters for example data                      #
#                                                                    #
######################################################################

### Set parameters of generating example dataset

# number of erotic trials performed per participant
erotic_trial_size_per_participant = 18

# probability of successful guess if M0 is true
M0_prob = 0.5

# ESP_user_percentage sets the percentage of ESP-users, or ESP-capable individuals 
# in the population 
# set this to 1 to simulate that everyone has the same level of ESP ability
# set this to 0 to simulate that noone has ESP ability
# set this to something in between to simulate that only a portion of the population can
# use ESP, and the others are just guessing randomly, 
# For example ESP_user_percentage = 0.03 and M1_prob = 0.65 simulates that only 3% 
# of the total poupulation can use ESP, and prediction success rate in this subgroup is 65%.
ESP_user_percentage = 1 

# probability of successful guess among ESP-users in the simulated sample
M1_prob = 0.51

# chance of refusing consent or being excluded due to one of the eligibility criteria (e.g. underaged)
chance_for_refuse_consent = 0.01

# chance of stopping prematurely in each trial, to simulate some missing data due to unexpected events
# setting this to 0.002 would generate roughly 2% missing data
# note that premature stopping cannot introduce bias in our confirmatory analyses
chance_for_stopping_session = 0.002



######################################################################
#                                                                    #
#                       Generate example data                        #
#                                                                    #
######################################################################


# simulate the performance of 10,000 potential participants
# the code in this section simulates "realistic" data, producing
# the same data structure as expected from real data collected in the study
# adjust the parameters in the "Set parameters for example data" section

list = list(NA)

# some laboratory info
lab_characteristics_sim = data.frame(exp_ID = paste("exp_", 1:10, sep =""),
                                     exp_ASGS = 1:10,
                                     lab_ID = rep(paste("lab_", 1:5, sep =""), each = 2),
                                     lab_ASGS = 10:1)

# Number of simulated sessions (1 regular session contains 36 trials, 18 erotic and 18 non-erotic)
# generating 10000 simulated participants takes about 4 mins on an i7 6600U 2.6 Ghz CPU 
num_sim = 6000

pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = num_sim, clear = FALSE, width= 60)

for(i in 1:num_sim){
  pb$tick() # adds a tick to the progress bar
  
  ESP_user = rbinom(1, 1, prob = ESP_user_percentage)
  prob = M0_prob + (M1_prob - M0_prob)*ESP_user
  subj_data = as.data.frame(matrix(NA, nrow = 1, ncol = 22))
  names(subj_data) = c("timestamp", "participant_ID", "experimenter_ID_code", "experimenter_ASGS_total_score",
                       "laboratory_ID_code", "sitePI_ASGS_total_score", "session_type", "consent_screen_answer", 
                       "refused_to_answer_sexual_orientation_question","age", "sex" , "final_consent", 
                       "ESP_Q_item_1", "ESP_Q_item_2", "ESP_Q_item_3", "SS_Q_item_1", "SS_Q_item_2", 
                       "trial_number", "guessed_side", "target_side", "reward_type", "sides_match")
  lab_chars = sample(1:nrow(lab_characteristics_sim),1)
  exp_ID = as.character(lab_characteristics_sim[lab_chars,"exp_ID"])
  exp_ASGS = lab_characteristics_sim[lab_chars,"exp_ASGS"]
  lab_ID = as.character(lab_characteristics_sim[lab_chars,"lab_ID"])
  lab_ASGS = lab_characteristics_sim[lab_chars,"lab_ASGS"]
  consent_screen_answer = sample(c("yes", "no"), prob = c(1-chance_for_refuse_consent, chance_for_refuse_consent), 1)
  refused_to_answer_sexual_orientation_question = "no"
  age = sample(c("18-29", "30-44"), 1)
  sex = sample(c("Female", "Male"), 1)
  final_consent = "yes"
  ESP_Q_item_1 = sample(c("Definitely Does not", "Probably does not", "Don't know", "Probably does", "Definitely does"), 1)
  ESP_Q_item_2 = sample(c("Definitely no", "Probably no", "Don't know", "Probably yes", "Definitely yes"), 1)
  ESP_Q_item_3 = sample(c("No, never", "Only a few times", "Occasionally", "Regurarly in the past", "Regurarly no"), 1)
  SS_Q_item_1 = sample(c("very untrue", "untrue", "between true and untrue", "true", "very true"), 1)
  SS_Q_item_2 = sample(c("very untrue", "untrue", "between true and untrue", "true", "very true"), 1)
  reward_type = c(rep(NA, 5), sample(rep(c("erotic", "neutral"), each = erotic_trial_size_per_participant), replace = F))
  
  for(j in 1:((erotic_trial_size_per_participant*2)+5)){
    subj_data[j, "timestamp"] = "simulation"
    subj_data[j, "participant_ID"] = paste("part_", i, sep = "")
    subj_data[j, "experimenter_ID_code"] = exp_ID
    subj_data[j, "experimenter_ASGS_total_score"] = exp_ASGS
    subj_data[j, "laboratory_ID_code"] = lab_ID
    subj_data[j, "sitePI_ASGS_total_score"] = lab_ASGS
    subj_data[j, "session_type"] = "live"
    subj_data[j, "consent_screen_answer"] = if(j == 1){NA} else {consent_screen_answer}
    subj_data[j, "refused_to_answer_sexual_orientation_question"] = if(j == 1 | j == 2){NA} else {refused_to_answer_sexual_orientation_question}
    subj_data[j, "age"] = if(j == 1 | j == 2){NA} else {age}
    subj_data[j, "sex"] = if(j == 1 | j == 2){NA} else {sex}
    subj_data[j, "final_consent"] = if(j == 1 | j == 2 | j == 3){NA} else {final_consent}
    subj_data[j, "ESP_Q_item_1"] = if(j == 1 | j == 2 | j == 3){NA} else {ESP_Q_item_1}
    subj_data[j, "ESP_Q_item_2"] = if(j == 1 | j == 2 | j == 3){NA} else {ESP_Q_item_2}
    subj_data[j, "ESP_Q_item_3"] = if(j == 1 | j == 2 | j == 3){NA} else {ESP_Q_item_3}
    subj_data[j, "SS_Q_item_1"] = if(j == 1 | j == 2 | j == 3 | j == 4){NA} else {SS_Q_item_1}
    subj_data[j, "SS_Q_item_2"] = if(j == 1 | j == 2 | j == 3 | j == 4){NA} else {SS_Q_item_2}
    subj_data[j, "trial_number"] = if(j == 1 | j == 2 | j == 3 | j == 4 | j == 5){NA} else {j-5}
    subj_data[j, "guessed_side"] = if(j == 1 | j == 2 | j == 3 | j == 4 | j == 5){NA} else {sample(c("left", "right"), 1)}
    flip = rbinom(1, 1, prob = prob)
    subj_data[j, "target_side"] = if(j == 1 | j == 2 | j == 3 | j == 4 | j == 5){NA} else { if(flip == 0){ if(subj_data[j, "guessed_side"] == "left"){"right"} else {"left"}} else {if(subj_data[j, "guessed_side"] == "left"){"left"} else {"right"}}}
    subj_data[j, "reward_type"] = if(j == 1 | j == 2 | j == 3 | j == 4 | j == 5){NA} else {reward_type[j]}
    subj_data[j, "sides_match"] = if(j == 1 | j == 2 | j == 3 | j == 4 | j == 5){NA} else {if(subj_data[j, "guessed_side"] == subj_data[j, "target_side"]){"true"} else {"false"}}
    list[[i]] = subj_data
    if(!is.na(subj_data[j, "consent_screen_answer"])){if(subj_data[j, "consent_screen_answer"] == "no") break}
    if(runif(1, min = 0, max = 1) < chance_for_stopping_session) break
  }
}

raw_data = do.call("rbind", list)

