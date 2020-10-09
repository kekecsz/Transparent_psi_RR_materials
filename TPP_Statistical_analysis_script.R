
# This script contains the statistical analysis that will be used for the confirmatory
# analysis in the Transparent Psi Project
# Please read the Statistical analysis section to make sense of the code

######################################################################
#                                                                    #
#                           Load packages                            #
#                                                                    #
######################################################################

library(HDInterval) # needed to calcluate HDI credible intervals in the Bayesian parameter estimation robustness test
library(ggplot2) # for plotting
library(emdist) # to calcluate earth mover's distance (EMD)
library(reshape2) # for melt()
library(plyr) # for rounding during plotting
library(lme4) # for glmer

######################################################################
#                                                                    #
#                        Load raw dataset                            #
#                                                                    #
######################################################################

# THIS SECTION WILL CONTAIN THE CODE IMPORTING THE ACTUAL DATA COLLECTED IN THE STUDY.

# to generate example data, use the code on the following link:
# https://github.com/kekecsz/Transparent_psi_RR_materials/blob/master/TPP_Statistical_analysis_script%20-%20generate%20example%20data.R

# or use one of the pregenerated example data files, which were used to produce the figures in the manuscript:

# data in which M0 was simulated to be true
raw_data = read.csv("https://raw.githubusercontent.com/kekecsz/Transparent_psi_RR_materials/master/TPP_example_data_M0.csv")
# data in which M1 was simulated to be true, the file is too big, so it is compressed.
# download the file from https://github.com/kekecsz/Transparent_psi_RR_materials/blob/master/TPP_example_data_M1.zip and open manually


######################################################################
#                                                                    #
#                            Functions                               #
#                                                                    #
######################################################################

######################################################################
#                  Bayes factor calculation functions                #
######################################################################


### Functions for Bayes factor caclulation using beta prior
# These functions are required to run the Bayes factor analysis
# The custom code is necessary because we use beta priors, and 
# the BayesFactor package by default does not have built in beta priors
# We thank Richard Morey for his help in developing these functions!


fullAlt_beta = Vectorize(function(p, y, N, alpha, beta){
  exp(dbinom(y, N, p, log = TRUE) + dbeta(p, alpha, beta, log = TRUE)) 
},"p")

normalize_beta = function(alpha, beta, interval){
  diff(pbeta(interval, alpha, beta))
}

restrictedAlt_beta = function(p,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  fullAlt_beta(p, y, N, alpha, beta) / normalize_beta(alpha, beta, interval) * (p>interval[1] & p<interval[2])
}

margLike_beta = function(y, N, y_prior, N_prior, interval){
  integrate(restrictedAlt_beta, interval[1], interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)[[1]]
}

BF01_beta = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y, N, null_prob) / margLike_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)
},"y")


######################################################################
#                       Other supporting functions                   #
######################################################################

### Function calculating the highest density interval using sampling
# We use hdi() from the library(HDInterval)
# this function is needed for the Bayesian parameter estimation robustness test

mode_HDI <- function(scale, density, crit_width = 0.95, n_samples = 1e5){
  samp <- sample(x = scale, size = n_samples, replace = TRUE, prob = density)
  hdi_result = hdi(samp, credMass=crit_width)
  result = c(scale[which(density == max(density))], # mode
             hdi_result[1], # lower bound
             hdi_result[2]) # upper bound
  
  # only needed for the names of the result
  Crit_lb = (1-crit_width)/2
  Crit_ub = crit_width + (1-crit_width)/2
  
  names(result) = c("mode", paste(Crit_lb*100, "%", sep = ""), paste(Crit_ub*100, "%", sep = ""))
  return(result)
}



# convert logit to probability
# this is used for conversion of the results of the
# logistic regression to the probability scale
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# rule for hypothesis testing inference for Bayesian proportion tests
BF_inference_function = function(BF){
  if(Inference_threshold_BF_low >= BF) {return("M1")
  } else if(Inference_threshold_BF_high <= BF) {return("M0")
  } else {return("Inconclusive")}
}

######################################################################
#                                                                    #
#                          Data analysis                             #
#                                                                    #
######################################################################
### Analysis parameters
# these are the analysis parameters specified in our protocol

# number of erotic trials performed per participant
erotic_trial_size_per_participant = 18

# probability of successful guess if M0 is true
M0_prob = 0.5

# interim analysis points (in total number of erotic trials performed)
when_to_check = c(37836, 62388, 86958, 111528, 136080)

# thresholds to infer support for M0 (high) or M1 (low)
Inference_threshold_BF_high = 25
Inference_threshold_BF_low = 1/Inference_threshold_BF_high

# this information is used both for calculating replication Bayes factor, and the Bayesian parameter estimation robustness test. 
# Here we use data from Bem's experiment 1, 828 successes within 1560 erotic trials
y_prior = 828 #number of successes in erotic trials in Bem's experiment 1
N_prior = 1560 # number of erotic trials in Bem's experiment 1


# smallest effect size of interest in the NHST equivalence tests 
# these are used in both the mixed model analysis in the primary analysis
# and the proportion test in the robustness analysis
minimum_effect_threshold_NHST = 0.01
# p threshold for the NHST tests
# these are used in both the mixed model analysis in the primary analysis
# and the proportion test in the robustness analysis
# although, in the primary analysis this is adjusted for multiple testing
# using Bonferroni's correction
Inference_threshold_NHST = 0.005

# in the Bayesian parameter estimation robustness test this will determine the region of practical 
#equivalence (ROPE) interval. The ROPE is interpreted similarly to SESOI, but not entireli the same. 
# See Kruschke, J. K., & Liddell, T. M. (2017). The Bayesian New Statistics: Hypothesis testing, 
# estimation, meta-analysis, and power analysis from a Bayesian perspective. 
# Psychonomic Bulletin & Review, 1-29. 
minimum_effect_threshold_Bayes_Par_Est = 0.006

# This threshold is used to set the HDI width to check against the ROPE in the Bayesian parameter 
# estimation robustness test, if ths parameter is set to 0.05, it means that we would 
# expect that 95% of the probability mass would be within the ROPE to support M0, or outside of it
# to support M1
Inference_threshold_robustness_Bayes_Par_Est = 0.05 


######################################################################
#                                                                    #
#                          Data management                           #
#                                                                    #
######################################################################

raw_data[,"sides_match"] = as.factor(tolower(as.logical(raw_data[,"sides_match"])))
raw_data[,"participant_ID"] = as.factor(raw_data[,"participant_ID"])

# sessions conducted with the test accounts or without lab_IDs are excluded
lab_IDs_to_exclude <- c("", "18155ef201564afbb81f6a8b74aa9a033eac51ec6595510eca9606938ffaced3", "ece83ceb8611d1926746e5bb3597ed1e8cb5d336521331b31961d5c0348883cf", "bd2dd15be34863e9efb77fbddfe744382a9c62c6a497e8bcf3097a47905b905b", "fff9cb9dcc3ac735fc25a59f424e98278a731c23ccd57276d292996c2ba7784f")
data_nontest <- raw_data[!(raw_data[,"laboratory_ID_code"] %in% lab_IDs_to_exclude), ]

# add a row_counter, which will be useful to distinguish data coming in after the stopping rule was met.
data_nontest[, "row_counter"] = 1:nrow(data_nontest)

# extract data from erotic trials 
data_nontest_trials = data_nontest[!is.na(data_nontest[, "trial_number"]),]
data_nontest_trials_erotic = data_nontest_trials[data_nontest_trials[, "reward_type"] == "erotic", ]
# drop unused factor levels
data_nontest_trials_erotic[,"participant_ID"] = droplevels(data_nontest_trials_erotic[,"participant_ID"])



######################################################################
#                                                                    #
#                    Primary confirmatory test                       #
#                                                                    #
######################################################################

# This section conducts the primary confirmatory analysis at each stopping point.
# It also cuts the data at the point where one of the stopping rules has been met.

results_table = data.frame(matrix(NA, nrow = 1, ncol = 7))
names(results_table) = c("Mixed_mod_CIlb", "Mixed_mod_CIub", "mixed_CI_width","BF_replication", "BF_uniform", "BF_BUJ", "checked_at")

# this is a counter to count the number of tests conducted using the mixed model
# due to sequential testing. This is used to adjust the p-value threshold 
# for the number of comparions made
comparisons_Mixed_NHST = 0

for(i in 1:length(when_to_check)){

  # determin current stopping point and next stopping point
  current_stopping_point = when_to_check[i]
  if(i < length(when_to_check)){next_stopping_point = when_to_check[i+1]} else {next_stopping_point = "last"}
  print(paste("analyzing at reaching", current_stopping_point, "erotic trials"))
  
  # sampling starting from the beggining of the full simulated dataset (from the first trial of the first participant) 
  # until reaching the next interim analysis point
  data_BF = data_nontest_trials_erotic[1:current_stopping_point,]
  last_row = data_BF[nrow(data_BF), "row_counter"]
  # number of successes and total N of trials
  successes = sum(as.logical(data_BF[,"sides_match"]))
  data_BF[,"sides_match_numeric"] = as.numeric(as.logical(data_BF[,"sides_match"]))
  total_N = current_stopping_point
  results_table[i, "checked_at"] = current_stopping_point
  
  #================================================================#
  #            Mixed effect logistic regression analysis           #
  #================================================================#
  
  # advance the counter to see how much adjustment needs to be made to the
  # NHST inference threshold due to multiple testing
  comparisons_Mixed_NHST = comparisons_Mixed_NHST + 2 # we add 2 at each sequential stopping point because we do two tests at each stop point, one for M0 and one for M1
  
  # build mixed logistic regression model and extract model coefficient and SE
  mod_mixed = glmer(sides_match_numeric ~ 1 + (1|participant_ID), data = data_BF, family = "binomial")
  estimate_mixed = summary(mod_mixed)$coefficients[1,1]
  se_mixed = summary(mod_mixed)$coefficients[1,2]
  
  # compute confidence interval on the probability scale, and save into results_table
  results_table[i,"mixed_CI_width"] = 1-(Inference_threshold_NHST/comparisons_Mixed_NHST)
  wald_ci_mixed_logit <- c(estimate_mixed - se_mixed* qnorm(1-((Inference_threshold_NHST/comparisons_Mixed_NHST)/2)),
                           estimate_mixed + se_mixed* qnorm(1-((Inference_threshold_NHST/comparisons_Mixed_NHST)/2)))
  wald_ci_mixed = logit2prob(wald_ci_mixed_logit)
  
  results_table[i, "Mixed_mod_CIlb"] = wald_ci_mixed[1]
  results_table[i, "Mixed_mod_CIub"] = wald_ci_mixed[2]
  

# Statistical inference based on the results of the mixed model analysis  
  
  minimum_effect = M0_prob+minimum_effect_threshold_NHST
  if(results_table[i, "Mixed_mod_CIub"] < minimum_effect){Mixed_NHST_inference = "M0"
    } else if(results_table[i, "Mixed_mod_CIlb"] > M0_prob){Mixed_NHST_inference = "M1"
    } else {Mixed_NHST_inference = "Inconclusive"}

  
  
  
  
  
  #================================================================#
  #        Calculating Bayes factors using different priors        #
  #================================================================#
  
  # as determined in the analysis plan, three different prior distributions are used for M1
  # to ensure the robustness of the statistical inference to different analytical choices
  # the same 
  
  ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
  
  BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_replication"] = round(BF_replication, 3)
  BF_replication_inference = BF_inference_function(BF_replication)
  
  
  ### Bayes factor with uniform prior
  # using a non-informative flat prior distribution with alpha = 1 and beta = 1
  
  BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_uniform"] = round(BF_uniform, 3)
  BF_uniform_inference = BF_inference_function(BF_uniform)
  
  ### Bayes factor with BUJ prior
  # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
  # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
  # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
  # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
  # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
  # This distribution has 90% of its probability mass under p = 0.712, which we determined 
  # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
  # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
  # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
  # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
  # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
  
  BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  results_table[i, "BF_BUJ"] = round(BF_BUJ, 3)
  BF_BUJ_inference = BF_inference_function(BF_BUJ)
  

  #================================================================#
  #                    Main analysis inference                     #
  #================================================================#
  
    # determine final inference (supported model) based on the inferences drawn
    # from the mixed model and the Bayes factors 
  if(all(c(Mixed_NHST_inference, BF_replication_inference, BF_uniform_inference, BF_BUJ_inference) == "M1")) {
    primary_analysis_inference = "M1"
    which_threshold_passed = as.character(Inference_threshold_BF_low)
    print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
    break} else if(all(c(Mixed_NHST_inference, BF_replication_inference, BF_uniform_inference, BF_BUJ_inference) == "M0")) {
      primary_analysis_inference = "M0"
      which_threshold_passed = as.character(Inference_threshold_BF_high)
      print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
      break} else if((next_stopping_point != "last") & (nrow(data_nontest_trials_erotic) < next_stopping_point)){
        primary_analysis_inference = "Ongoing"
        which_threshold_passed = "Ongoing"
        print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))
        break} else {
          primary_analysis_inference = "Inconclusive"
          which_threshold_passed = paste("neither ", Inference_threshold_BF_low, " or ", Inference_threshold_BF_high, sep = "")
          print(paste("main analysis inference at latest stopping point =", primary_analysis_inference))}
  
}

##################################################
#   Results: Sample and study characteristics    #
##################################################

# The code in this segment calculates all sample and study characteristics 
# mentioned in the results sescion of the manuscript

data_nontest_untilstudystop = data_nontest[1:which(data_nontest[, "row_counter"] == last_row),]
data_nontest_trials_erotic_untilstudystop = data_nontest_trials_erotic[1:which(data_nontest_trials_erotic[, "row_counter"] == last_row),]

data_nontest[,"participant_ID"] = droplevels(data_nontest[,"participant_ID"])
N_participants_started_session_total = length(unique(data_nontest[, "participant_ID"]))

data_nontest_untilstudystop[,"participant_ID"] = droplevels(data_nontest_untilstudystop[,"participant_ID"])
N_participants_started_session_untilstudystop = length(unique(data_nontest_untilstudystop[, "participant_ID"]))

N_participants_started_session_afterstudystop = N_participants_started_session_total - N_participants_started_session_untilstudystop

data_nontest_trials_erotic_untilstudystop[,"participant_ID"] = droplevels(data_nontest_trials_erotic_untilstudystop[,"participant_ID"])
N_participants_started_trials_untilstudystop = length(unique(data_nontest_trials_erotic_untilstudystop[, "participant_ID"]))

data_BF[,"participant_ID"] = droplevels(data_BF[,"participant_ID"])
N_participants_data_included_in_main_analysis = length(unique(data_BF[, "participant_ID"]))

proportion_participants_novaliddata_untilstudystop = (N_participants_started_session_untilstudystop - N_participants_data_included_in_main_analysis)/N_participants_started_session_untilstudystop

N_erotic_trials = total_N



data_BF_split_pre = data_BF
data_BF_split_pre[,"sides_match"] = as.character(data_BF_split_pre[,"sides_match"])
data_BF_split = split(data_BF_split_pre, f = data_BF_split_pre[,"participant_ID"])
first_rows_of_each_participant = do.call(rbind,lapply(data_BF_split, function(x) x[1,]))
 
age_range_of_most_participants = names(which.max(table(first_rows_of_each_participant[,"age"])))
age_range_of_most_participants_proportion = round(table(first_rows_of_each_participant[,"age"])[which.max(table(first_rows_of_each_participant[,"age"]))]/sum(table(first_rows_of_each_participant[,"age"])), 4)

N_sex_women = sum(first_rows_of_each_participant[,"sex"] == "Female")
sex_women_proportion = N_sex_women/N_participants_data_included_in_main_analysis
N_sex_men = sum(first_rows_of_each_participant[,"sex"] == "Male")
sex_men_proportion = N_sex_men/N_participants_data_included_in_main_analysis

first_rows_of_each_participant[,"ESP_Q_item_1_num"] = as.numeric(ordered(factor(first_rows_of_each_participant[,"ESP_Q_item_1"]), levels = c("Definitely Does not", "Probably does not", "Don't know", "Probably does", "Definitely does")))
ESP_Q_mean = mean(first_rows_of_each_participant[,"ESP_Q_item_1_num"])
ESP_Q_SD = sd(first_rows_of_each_participant[,"ESP_Q_item_1_num"])

first_rows_of_each_participant[,"SS_Q_item_1_num"] = as.numeric(ordered(factor(first_rows_of_each_participant[,"SS_Q_item_1"]), levels = c("very untrue", "untrue", "between true and untrue", "true", "very true")))
# reverst scale as this question is reverse scored 
first_rows_of_each_participant[,"SS_Q_item_1_num"] = 6 - first_rows_of_each_participant[,"SS_Q_item_1_num"]
first_rows_of_each_participant[,"SS_Q_item_2_num"] = as.numeric(ordered(factor(first_rows_of_each_participant[,"SS_Q_item_2"]), levels = c("very untrue", "untrue", "between true and untrue", "true", "very true")))
first_rows_of_each_participant[,"SS_Q_average_score"] = apply(cbind(first_rows_of_each_participant[,"SS_Q_item_1_num"], first_rows_of_each_participant[,"SS_Q_item_2_num"]), 1, mean)
SS_Q_mean = mean(first_rows_of_each_participant[,"SS_Q_average_score"])
SS_Q_SD = sd(first_rows_of_each_participant[,"SS_Q_average_score"])

N_guessed_side_left = sum(data_BF[,"guessed_side"] == "left")
guessed_side_left_proportion = N_guessed_side_left/length(data_BF[,"guessed_side"] == "left" | data_BF[,"guessed_side"] == "right")

N_target_side_left = sum(data_BF[,"target_side"] == "left")
target_side_left_proportion = N_target_side_left/length(data_BF[,"target_side"] == "left" | data_BF[,"target_side"] == "right")

data_nontest_trials[, "participant_ID"] = droplevels(data_nontest_trials[, "participant_ID"])
data_nontest_trials_analyzedparticipants = data_nontest_trials[data_nontest_trials[, "participant_ID"] %in% unique(data_BF[, "participant_ID"]),]

data_nontest_trials_analyzedparticipants[, "participant_ID"] = droplevels(data_nontest_trials_analyzedparticipants[, "participant_ID"])
data_nontest_trials_analyzedparticipants_split = split(data_nontest_trials_analyzedparticipants, f = data_nontest_trials_analyzedparticipants[,"participant_ID"])


N_erotic_trials_per_participant = do.call(rbind,lapply(data_nontest_trials_analyzedparticipants_split, function(x) sum(x[,"reward_type"] == "erotic")))
N_sessions_terminated = nrow(N_erotic_trials_per_participant) - sum(N_erotic_trials_per_participant == 18)
sessions_terminated_proportion = N_sessions_terminated/nrow(N_erotic_trials_per_participant)
N_missing_erotic_trials = sum(18 - N_erotic_trials_per_participant)



##################################################
#            Results: Primary analysis           #
##################################################

# proportion of successful guesses
success_proportion = successes/total_N
success_proportion
# total number of erotic trials
N_erotic_trials
# statistical inference
primary_analysis_inference
# which threshold was passed
which_threshold_passed

## Figure 1
# The Bayes factors calculated at each new experimental trial are displayed on this figure. 

# final Bayes factor results to set the scale of the plot
fig_1_plotdata = data.frame(cbind(c("BF-replication", "BF-uniform", "BF-BUJ"), round(c(BF_replication, BF_uniform, BF_BUJ), 3)))
names(fig_1_plotdata) = c("BF_type", "BF_value")
fig_1_plotdata[,"BF_value"] = as.numeric(as.character(fig_1_plotdata[,"BF_value"]))

# setting limits and breaks and label positions
fig_1_y_axis_breaks = if(min(fig_1_plotdata[,"BF_value"]) < 0.01 | max(fig_1_plotdata[,"BF_value"]) > 100){
  c(0.002, 0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100, 500)
} else {c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100)}
fig_1_y_axis_limits = if(min(fig_1_plotdata[,"BF_value"]) < 0.01 | max(fig_1_plotdata[,"BF_value"]) > 100){c(0.001, 1000)} else {c(0.005, 200)}
fig_1_y_axis_text_position = if(min(fig_1_plotdata[,"BF_value"]) < 0.01 | max(fig_1_plotdata[,"BF_value"]) > 100){c(200, 1, 0.005)} else {c(100, 1, 0.01)}

# calculating cumulative successes
cumulative_successes_trial_N_table = data.frame(cbind(cumsum(as.logical(data_BF[,"sides_match"])), 1:nrow(data_BF)))
names(cumulative_successes_trial_N_table) = c("successes", "total_N")

# calculating Bayes factors for each new experimental trial with all three priors
# calcluating these takes about 3 minutes on an i7-6600 2.6GHz GPU, dependening on how many data points are there
BF_replication_cumulative <- apply(cumulative_successes_trial_N_table, 1, function(x) BF01_beta(y = x["successes"], N = x["total_N"], y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob)) 
BF_uniform_cumulative <- apply(cumulative_successes_trial_N_table, 1, function(x) BF01_beta(y = x["successes"], N = x["total_N"], y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob)) 
BF_BUJ_cumulative <- apply(cumulative_successes_trial_N_table, 1, function(x) BF01_beta(y = x["successes"], N = x["total_N"], y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob))

# fitting smoothing spline to get a smooth curve
fit_replication = smooth.spline(1:nrow(data_BF), BF_replication_cumulative, df = 80)
fit_uniform = smooth.spline(1:nrow(data_BF), BF_uniform_cumulative, df = 80)
fit_BUJ = smooth.spline(1:nrow(data_BF), BF_BUJ_cumulative, df = 80)

BF_replication_cumulative_spline = predict(fit_replication)$y
BF_uniform_cumulative_spline = predict(fit_uniform)$y
BF_BUJ_cumulative_spline = predict(fit_BUJ)$y

# put all data needed for plotting in a data frame
fig_1_plotdata_full = data.frame(cbind(c(BF_replication_cumulative_spline, BF_uniform_cumulative_spline, BF_BUJ_cumulative_spline), rep(1:nrow(data_BF), 3), rep(c("BF_replication", "BF_uniform", "BF_BUJ"), each = nrow(data_BF))))
names(fig_1_plotdata_full) = c("BF_value", "total_N", "BF_type")
fig_1_plotdata_full[,"total_N"] = as.numeric(as.character(fig_1_plotdata_full[,"total_N"]))
fig_1_plotdata_full[,"BF_value"] = as.numeric(as.character(fig_1_plotdata_full[,"BF_value"]))

# Bayes factors at each predetermined sequential stopping point
BF_table = results_table[,c("BF_replication", "BF_uniform", "BF_BUJ", "checked_at")]
BF_table_melted = melt(BF_table[nrow(BF_table),], id=c("checked_at"))

# generate figure 1
figure_1 <- ggplot(fig_1_plotdata_full, aes(y = BF_value, x = total_N, group = BF_type))+
  geom_line(aes(linetype = BF_type), size = 1.2)+
  scale_linetype_manual(name="Prior", labels=c("BUJ","Replication","Uniform"), values=c("solid", "dashed", "twodash"))+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf), alpha = 0.4, fill=c("grey60"))+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high), alpha = 0.2, fill=c("grey80"))+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low), alpha = 0.4, fill=c("grey60"))+
  scale_y_log10(limits = fig_1_y_axis_limits, breaks = fig_1_y_axis_breaks)+
  scale_x_continuous(breaks = c(0, round_any(BF_table[1,"checked_at"]*1/3, 1000), round_any(BF_table[1,"checked_at"]*2/3, 1000), BF_table[,"checked_at"]))+
  geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
  geom_vline(xintercept = BF_table[,"checked_at"], linetype = "dotted")+
  geom_point(data = BF_table_melted, aes(y = value, x = checked_at, group = variable), size = 3.5, shape = 21, fill = "white")+
  annotate("text", x=-500, y=c(fig_1_y_axis_text_position), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270)+
  xlab("Number of experimental trials")+
  ylab("Bayes factor")+
  theme_bw()+
  theme(legend.position="bottom")
  


# ggplot returns this warning:
# "Transformation introduced infinite values in continuous y-axis"
# but this is not important so we suppress it
suppressWarnings(print(figure_1))




##################################################
#                 Robustness tests               #
##################################################

#================================================================#
#               Robustness test of BF results with NHST          #
#================================================================#

# robustness of BF results is tested with NHST proportion tests
# here we perform both an equivalence test and an equality test to draw statistical inference

# equivalence test
equivalence_test_p = prop.test(x = successes, n = total_N,p = M0_prob+minimum_effect_threshold_NHST, alternative = "less")$p.value

# equality test
equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value

# making statistical inference
if(Inference_threshold_NHST > equivalence_test_p){inference_robustness_NHST = "M0"
  } else if(Inference_threshold_NHST > equality_test_p){inference_robustness_NHST = "M1" 
  } else {inference_robustness_NHST = "Inconclusive"}


#=======================================================================#
#   Robustness test of BF results with Bayesian parameter estimation    #
#=======================================================================#

# robustness of BF results is tested by calculating HDI of the posterior distribution and checking its relation to
# the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
# (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
# analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 

# calculate posterior distribution using beta distribution updating
prior_alpha = y_prior + 1
prior_beta = N_prior-y_prior+1

posterior_alpha = prior_alpha + successes
posterior_beta = prior_beta + total_N - successes

scale = seq(0, 1, length = 10001)
posterior_density = dbeta(scale, posterior_alpha, posterior_beta)

# calculate HDI for the posterior distribution
# (here we calculate the upper and lower bound of the 90% of the probability mass
# because we use a one-tailed test. This means that the total probability mass below
# the upper bound of the 90% HDI will be 95%)
hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)

# parameters for decision making
HDI_lb = hdi_result[2]
HDI_ub = hdi_result[3]

# making inference decision
ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
} else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
} else {inference_robustness_Bayes_Par_Est = "Inconclusive"}

# probability that the parameter falls outside of the ROPE
Probability_parameter_higher_than_ROPE = sum(posterior_density[scale>ROPE])/sum(posterior_density)

## Figure 2
# Figure 2 displays the Confidence interval computed based on the 
# final mixed model primary analysis and the posterior distribution of the 
# parameter based on the Bayesian Parameter Estimation

fig_2_sample <- as.data.frame(sample(x = scale, size = 1000000, replace = TRUE, prob = posterior_density))
names(fig_2_sample) = "value"
fig_2_sample_within_HDI = as.data.frame(fig_2_sample[fig_2_sample > (hdi_result[2]-0.0001) & fig_2_sample < (hdi_result[3]+0.0001),])
names(fig_2_sample_within_HDI) = "value"

figure_2_pre = ggplot(fig_2_sample, aes(x = value))+
  geom_density(aes(y = ..scaled..), color = "white")

fig_2_plotdata <- ggplot_build(figure_2_pre)$data[[1]]
fig_2_segment_data <- data.frame(x1 = c(hdi_result[2], hdi_result[3]), y1 = c(0, 0), xend=c(hdi_result[2], hdi_result[3]),
                                 yend=c(approx(x = fig_2_plotdata$x, y = fig_2_plotdata$y, xout = hdi_result[2])$y, approx(x = fig_2_plotdata$x, y = fig_2_plotdata$y, xout = hdi_result[3])$y), 
                                 lty = c("solid"))

mixed_NHST_final_CI = as.numeric(results_table[nrow(results_table),c("Mixed_mod_CIlb", "Mixed_mod_CIub")])
fig_2_propCI_segment_data <- data.frame(x1 = c(mixed_NHST_final_CI[1], mixed_NHST_final_CI[1], mixed_NHST_final_CI[2]), y1 = c(-0.04, -0.02, -0.02), xend=c(mixed_NHST_final_CI[2], mixed_NHST_final_CI[1], mixed_NHST_final_CI[2]),
                                 yend=c(-0.04, -0.06, -0.06), lty = c("solid"))

fig_2_ROPE_equlim_segment_data <- data.frame(x1 = c(ROPE, M0_prob+minimum_effect_threshold_NHST), y1 = c(0, 0), xend=c(ROPE, M0_prob+minimum_effect_threshold_NHST),
                                        yend=c(Inf, -Inf))


figure_2 = figure_2_pre + 
            geom_area(data = subset(fig_2_plotdata, x >= hdi_result[2] & x <= hdi_result[3]), aes(x=x, y=y), fill="light gray")+
            geom_line(data = fig_2_plotdata, aes(x=x, y=y))+
  geom_hline(yintercept = 0)+
            geom_segment(data = fig_2_segment_data, aes(x=x1, xend=xend, y=y1, yend=yend))+
            geom_segment(data = fig_2_propCI_segment_data, aes(x=x1, xend=xend, y=y1, yend=yend))+
            geom_segment(data = fig_2_ROPE_equlim_segment_data, aes(x=x1, xend=xend, y=y1, yend=yend), linetype=c("dashed", "dotted"))+
            annotate("text", x = hdi_result[1], y = 0.1, label = paste("90% HDI: [", round(hdi_result[2], 3), ", ", round(hdi_result[3], 3), "]", sep = ""), size = 4, fontface = 2) +
            annotate("text", x = mean(c(mixed_NHST_final_CI[1], mixed_NHST_final_CI[2])), y = -0.08, label = paste(results_table[nrow(results_table),"mixed_CI_width"] ," CI: [", round(mixed_NHST_final_CI[1], 3), ", ", round(mixed_NHST_final_CI[2], 3), "]", sep = ""), size = 4, fontface = 2) +
            geom_vline(xintercept = M0_prob, size = 1.3)+
            xlab("successful guess probability")+
            ylab("scaled density")+
            ylim(-0.1, 1)+
            theme_bw()+
            theme(panel.border = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
          #        legend.position = "bottom",
                  axis.line = element_line(colour = "black"),
                  axis.text.x = element_text(color = "black", face = "bold", size = 12, margin = margin(1,0,0,0,"mm"), vjust = 0.5),
                  axis.text.y = element_text(face = "bold", color = "black", size = 12),
                  axis.title = element_text(size = 16))
figure_2

#=======================================================================#
#      Determine final inference of all robustness tests combined       #
#=======================================================================#

# the main analysis inference is only robust if all robustness tests came to the same 
# inference as the final inference of the primary analysis
inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}
Robust = if(primary_analysis_inference == "Inconclusive"){"NA, main inference inconclusive"
} else if(primary_analysis_inference == "Ongoing"){"NA, main analysis ongoing"
} else if(primary_analysis_inference == inference_robustness){"robust"
} else {"not robust"}
Robust



######################################################################
#                         Exploratory analysis                       #
######################################################################
# EXPLORATORY ANALYSIS RESULTS WILL NOT AFFECT THE CONCLUSIONS OF OUR STUDY

# This exploratory analysis aims to explore the empirical distribution of
# successful guess rate in the population, by contrasting it to a theoretical
# distribution expected if the true successful guess chance is homogeneously
# 50% in the population. This may provide information for future research on
# potential irregulities in the distribution of guess rates.


#=======================================================================#
#           Comparison of expected and observed distributions           #
#=======================================================================#

# calculate proportion of successful guesses for each participant in the observed data
data_BF_split_finishedalltrials = data_BF_split[which(sapply(data_BF_split, nrow) == 18)]
success_proportions_empirical_finishedalltrials = sapply(data_BF_split_finishedalltrials, function(x) mean(as.logical(x[,"sides_match"])))

# samples 1,000,000 participants from a population with a 50% successfull guess chance
# homogeneous in the population
# we call this the theoretical sample, because it approximates the theoretical null model
sim_null_participant_num = 1000000
success_proportions_theoretical <- rbinom(sim_null_participant_num, size = erotic_trial_size_per_participant, prob=M0_prob)/erotic_trial_size_per_participant

# determine possible values of success rates
possible_success_rates = 0
for(i in 1:erotic_trial_size_per_participant){
  possible_success_rates[i+1] = round(1/(erotic_trial_size_per_participant/i), 2)
}
possible_success_rates_char = as.character(possible_success_rates)
success_proportions_theoretical_char_rounded = as.character(round(success_proportions_theoretical, 2))
success_proportions_empirical_finishedalltrials_char_rounded = as.character(round(success_proportions_empirical_finishedalltrials, 2))

# determine the distribution in the theoretical sample
success_rates_theoretical = NA
for(i in 1:length(possible_success_rates)){
  success_rates_theoretical[i] = sum(success_proportions_theoretical_char_rounded == possible_success_rates_char[i])
}
success_rates_theoretical_prop = matrix(success_rates_theoretical/sum(success_rates_theoretical))

# determine the distribution in the empirical sample 
success_rates_empirical = NA
for(i in 1:length(possible_success_rates)){
  success_rates_empirical[i] = sum(success_proportions_empirical_finishedalltrials_char_rounded == possible_success_rates_char[i])
}
success_rates_empirical_prop = matrix(success_rates_empirical/sum(success_rates_empirical))



# Display the two distributions overlayed on Figure 3
fig_3_histogram_plot_data = as.data.frame(c(success_rates_theoretical_prop, success_rates_empirical_prop))
fig_3_histogram_plot_data = cbind(fig_3_histogram_plot_data, factor(c(rep("Expected if M0 is true", length(success_rates_theoretical_prop)), rep("Observed", length(success_rates_empirical_prop)))))
fig_3_histogram_plot_data = cbind(fig_3_histogram_plot_data, factor(rep(possible_success_rates_char, 2)))
names(fig_3_histogram_plot_data) = c("proportion", "group", "success")

figure_3 =  ggplot(fig_3_histogram_plot_data, aes(y = proportion, x = success, group = group))+
  geom_bar(aes(fill = group), alpha = 0.5, stat = "identity", position = "identity")+
  scale_fill_manual(values = c("darkgrey", "black")) +
  xlab("Successful guess rate") +
  ylab("Proportion") +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(angle = 90, color = "black", face = "bold", size = 12, margin = margin(1,0,0,0,"mm"), vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.title = element_text(size = 16))
figure_3

# earth mover's distance
EMD = emd2d(success_rates_theoretical_prop, success_rates_empirical_prop)


# Exploring the difference in success rate between those who did and those who did not finish all experimental trials.
data_BF_split_didnotfinishalltrials = data_BF_split[which(sapply(data_BF_split, nrow) != 18)]
success_proportions_empirical_didnotfinishalltrials = sapply(data_BF_split_didnotfinishalltrials, function(x) mean(as.logical(x[,"sides_match"])))

mean_success_rate_finishedalltrials = mean(success_proportions_empirical_finishedalltrials)
success_rate_finishedalltrials_SE = sd(success_proportions_empirical_finishedalltrials)/sqrt(length(success_proportions_empirical_finishedalltrials))
success_rate_finishedalltrials_CI_lb = mean_success_rate_finishedalltrials - 1.96*success_rate_finishedalltrials_SE
success_rate_finishedalltrials_CI_ub = mean_success_rate_finishedalltrials + 1.96*success_rate_finishedalltrials_SE

mean_success_rate_didnotfinishalltrials = mean(success_proportions_empirical_didnotfinishalltrials)
success_rate_didnotfinishalltrials_SE = sd(success_proportions_empirical_didnotfinishalltrials)/sqrt(length(success_proportions_empirical_didnotfinishalltrials))
success_rate_didnotfinishalltrials_CI_lb = mean_success_rate_didnotfinishalltrials - 1.96*success_rate_didnotfinishalltrials_SE
success_rate_didnotfinishalltrials_CI_ub = mean_success_rate_didnotfinishalltrials + 1.96*success_rate_didnotfinishalltrials_SE





