
# This script contains the code that was used to perform the
# analysis of operational characteristics in the Transparent Psi Project
# Please read the Statistical analysis section and the Operational characteristics
# section of the protocol to make sense of the parameters in the code


######################################################################
#                                                                    #
#                          Changeable parameters                     #
#                                                                    #
######################################################################

# Set parameters of the simulations below. Keeping the parameters as is will perform the 
# method 2 operational characteristics analysis mentioned in the protocol.
# The parameters are settable because this code was used to fine-tune our analysis and sampling plan.

# Note that you can specify several different parameters to try in one run
# For example if True_prob = c(0.5, 0.51) and BFs_to_try = c(25, 35) the code will perform simulations
# with four different parameters: 
# 1) True_prob = 0.5 and Inference_threshold_BF_low = 1/25, Inference_threshold_BF_high = 25; 
# 2) True_prob = 0.51 and Inference_threshold_BF_low = 1/25, Inference_threshold_BF_high = 25;
# 3) True_prob = 0.50 and Inference_threshold_BF_low = 1/35, Inference_threshold_BF_high = 35;
# 4) True_prob = 0.51 and Inference_threshold_BF_low = 1/35, Inference_threshold_BF_high = 35.
# The results are output in a table with the results of each of these simulations in separate rows 
# stored in the output_frame object

### Set the number of simulated samples to analyze
# we used 50000 simulations for each simulated true population effect size
# 50000 iterations with 10 stops per simulation runs for 100 minutes on i7 6600U 2.6GHz processor for each different true effec size simulated
iterations = 10000

### Set Model 0 probability of successful guesses. In our study this is 0.5
M0_prob = 0.50

### True population probability of success to try in the simulation
# this parameter can be set to a specific true population proportion such as 0.51, or as "sample"
# if set to "sample" it will sample true effect sizes from a prior distribution calculated from previous results (see below the y_and_N_prior_to_try parameter)
# setting this to c(0.5, "sample") will perform the method 2 analysis mentioned in ur protocol. With 50000 iterations this would take about 3 hours with an i 7 processor.
# setting this to c(0.45, 0.48, 0.49, 0.495, seq(0.498, 0.512, 0.001), 0.515, 0.520, 0.530, 0.560) will perform the methid 1 analysis mentioned in the protocol
# WARNING: performing the method 1 analysis with all effect sizes at once may run for several days on a regular machine
True_prob = c(0.5, 0.51)

### Set stopping parameters
# Minimal and maximal sample sizes. Sample size is given in total number of completed erotic trials.
# In the previous version of the proposal we proposed 30060 (minimum sample size), 37836, 45612, 53406, 61182, 68958, 76734, 84528, 92304, and a 100080 (maximum sample size)
min_trialnumber_to_try = c(37836)
max_trialnumber_to_try = 136080
# Number of interim analysis points. Analyses will be performed at equal sample size distances.
# in our study we use 10 interim analysis points.
times_to_stop_to_try = 5
# set BF threshold for inference decisions. Here the upper BF threshold needs to be specified as BF M0 vs M1.
# for example if this is set to 25, data collection will stop at eaching 25 BF or 1/25 BF.
BFs_to_try = c(25)

### Give information from previous study
# this information is used both for the sampling method (method 2 in the protocol) of analysis of 
# operational characteristics and for calculating replication Bayes factor, and the Bayesian parameter estimation robustness test. 
# Number of trials where prediction of target was successful and total number of trials needs to be set separated by a comma.
# Here we use data from Bem's experiment 1, "828, 1560" meaning 828 successes within 1560 trials
y_and_N_prior_to_try = "828, 1560" 

### Set parameters of robustness tests
# set the types of robustness test to perform after data collection stopped, 
# the type of test needs to be included in the "" separated by a comma, for example: "NHST_probtest, Bayes_Par_Est, Mixed_NHST"
# will perform an NHST proportion test, a Bayesian parameter estimation, and a mixed model analysis as robustness test
# the proportion test uses the SESIO and p threshold set below, the Bayesian parameter estimation uses the minimum_effect_threshold_Bayes_Par_Est_to_try as ROPE set below
# if set to NA no robustness test is conducted after data collection stopped
robustness_test_with = "NHST_probtest, Bayes_Par_Est"

# smallest effect size of interest in the NHST equivalence test 
# (this is used both in the NHST proportion test robustness test and in the NHST as the primary analysis, if any, see below)
SESOI = 0.01 
# p threshold for the NHST proportion test robustness test
Inference_threshold_robustness_NHST = 0.005 

# in the Bayesian parameter estimation robustness test this will determine the region of practical 
#equivalence (ROPE) interval. The ROPE is interpreted similarly to SESOI, but not entireli the same. 
# See Kruschke, J. K., & Liddell, T. M. (2017). The Bayesian New Statistics: Hypothesis testing, 
# estimation, meta-analysis, and power analysis from a Bayesian perspective. 
# Psychonomic Bulletin & Review, 1-29. 
minimum_effect_threshold_Bayes_Par_Est_to_try = 0.006
# this threshold is used to set the HDI width to check against the ROPE in the Bayesian parameter 
# estimation robustness test, if ths parameter is set to 0.05 for example, it means that we would 
# expect that 95% of the probability mass would be within the ROPE to accept a hypothesis
Inference_threshold_robustness_Bayes_Par_Est = 0.05 

### One- or two-sided test
# set whether tests are two or one sided and in which direction of one sided. 
# Can be set to "two.sided", "less", or "greater". In our study we use "greater".
alternative_hip_side = "greater"

### Determine the primary analysis
# This can be set to "prop_BF_replication", "prop_BF_uniform", "prop_BF_BUJ", "prop_NHST", "Mixed_NHST".
# if multiple analyses are specified, ALL analyses are performed, and support for M0 or M1 is only 
# declared in the primary analyses, if ALL analyses yield the same conclusion, otherwise. In this case,
# data collection is stopped at the interim stopping point. If there is no consensus in supporting either
# M0 or M1, data collection continuous and all analyses are repeated at the next stopping point. In case of
# NHST analyses, critical p-value threshold is corrected for multiple comparisons at every interim stopping
# point to compensate for inflated inferential error rates due to "data peaking".
Primary_analysis = c("prop_BF_replication", "prop_BF_uniform", "prop_BF_BUJ", "Mixed_NHST")

### Set parameters of NHST test if Primary_analysis is set to "both" or "NHST"
# if doing NHST as a main analysis what is the p-threshold of accepting the hypotheses
# the code authomatically adjust for multiple comparisons due to sequential analysis,
# uses Bonferroni correction, so this does not have to be factored in.  
p_threshold_to_try = 0.005

### determine the extent of personal differences
# there may be peronal differences in the probability of successful guesses
# by setting this parameter to 0 (zero), we simulate no personal differences
# by setting this parameter to higher than zero, the code will determin a "talen" value for each person,
# which is computed by adding a random number with zero mean and this parameter as the SD
# to True_prob, which is the true mean "talent" (probability of a successful guesses) in the population
# this is done by the talent_randomizer custom function.
# setting this to 0.15 will result in 5% of the simulated sample being extremely lucky (75%+ chance of success)
# and 5% of the population being extremely unlucky (25%- chance of success), and 10% of the population having
# 70%+ success rate, while 10% of the population having 30%- success rate.
SD_personal_differences_to_try = c(0) #0 to simulate no intercorrelation, or 0.15 to simulate reasonably high intercorrelation








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
# we thank Richard Morey for his help in developing these functions


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
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, rel.tol = 1e-15)[[1]]
}

BF01_beta = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y, N, null_prob) / margLike_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)
},"y")



######################################################################
#                       Other supporting functions                   #
######################################################################

### Function to find mode

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


### Function to round to nearest x function
# this is used to make sure that interim analysis points are at the end of full study sessions
# this is just to make interim analysis points in the simulation more logical, 
# in the real study interim analysis points won't be fixed to completed study sessions

mround <- function(x,base){ 
  base*round(x/base) 
} 


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


# randomize psi talent if personal differences are simulated
talent_randomizer <- function(True_prob, SD_personal_differences){
  talent = True_prob + rnorm(mean = 0, sd = SD_personal_differences, n = 1)
  if(talent<0.001){talent = 0.001}else if(talent>0.999){talent = 0.999}
  return(talent)
}


# convert logit to probability
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}



######################################################################
#                         Simulation function                        #
######################################################################

# function performing sampling, data analysis, determining inference, and robustness

simul <- function(M0_prob, True_prob, True_prob_samp, max_num_trials,
                  y_prior, N_prior, alternative_hip_side, 
                  Inference_threshold_BF_low, Inference_threshold_BF_high, when_to_check, 
                  Primary_analysis, minimum_effect_threshold_NHST, minimum_effect_threshold_Bayes_Par_Est,
                  Inference_threshold_NHST,robustness_test_with,
                  Inference_threshold_robustness_NHST, Inference_threshold_robustness_Bayes_Par_Est,
                  SD_personal_differences){
  
  # Decide whether we use a fixed effect size or effect size sampled from the prior_density_distribution
  if(True_prob == "sample"){
    True_prob = sample(True_prob_samp, 1)
  }
    
  
  # one or two sided tests
  # if alternative_hip_side == "greater" or "less", a one-sided test will be performed 
  if(alternative_hip_side == "greater"){
    alter_Interval_proptest = c(0.5,1)
    equivalence_test_side = "less"
  } else if(alternative_hip_side == "less"){
    alter_Interval_proptest = c(0,0.5)
    equivalence_test_side = "greater"
  } else {
    alter_Interval_proptest = c(0,1)
    equivalence_test_side = "two.sided"
  }

  # vector containing the names of the robustness tests used in this simulation
  if(!is.na(robustness_test_with)){
    robustness_test_types = unlist(strsplit(as.character(robustness_test_with), split=", "))
  }
  
  
  ##########################################################
  #         Outputs of the simulation functions            #
  ##########################################################

  # BF test outputs
  BF_replication = NA # replication Bayes Factor 
  BF_uniform  = NA # Bayes Factor with uniform prior
  BF_BUJ = NA # Bayes Factor with BUJ prior (see below)
  # True_prob is also part of the output
  data_length_at_stop = NA # the number of trials where the simulated study ended in the Bayes Factor analysis as a main analysis 
  Final_primary_inference = NA # inference based on results in the Bayes Factor analysis as a main analysis 
  
  #NHST test outputs
  p_greater_then_nullprob = NA # p-value from the NHST equality test as a main analysis
  p_less_than_minimum_effect_threshold = NA # p-value from the NHST equivalence test as a main analysis
  inference_NHST = NA # inference based on results of the NHST test as a main analysis
  
  # robustness test outputs
  robustness_NHST_probtest_p_greater_then_nullprob = NA # p-value from the NHST equality test as a robustness analysis
  robustness_NHST_probtest_p_less_than_minimum_effect_threshold = NA # p-value from the NHST equivalence test as a robustness analysis
  inference_robustness = NA # inference based on results of all robustness analyses combined
  inference_robustness_NHST = NA # inference based on results of the NHST test as a robustness analysis
  inference_robustness_Bayes_Par_Est = NA # inference based on results of the Bayesian Parameter Estimation robustness test
  inference_robustness_mixed = NA # inference based on results of the NHST test as a robustness analysis
  Robust = NA # whether the test inferences were the same in the main and robustness analyses
  

  
  ##########################################################
  #                 Draw simulated sample                  #
  ##########################################################
  # generate data with True_prob, the simulated data asumes each trial to be independent
  # see reasoning on this decision in the Treating each trial as independent in data analysis 
  # subsection in the Additional Methodological Considerations section of the protocol
  # also, see demonstration that participant effects (clustering) does not affect the results of this analysis
  # via this R demo: https://github.com/kekecsz/Transparent_Psi_Project_scripts/blob/master/Demonstrate%20insensitivity%20to%20clustering%20of%20the%20effect.R
  
  # only generate personal differences if True_prob != 0.5
  
  if(any(robustness_test_types == "Mixed_NHST")|any(Primary_analysis == "Mixed_NHST")){
    if((SD_personal_differences != 0) & (True_prob != 0.5)){
      data_all_M1 <- as.vector(replicate(mround(max_num_trials, base = 18)/18, rbinom(18, size = 1, prob=talent_randomizer(True_prob = True_prob, SD_personal_differences = SD_personal_differences))))
    } else {
      data_all_M1 <- rbinom(mround(max_num_trials, base = 18), size = 1, prob=True_prob)
    }
    data_all_M1_dataframe = as.data.frame(cbind(data_all_M1, rep(1:(mround(max_num_trials, base = 18)/18), each = 18)))
    data_all_M1_dataframe[,2] = as.factor(data_all_M1_dataframe[,2])
    names(data_all_M1_dataframe) = c("outcome","ID")
  } else {
    if((SD_personal_differences != 0) & (True_prob != 0.5)){
      data_all_M1 <- as.vector(replicate(max_num_trials/18, rbinom(18, size = 1, prob=talent_randomizer(True_prob = True_prob, SD_personal_differences = SD_personal_differences))))
    } else {
      data_all_M1 <- rbinom(max_num_trials, size = 1, prob=True_prob)
    }
  }
  

    
    ##########################################################
    # Bayes Factor test as a main/primary hypothesis testing #
    ##########################################################
    
    ### Analysis on the dataset in which M1 is true
    # compute Bayes factor on the simulated data where M1 is true
    # Bayes factor is computed at the segments of the dataset specified by when_to_check. 
    # if the BF thresholds for stopping are reached (set in the parameters Inference_threshold_BF_low and Inference_threshold_BF_high), the loop stops

  comparisons_prop_NHST = 0
  comparisons_Mixed_NHST = 0
  
    for(i in when_to_check){
      
      # dataset used in the Bayes factor calculation, it includes data from the beggining the
      # data_all_M1 full simulated dataset, the size of the sample it contains increases in the
      # loop until the BF threshold or the maxumum sample size is reached 
      data_BF = data_all_M1[1:i]
      data_dataframe = data_all_M1_dataframe[1:i,]
      
      # number of successes and total N of trials
      successes = sum(data_BF)
      total_N = length(data_BF)
      
      Primary_inferences = NA
      Primary_inferences_num = 0
      
      if(any(Primary_analysis == "prop_BF_replication")){
        Primary_inferences_num = Primary_inferences_num + 1
      ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
      
      BF <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = alter_Interval_proptest, null_prob = M0_prob) #numbers higher than 1 support the null
      
      if(Inference_threshold_BF_low >= BF) {Primary_inferences[Primary_inferences_num] = "M1"
      } else if(Inference_threshold_BF_high <= BF) {Primary_inferences[Primary_inferences_num] = "M0"
      } else {Primary_inferences[Primary_inferences_num] = "Inconclusive"}
      }
      
      if(any(Primary_analysis == "prop_BF_uniform")){
        Primary_inferences_num = Primary_inferences_num + 1
        ### Bayes factor with uniform prior
        # using a non-informative flat prior distribution with alpha = 1 and beta = 1
        
        BF <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = alter_Interval_proptest, null_prob = M0_prob) #numbers higher than 1 support the null

        if(Inference_threshold_BF_low >= BF) {Primary_inferences[Primary_inferences_num] = "M1"
        } else if(Inference_threshold_BF_high <= BF) {Primary_inferences[Primary_inferences_num] = "M0"
        } else {Primary_inferences[Primary_inferences_num] = "Inconclusive"}
      }
      
      if(any(Primary_analysis == "prop_BF_BUJ")){
        Primary_inferences_num = Primary_inferences_num + 1
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
        
        BF <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = alter_Interval_proptest, null_prob = M0_prob) #numbers higher than 1 support the null
        
        if(Inference_threshold_BF_low >= BF) {Primary_inferences[Primary_inferences_num] = "M1"
        } else if(Inference_threshold_BF_high <= BF) {Primary_inferences[Primary_inferences_num] = "M0"
        } else {Primary_inferences[Primary_inferences_num] = "Inconclusive"}
      }
            
      

      if(any(Primary_analysis == "prop_NHST")){
        Primary_inferences_num = Primary_inferences_num + 1      
        comparisons_prop_NHST = comparisons_prop_NHST + 2 #because we do two tests at each stop point
        
        # equality test
        p_greater_then_nullprob = prop.test(x = successes, n = total_N, p = M0_prob, alternative = alternative_hip_side)$p.value
        p_greater_then_nullprob = p.adjust(p = p_greater_then_nullprob, method = "bonferroni", n = comparisons_prop_NHST)
        
        # equivalence test
        if(equivalence_test_side == "less"){
          p_less_than_minimum_effect_threshold = prop.test(x = successes, n = total_N, p = M0_prob+minimum_effect_threshold_NHST, alternative = equivalence_test_side)$p.value
          p_less_than_minimum_effect_threshold = p.adjust(p = p_less_than_minimum_effect_threshold, method = "bonferroni", n = comparisons_prop_NHST)
        } else if(equivalence_test_side == "greater"){
          p_less_than_minimum_effect_threshold = prop.test(x = successes, n = total_N, p = M0_prob-minimum_effect_threshold_NHST, alternative = equivalence_test_side)$p.value
          p_less_than_minimum_effect_threshold = p.adjust(p = p_less_than_minimum_effect_threshold, method = "bonferroni", n = comparisons_prop_NHST)
        } else if(equivalence_test_side == "two.sided"){
          p_equivalence_upper_bound_M1_true = prop.test(x = successes, n = total_N, p = M0_prob+minimum_effect_threshold_NHST, alternative = "less")$p.value
          p_equivalence_upper_bound_M1_true = p.adjust(p = p_equivalence_upper_bound_M1_true, method = "bonferroni", n = comparisons_prop_NHST)
          p_equivalence_lower_bound_M1_true = prop.test(x = successes, n = total_N, p = M0_prob-minimum_effect_threshold_NHST, alternative = "greater")$p.value
          p_equivalence_lower_bound_M1_true = p.adjust(p = p_equivalence_lower_bound_M1_true, method = "bonferroni", n = comparisons_prop_NHST)
          p_less_than_minimum_effect_threshold = max(c(p_equivalence_upper_bound_M1_true, p_equivalence_lower_bound_M1_true))
        } 
        
        
        if(Inference_threshold_NHST > p_less_than_minimum_effect_threshold){Primary_inferences[Primary_inferences_num] = "M0"
        } else if(Inference_threshold_NHST > p_greater_then_nullprob){Primary_inferences[Primary_inferences_num] = "M1"
        } else {Primary_inferences[Primary_inferences_num] = "Inconclusive"}
      }
        
        if(any(Primary_analysis == "Mixed_NHST")){
          Primary_inferences_num = Primary_inferences_num + 1      
          comparisons_Mixed_NHST = comparisons_Mixed_NHST + 2 #because we do two tests at each stop point
          
          
          mod_mixed = glmer(outcome ~ 1 + (1|ID), data = data_dataframe, family = "binomial")
          estimate_mixed = summary(mod_mixed)$coefficients[1,1]
          se_mixed = summary(mod_mixed)$coefficients[1,2]
          
          wald_ci_mixed_logit <- c(estimate_mixed - se_mixed* qnorm(1-((Inference_threshold_robustness_NHST/comparisons_Mixed_NHST)/2)),
                                   estimate_mixed + se_mixed* qnorm(1-((Inference_threshold_robustness_NHST/comparisons_Mixed_NHST)/2)))
          wald_ci_mixed = logit2prob(wald_ci_mixed_logit)
          
          wald_ci_mixed_lb = wald_ci_mixed[1]
          wald_ci_mixed_ub = wald_ci_mixed[2]
          
          if(alternative_hip_side == "less"){
            minimum_effect = M0_prob-minimum_effect_threshold_NHST
            if(wald_ci_mixed_lb > minimum_effect){Primary_inferences[Primary_inferences_num] = "M0"
            } else if(wald_ci_mixed_ub < M0_prob){Primary_inferences[Primary_inferences_num] = "M1"
            } else {Primary_inferences[Primary_inferences_num] = "Inconclusive"}
          } else if (alternative_hip_side == "greater"){
            minimum_effect = M0_prob+minimum_effect_threshold_NHST
            if(wald_ci_mixed_ub < minimum_effect){Primary_inferences[Primary_inferences_num] = "M0"
            } else if(wald_ci_mixed_lb > M0_prob){Primary_inferences[Primary_inferences_num] = "M1"
            } else {Primary_inferences[Primary_inferences_num] = "Inconclusive"}
          } else if (alternative_hip_side == "two.sided"){
            minimum_effect_ub = M0_prob+minimum_effect_threshold_NHST
            minimum_effect_lb = M0_prob-minimum_effect_threshold_NHST
            if(wald_ci_mixed_ub < minimum_effect_ub & wald_ci_mixed_lb > minimum_effect_lb){Primary_inferences[Primary_inferences_num] = "M0"
            } else if(wald_ci_mixed_lb > M0_prob | wald_ci_mixed_ub < M0_prob){Primary_inferences[Primary_inferences_num] = "M1"
            } else {Primary_inferences[Primary_inferences_num] = "Inconclusive"}
          }
          
        }
        
        
        #================================================================#
        #                    Main analysis inference                     #
        #================================================================#
        
        if(all(Primary_inferences == "M1")){Final_primary_inference = "M1"
        break
        } else if(all(Primary_inferences == "M0")){Final_primary_inference = "M0"
        break
        } else {Final_primary_inference = "Inconclusive"}
        
      }
      
      # store sample size at last stopping point
      data_length_at_stop = total_N
      
      
      ##################################################
      #                 Robustness tests               #
      ##################################################
      
      #================================================================#
      #               Robustness test of BF results with NHST          #
      #================================================================#
      
      # robustness of BF results is tested with NHST proportion tests
      # here we perform both an equality test and an equivalence test to draw statistical inference
      
      if(any(robustness_test_types == "NHST_probtest")){
        # equality test
        robustness_NHST_probtest_p_greater_then_nullprob = prop.test(x = successes, n = total_N, p = M0_prob, alternative = alternative_hip_side)$p.value
        
        # equivalence test
        if(equivalence_test_side == "less"){
          robustness_NHST_probtest_p_less_than_minimum_effect_threshold = prop.test(x = successes, n = total_N,
                                                                                    p = M0_prob+minimum_effect_threshold_NHST, 
                                                                                    alternative = equivalence_test_side)$p.value
        } else if(equivalence_test_side == "greater"){
          robustness_NHST_probtest_p_less_than_minimum_effect_threshold = prop.test(x = successes, 
                                                                                    n = total_N, p = M0_prob-minimum_effect_threshold_NHST, 
                                                                                    alternative = equivalence_test_side)$p.value
        } else if(equivalence_test_side == "two.sided"){
          robustness_NHST_probtest_p_less_than_minimum_effect_threshold_upper = prop.test(x = successes, n = total_N, 
                                                                                          p = M0_prob+minimum_effect_threshold_NHST, alternative = "less")$p.value
          robustness_NHST_probtest_p_less_than_minimum_effect_threshold_lower = prop.test(x = successes, n = total_N, 
                                                                                          p = M0_prob-minimum_effect_threshold_NHST, 
                                                                                          alternative = "greater")$p.value
          robustness_NHST_probtest_p_less_than_minimum_effect_threshold = max(c(robustness_NHST_probtest_p_less_than_minimum_effect_threshold_upper,
                                                                                robustness_NHST_probtest_p_less_than_minimum_effect_threshold_lower))
        }
        
        # making inference decision
        if(!is.na(Inference_threshold_robustness_NHST)){
          if(Inference_threshold_robustness_NHST > robustness_NHST_probtest_p_less_than_minimum_effect_threshold){inference_robustness_NHST = "M0" 
          } else if(Inference_threshold_robustness_NHST > robustness_NHST_probtest_p_greater_then_nullprob){inference_robustness_NHST = "M1"
          } else {inference_robustness_NHST = "Inconclusive"}
        }
      }
      
      
      
      #=======================================================================#
      #   Robustness test of BF results with Bayesian parameter estimation    #
      #=======================================================================#
      
      # robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
      # the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
      # (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
      # analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 
      
      if(any(robustness_test_types == "Bayes_Par_Est")){
        
        # calculate posterior distribution using beta distribution
        prior_alpha = y_prior + 1
        prior_beta = N_prior-y_prior+1
        
        posterior_alpha = prior_alpha + successes
        posterior_beta = prior_beta + total_N - successes
        
        scale = seq(0, 1, length = 1001)
        posterior_density = dbeta(scale, posterior_alpha, posterior_beta)
        
        # calculate HDI for the posterior distribution
        if(alternative_hip_side == "less" | alternative_hip_side == "greater"){
          hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)
        } else {hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est, n_samples = 1e6)}
        
        # parameters for decision making
        HDI_lb = hdi_result[2]
        HDI_ub = hdi_result[3]
        
        # making inference decision
        if(alternative_hip_side == "less"){
          ROPE = M0_prob-minimum_effect_threshold_Bayes_Par_Est
          if(HDI_lb > ROPE){inference_robustness_Bayes_Par_Est = "M0"
          } else if(HDI_ub < ROPE){inference_robustness_Bayes_Par_Est = "M1"
          } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
        } else if (alternative_hip_side == "greater"){
          ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
          if(HDI_ub < ROPE){inference_robustness_Bayes_Par_Est = "M0"
          } else if(HDI_lb > ROPE){inference_robustness_Bayes_Par_Est = "M1"
          } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
        } else if (alternative_hip_side == "two.sided"){
          ROPE_ub = M0_prob+minimum_effect_threshold_Bayes_Par_Est
          ROPE_lb = M0_prob-minimum_effect_threshold_Bayes_Par_Est
          if(HDI_ub < ROPE_ub & HDI_lb > ROPE_lb){inference_robustness_Bayes_Par_Est = "M0"
          } else if(HDI_lb > ROPE_ub | HDI_ub < ROPE_lb){inference_robustness_Bayes_Par_Est = "M1"
          } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
        }
        
      }
      
      
      #====================================================================================#
      #   Robustness test of using mixed models to account for personal differences        #
      #====================================================================================#
      
      if(any(robustness_test_types == "Mixed_NHST")){
        
        
        mod_mixed = glmer(outcome ~ 1 + (1|ID), data = data_dataframe, family = "binomial")
        estimate_mixed = summary(mod_mixed)$coefficients[1,1]
        se_mixed = summary(mod_mixed)$coefficients[1,2]
        
        wald_ci_mixed_logit <- c(estimate_mixed - se_mixed* qnorm(1-(Inference_threshold_robustness_NHST/2)),
                                 estimate_mixed + se_mixed* qnorm(1-(Inference_threshold_robustness_NHST/2)))
        wald_ci_mixed = logit2prob(wald_ci_mixed_logit)
        
        wald_ci_mixed_lb = wald_ci_mixed[1]
        wald_ci_mixed_ub = wald_ci_mixed[2]
        
        if(alternative_hip_side == "less"){
          minimum_effect = M0_prob-minimum_effect_threshold_NHST
          if(wald_ci_mixed_lb > minimum_effect){inference_robustness_mixed = "M0"
          } else if(wald_ci_mixed_ub < M0_prob){inference_robustness_mixed = "M1"
          } else {inference_robustness_mixed = "Inconclusive"}
        } else if (alternative_hip_side == "greater"){
          minimum_effect = M0_prob+minimum_effect_threshold_NHST
          if(wald_ci_mixed_ub < minimum_effect){inference_robustness_mixed = "M0"
          } else if(wald_ci_mixed_lb > M0_prob){inference_robustness_mixed = "M1"
          } else {inference_robustness_mixed = "Inconclusive"}
        } else if (alternative_hip_side == "two.sided"){
          minimum_effect_ub = M0_prob+minimum_effect_threshold_NHST
          minimum_effect_lb = M0_prob-minimum_effect_threshold_NHST
          if(wald_ci_mixed_ub < minimum_effect_ub & wald_ci_mixed_lb > minimum_effect_lb){inference_robustness_mixed = "M0"
          } else if (wald_ci_mixed_lb > M0_prob | wald_ci_mixed_ub < M0_prob){inference_robustness_mixed = "M1"
          } else {inference_robustness_mixed = "Inconclusive"}
        }
        
      }
 

  
  
  
  #=======================================================================#
  #      Determine final inference of all robustness tests combined       #
  #=======================================================================#
  
  # the main analysis inference is only robust if all robustness tests came to the same inference as the main test
  if(!is.na(robustness_test_with)){
    inferences = na.omit(c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est, inference_robustness_mixed))
    inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}
  
    Robust = if(Final_primary_inference == inference_robustness){"robust"} else {"not robust"}
  }

  
  
    
  ##################################################
  #           Output simulation results            #
  ##################################################
  
    output <- c(BF_replication,
                BF_uniform,
                BF_BUJ,
                True_prob,
                data_length_at_stop,
                p_greater_then_nullprob,
                p_less_than_minimum_effect_threshold,
                robustness_NHST_probtest_p_greater_then_nullprob,
                robustness_NHST_probtest_p_less_than_minimum_effect_threshold,
                Final_primary_inference,
                inference_robustness,
                inference_NHST,
                inference_robustness_NHST,
                inference_robustness_Bayes_Par_Est,
                inference_robustness_mixed,
                Robust)
  
  return(output)
}


######################################################################
#                                                                    #
#                           Load packages                            #
#                                                                    #
######################################################################

library(parallel) # for parallel processing on multiple cores
cl <- makeCluster(detectCores()) # set up working clusters (run task on all but one processor cores)
library(pbapply) # to progess bar
## load libraries on workers
clusterEvalQ(cl, library(HDInterval)) # needed to calcluate HDI credible intervals in the Bayesian parameter estimation robustness test
clusterEvalQ(cl, library(lme4)) #for mixed model robustness test




start_time <- Sys.time() # monitoring running time

######################################################################
#                                                                    #
#                          Running the simulaton                     #
#                                                                    #
######################################################################

######################################################################
#                   Simulation and output                            #
######################################################################
# this part of the code is designated to set up the simulation while allowing for
# multiple parameters to be tested in the same simulation set in the changable parameters above.
# for example, this code allows to test operational charateristics with both 0.5, 0.505 and 0.51 true 
# population effect sizes in one run, and output results in an output table



# create a dataframe which will house the final results of the simulation series
output_frame <- expand.grid(iterations, M0_prob, True_prob, alternative_hip_side, p_threshold_to_try, SESOI, minimum_effect_threshold_Bayes_Par_Est_to_try, robustness_test_with, BFs_to_try, y_and_N_prior_to_try, min_trialnumber_to_try, max_trialnumber_to_try, times_to_stop_to_try, Inference_threshold_robustness_NHST, Inference_threshold_robustness_Bayes_Par_Est, SD_personal_differences_to_try)
names(output_frame) = c("iterations", "M0_prob", "True_prob", "alternative_hip_side", "p_threshold_to_try", "SESOI", "minimum_effect_threshold_Bayes_Par_Est_to_try", "robustness_test_with", "Inference_threshold_BF_high", "y_and_N_prior_to_try", "min_trial_n", "max_trial_n", "times_to_stop", "Inference_threshold_robustness_NHST", "Inference_threshold_robustness_Bayes_Par_Est", "SD_personal_differences")
output_frame[,c("Supported_M1", "Supported_M0", "Inconclusive", 
                "Mode_num_trials_true", "Stopped_at_min_trial_n", "Stopped_at_max_trial_n", "Stopped_75percent_of_studies_at_or_before_this_n",
                "Robust", "Primary_analysis")] <- NA

BFs_output = list(NA)

output_frame[, "Primary_analysis"] = paste(Primary_analysis, collapse = " ")

# This for loop runs the simulation series with all the possible combinations of the parameters to try given above
for(i in 1:nrow(output_frame)){
  
  M0_prob = output_frame[i,"M0_prob"]

  # set BF optional stopping thresholds which will be also used in the power analysis
  Inference_threshold_BF_high = output_frame[i,"Inference_threshold_BF_high"]
  Inference_threshold_BF_low = 1/Inference_threshold_BF_high

  # set prior distribution and scale if applicable
  y_prior = as.numeric(unlist(strsplit(as.character(output_frame[i,"y_and_N_prior_to_try"]), split=",")))[1]
  N_prior = as.numeric(unlist(strsplit(as.character(output_frame[i,"y_and_N_prior_to_try"]), split=",")))[2]
  
  # set whether tests are two or one sides and in which direction of ine sided
  alternative_hip_side = as.character(output_frame[i,"alternative_hip_side"])
  
  # set true probability of success
  if(output_frame[i,"True_prob"] == "sample"){
    
    True_prob = output_frame[i,"True_prob"]
    
    scale_samp = runif(100000,0,1)
    
    # one or two sided tests
    # if alternative_hip_side == "greater" or "less", a one-sided test will be performed 
    if(alternative_hip_side == "greater"){
      interval_samp = c(0.5,1)
    } else if(alternative_hip_side == "less"){
      interval_samp = c(0,0.5)
    } else {
      interval_samp = c(0,1)
    }
    
    alpha_samp = y_prior+1
    beta_samp = N_prior-y_prior+1
    
    normalizing_constant = diff(pbeta(interval_samp, alpha_samp, beta_samp))
    range_samp = scale_samp>interval_samp[1] & scale_samp<interval_samp[2]
    
    prior_density_samp = range_samp*dbeta(scale_samp, alpha_samp, beta_samp)/normalizing_constant
    maxDens = max(prior_density_samp)
    
    accepted = ifelse(runif(100000,0,1) < prior_density_samp / maxDens, TRUE, FALSE)
    # Plot the result
    
    True_prob_samp = scale_samp[accepted]
    
  } else {
    True_prob = as.numeric(as.character(output_frame[i,"True_prob"]))
    True_prob_samp = NA
  }
  
    
  #set p_threhold for the main NHST analyses if any
  Inference_threshold_NHST = output_frame[i,"p_threshold_to_try"]
  SESOI = output_frame[i,"SESOI"]
  
  # set minimum and maximum sample size
  min_trialnumber = output_frame[i,"min_trial_n"]
  max_trialnumber = output_frame[i,"max_trial_n"]
  
  # set the final stopping points
  when_to_check = mround(x = seq(min_trialnumber, max_trialnumber, (max_trialnumber - min_trialnumber)/(output_frame[i,"times_to_stop"]-1)), base = 18)
  
  # set inference criteria for robustness test
  Inference_threshold_robustness_NHST =  output_frame[i,"Inference_threshold_robustness_NHST"]
  Inference_threshold_robustness_Bayes_Par_Est =  output_frame[i,"Inference_threshold_robustness_Bayes_Par_Est"]
  
  minimum_effect_threshold_Bayes_Par_Est = output_frame[i,"minimum_effect_threshold_Bayes_Par_Est_to_try"]
  
  SD_personal_differences =  output_frame[i,"SD_personal_differences"]
     
  ######################################################################
  #                      Run the simulation                            #
  ######################################################################
  
  # make objects available for working clusters (all cores)
  clusterExport(cl,c(
    "fullAlt_beta",
    "normalize_beta",
    "restrictedAlt_beta",
    "margLike_beta",
    "BF01_beta",
    "Mode",
    "mround",
    "mode_HDI",
    "talent_randomizer",
    "logit2prob",
    "simul",
    
    "iterations",
    "M0_prob",
    "True_prob",
    "True_prob_samp",
    "max_trialnumber",
    "y_prior",
    "N_prior",
    "alternative_hip_side",
    "Inference_threshold_BF_low",
    "Inference_threshold_BF_high",
    "when_to_check",
    "Primary_analysis",
    "SESOI",
    "minimum_effect_threshold_Bayes_Par_Est",
    "Inference_threshold_NHST",
    "robustness_test_with",
    "Inference_threshold_robustness_NHST",
    "Inference_threshold_robustness_Bayes_Par_Est",
    "SD_personal_differences"
  ))
  
  
  ######################################################################
  #                         Simulation function                        #
  ######################################################################
  
  # function performing sampling, data analysis, determining inference, and robustness
  
  print(paste("Task ", i, " of ", nrow(output_frame), sep = ""))
  out <- pbsapply(cl = cl, X = 1:iterations, FUN = function(i,...) { simul(M0_prob = M0_prob,
                                                                           True_prob = True_prob,
                                                                           True_prob_samp = True_prob_samp,
                                                                           max_num_trials = max_trialnumber,
                                                                           y_prior = y_prior,
                                                                           N_prior = N_prior,
                                                                           alternative_hip_side = alternative_hip_side,
                                                                           Inference_threshold_BF_low = Inference_threshold_BF_low,
                                                                           Inference_threshold_BF_high = Inference_threshold_BF_high,
                                                                           when_to_check = when_to_check,
                                                                           Primary_analysis = Primary_analysis,
                                                                           minimum_effect_threshold_NHST = SESOI,
                                                                           minimum_effect_threshold_Bayes_Par_Est = minimum_effect_threshold_Bayes_Par_Est,
                                                                           Inference_threshold_NHST = Inference_threshold_NHST,
                                                                           robustness_test_with = robustness_test_with,
                                                                           Inference_threshold_robustness_NHST = Inference_threshold_robustness_NHST,
                                                                           Inference_threshold_robustness_Bayes_Par_Est = Inference_threshold_robustness_Bayes_Par_Est,
                                                                           SD_personal_differences = SD_personal_differences) } )
  
  
  

  ######################################################################
  #          Prepare output table with results of simulation           #
  ######################################################################
  
  # output of the simulations are stored
  sim_output <- as.data.frame(t(out))
  names(sim_output) = c("BF_replication",
                        "BF_uniform",
                        "BF_BUJ",
                        "True_prob",
                        "data_length_at_stop",
                        "p_greater_then_nullprob",
                        "p_less_than_minimum_effect_threshold",
                        "robustness_NHST_probtest_p_greater_then_nullprob",
                        "robustness_NHST_probtest_p_less_than_minimum_effect_threshold",
                        "Final_primary_inference",
                        "inference_robustness",
                        "inference_NHST",
                        "inference_robustness_NHST",
                        "inference_robustness_Bayes_Par_Est",
                        "inference_robustness_mixed",
                        "Robust")
  

  
  BFs_output[[i]] = as.data.frame(cbind(as.numeric(as.character(sim_output$BF_replication)), as.numeric(as.character(sim_output$BF_uniform)), as.numeric(as.character(sim_output$BF_BUJ)), as.numeric(as.character(sim_output$True_prob)), as.numeric(as.character(sim_output$data_length_at_stop))))
  names(BFs_output[[i]]) = c("BF_replication", "BF_uniform", "BF_BUJ", "True_prob", "data_length_at_stop")
  
  sim_output[,"data_length_at_stop"] = as.numeric(as.character(sim_output[,"data_length_at_stop"]))
  output_frame[i,"Robust"] = sum(sim_output[,"Robust"]== "robust")/nrow(sim_output)
  

    output_frame[i,"Supported_M1"] = sum(sim_output[,"Final_primary_inference"] == "M1")/nrow(sim_output)
    output_frame[i,"Supported_M0"] = sum(sim_output[,"Final_primary_inference"] == "M0")/nrow(sim_output)
    output_frame[i,"Inconclusive"] = sum(sim_output[,"Final_primary_inference"] == "Inconclusive")/nrow(sim_output)
    output_frame[i,"Mode_num_trials_true"] = Mode(sim_output[,"data_length_at_stop"])
    output_frame[i,"Stopped_at_min_trial_n"] = length(which(sim_output[,"data_length_at_stop"] == output_frame[i,"min_trial_n"]))/nrow(sim_output)
    output_frame[i,"Stopped_at_max_trial_n"] = length(which(sim_output[,"data_length_at_stop"] == output_frame[i,"max_trial_n"]))/nrow(sim_output)
    output_frame[i,"Stopped_75percent_of_studies_at_or_before_this_n"] = quantile(sim_output[,"data_length_at_stop"])["75%"]

}

#stop the cluster
stopCluster(cl)




### inspect results in a table
output_frame
