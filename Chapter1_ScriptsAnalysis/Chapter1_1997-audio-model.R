#! usr/bin/R
# Chapter1_1997-audio-model.R
# H. Muller
# 2019-01-18

set.seed(123) # set this seed to reproduce the reported results


#############################
###     Load packages     ###
#############################

wants <- c('lme4', 'dplyr', 'ggplot2', 'tidyverse', 'tidyr')
has <- wants %in% rownames(installed.packages())
if(any(!has))install.packages(wants[!has])
lapply(wants, library, character.only = TRUE)
beep()


############################
###   Define functions   ###
############################

# Compute the time needed for a word to be recognized vie the look-up route
RT_lookup <- Vectorize(function(sg_t_act, epsilon) {
  act_vec <- rnorm(500, sg_t_act, sg_t_act/4)
  eps_vec <- rnorm(500, epsilon, epsilon/4)
  RT = sum(act_vec + eps_vec) / 500
  # RT = sum(act_vec) / 500
  my_out <- paste0(as.character(RT), '-lkp')
  return(my_out)
}) 

# Compute the time needed for a word to be recognized vie the decomposition route
RT_decompose <- Vectorize(function(sg_t_act, pl_t_act, epsilon, t_parsing) {
  act1_vec <- rnorm(500, pl_t_act, pl_t_act/4)
  act2_vec <- rnorm(500, sg_t_act, sg_t_act/4)
  eps_vec <- rnorm(500, epsilon, epsilon/4)
  par_vec <- rnorm(500, t_parsing, t_parsing/4)
  my_rt1 = sum(eps_vec + act1_vec) / 500
  # my_rt1 = act1_vec
  my_rt2 = sum(eps_vec + act2_vec + par_vec) / 500
  # my_rt2 = act2_vec + par_vec
  RT = pmin(my_rt1, my_rt2)
  if(RT==my_rt1){
    route <- 'lkp'
  } else {
    route <- 'dcmp'
  }
  my_out <- paste0(as.character(RT), '-', route)
  return(my_out)
}) 

# Scale the predicted reaction times in model time to the observed times
scale_RTmodel <- function(t_real, t_model) {
  rt_min = min(t_real)
  rt_max = max(t_real)
  division = (t_model - min(t_model)) / (max(t_model) - min(t_model))
  rt_scale = rt_min + division * (rt_max - rt_min)
}

# Scale the parsing penalty in model time to the observed timed
scale_RTmodel_p <- function(t_real, t_model, p) {
  rt_min = min(t_real)
  rt_max = max(t_real)
  division = (p - min(t_model)) / (max(t_model) - min(t_model))
  rt_scale = rt_min + division * (rt_max - rt_min)
}

# Compute the model's log-likelihood
my_loglik <- function(y, yhat, number, pp) {
  
  a <- ifelse(number=='sing',1,0)
  
  phi1 <- dnorm(y, yhat, sqrt(0.25))
  phi2 <- dnorm(y, yhat, sqrt(0.25))
  phi3 <- dnorm(y, yhat, sqrt(0.5))
  
  Phi2 <- pnorm(y, yhat, sqrt(0.25))
  Phi3 <- pnorm(y, yhat, sqrt(0.5))
  
  f1 <- phi1
  f2 <- phi2*Phi3 + phi3*Phi2
  
  the_loglik <- sum(a*log(f1)) + sum((1-a)*log(f2))
  return(the_loglik)
}

# Compute the model's AIC
get_aic <- function(dual_loglik) {
  my_aic <- 2*3 - 2*dual_loglik
}

# Compute the model's MSE
get_mse <- function(observation, prediction) {
  n <- length(observation)
  error <- observation-prediction
  error_sq <- error*error
  mse <- sum(error_sq) / n
  return(mse)
}

##########################
###     Read files     ###
##########################

setwd("~/WORKSTATION/Dissertation/Chapter1") # this line needs to be adjusted

en_words <- read.delim("./Chapter1_DataProcessed/2022_Muller_Baldey_Processed.txt", sep = "\t")


##############################
###   Compute components   ###
##############################


# Each singular's and plural's resting activation level relies on the log transformed frequency of the stem or the plural respectively
en_words$sg_rest_act <- en_words$logSTMfreq
en_words$pl_rest_act <- en_words$logPLfreq

# The activation need until activation threshold is reached, is a function of the resting activation level
en_words$sg_t_act <- 1 / (1 + en_words$sg_rest_act)
en_words$pl_t_act <- 1 / (1 + en_words$pl_rest_act)

# add index
en_words$index <- 1:nrow(en_words)


#########################
###   Fit the model   ###
#########################

# build grid for grid search
en_words$epsilon = 0
df_eval <- data.frame()
p_grid = seq(0,1,0.001)

# Calculate log-likelihood for every p
for (i in 1:length(p_grid)) {
  
  # set t_parsing (proof of concept)
  t_parsing = p_grid[i]
  en_words$t_parsing <- t_parsing
  
  # compute model-time RTs
  RTmodelRoute <- ifelse(en_words$number == 'sing', 
                         RT_lookup(en_words$sg_t_act, en_words$epsilon), 
                         RT_decompose(en_words$sg_t_act, en_words$pl_t_act, en_words$epsilon, en_words$t_parsing))
  RTmodel <- as.numeric(data.frame(strsplit(RTmodelRoute, "-"))[1,])
  route <- as.character(data.frame(strsplit(RTmodelRoute, "-"))[2,])
  
  en_words$RTmodel <- RTmodel
  en_words$route <- route
  
  # scale into real time
  en_words$RTmodel_s <- scale_RTmodel(en_words$logRT, en_words$RTmodel)
  en_words$logRT_sdown <- scale_RTmodel(en_words$RTmodel, en_words$logRT)
  t_parsing_s <- scale_RTmodel_p(en_words$RT, en_words$RTmodel, t_parsing)
  
  # compare model output to measured RT
  dual_loglik <- my_loglik(en_words$logRT_sdown, en_words$RTmodel, en_words$number, p=t_parsing)
  my_aic <- get_aic(dual_loglik)
  my_mse <- get_mse(en_words$logRT, en_words$RTmodel_s)
  results <- data.frame('p'=p_grid[i], 'p_s'=t_parsing_s, 'mse'=my_mse, 'AIC'=my_aic, 'loglik'=dual_loglik)
  df_eval <- rbind(df_eval, results)
  df_eval <- df_eval[order(df_eval$AIC, decreasing = FALSE),]
}
print(df_eval)
