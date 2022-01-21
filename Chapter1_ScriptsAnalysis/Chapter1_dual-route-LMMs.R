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

# Caclulate Mean Squared Error (MSE)
get_mse <- function(observation, prediction) {
  n <- length(observation)
  error <- observation-prediction
  error_sq <- error*error
  mse <- sum(error_sq) / n
  return(mse)
}

# Calculate the parsing penalty in observed time
get_pp <- function(lmm) {
  estimates <- fixef(lmm)
  intercept <- estimates[['(Intercept)']]
  if ('penalty:pnlty' %in% names(estimates)) {
    pp <- estimates[['penalty:pnlty']]
    pp_rescaled <- exp(intercept+pp)-exp(intercept) 
  } else {
    pp_rescaled <- NA
  }
  return(pp_rescaled)
}


##########################
###     Read files     ###
##########################

setwd("~/WORKSTATION/Dissertation/Chapter1") # this line needs to be adjusted

en_words <- read.delim("./Chapter1_DataProcessed/2022_Muller_Baldey_Processed.txt", sep = "\t")


##############################
###   Fit dual-route-LMM   ###
##############################

# make df for storing results
en_results <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("AIC", "MSE", 'PP', "type", "decomped_per"))

# make model
two_routes <- formula('logRT ~ logSTMfreq_c:decomp + penalty:pnlty + logRTprev_c + logword_duration_c + logtrial_c + (1 | word) + (1 | subject)')
my_lookup <- formula('logRT ~ logSTMorFRMfreq_c:lookup + logRTprev_c + logword_duration_c + logtrial_c + (1 | word) + (1 | subject)')

# get data
my_data <- en_words[complete.cases(en_words[,c('logRT', 'logSTMfreq_c', 'logFRMfreq_c')]),]
my_data2 <- my_data

### full decomposition model
en_decomp <- lmer(two_routes, data = my_data2)
AICen_decomp <- AIC(en_decomp)
mse_en_decomp <- get_mse(my_data2$logRT, predict(en_decomp))
pp_en_decomp <- get_pp(en_decomp)
newrow <- c(AICen_decomp, mse_en_decomp, pp_en_decomp, 'full_decomp', 1)
en_results[1, ] <- newrow

### full look-up model
en_lookup <- lmer(my_lookup, data = my_data2)
mse_en_lookup <- get_mse(my_data2$logRT, predict(en_lookup))
AICen_lookup <- AIC(en_lookup)
en_results <- rbind(en_results, c(AICen_lookup, mse_en_lookup, NA, 'full_listing', 0))

rm(AICen_dual_add, new_aic, threshold, final_data_add, en_dual_add_final, decomped_dual_add, decomped_dual_add_per)


##### Additive search algorithm #####

### get squared residuals
my_resid_add <- resid(en_lookup)
my_data2$my_resid_add <- my_resid_add
resid_mean_add <- with(my_data2, ave(my_resid_add, word), FUN= mean)
resid_s_mean_add <- resid_mean_add*resid_mean_add
my_data2$resid_s_mean_add <- resid_s_mean_add
minres_add <- min(resid_s_mean_add)
maxres_add <- max(resid_s_mean_add)

### get predictions
my_data$pred <- predict(en_lookup)

### Get AIC
AICen_dual_add <- 99999
pp_en_dual_add <- NA

### build iterator along unique residuals
iter <- unique(my_data2[my_data2$only_suffix == 'en', 'resid_s_mean_add'])

### fit model
for (i in 1:length(iter)) {
  t <- iter[i]
  my_data2 <- my_data
  my_data2$resid_s_mean_add <- resid_s_mean_add
  my_data2$pnlty <- ifelse(my_data2$only_suffix != 'en', 0,
                           ifelse(my_data2$resid_s_mean_add > t, 1, 0))
  my_data2$lookup <- ifelse(my_data2$only_suffix != 'en', 1,
                            ifelse(my_data2$resid_s_mean_add > t, 0, 1))
  my_data2$decomp <- ifelse(my_data2$only_suffix != 'en', 0,
                            ifelse(my_data2$resid_s_mean_add > t, 1, 0))
  en_dual_add <- lmer(two_routes, data = my_data2)
  
  ### check if new prediction shorter
  my_data2$pred_new <- predict(en_dual_add)
  my_data2$pnlty <- ifelse(my_data$only_suffix == '-', 0,
                           ifelse(my_data2$pnlty == 0, 0,
                                  ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  my_data2$lookup <- ifelse(my_data$only_suffix == '-', 1,
                            ifelse(my_data2$lookup == 1, 1,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 0, 1)))
  my_data2$decomp <- ifelse(my_data$only_suffix == '-', 0,
                            ifelse(my_data2$decomp == 0, 0,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  
  en_dual_add <- lmer(two_routes, data = my_data2)
  
  new_aic <- AIC(en_dual_add)

  if (new_aic < AICen_dual_add) {
    mse_en_dual_add <- get_mse(my_data2$logRT, predict(en_dual_add))
    pp_en_dual_add <- get_pp(en_dual_add)
    AICen_dual_add <- new_aic
    threshold <- t
    final_data_add <- my_data2
    final_data_add$resid <- resid(en_dual_add)
    en_dual_add_final <- en_dual_add
  }
}
decomped_dual_add <- final_data_add[final_data_add$pnlty == 1, 'word']
decomped_dual_add_per <- length(decomped_dual_add) / (length(my_data2[my_data2$only_suffix == 'en', 'word'])/100)
en_results <- rbind(en_results, c(AICen_dual_add, mse_en_dual_add, pp_en_dual_add, 'dual_add', decomped_dual_add_per))


##### Subtractive search algorithm #####

rm(AICen_dual_rm, new_aic, threshold, final_data_rm, en_dual_rm_final, decomped_dual_rm, decomped_dual_rm_per)

### get squared residuals
my_resid_rm <- resid(en_lookup)
my_data2$my_resid_rm <- my_resid_rm
resid_mean_rm <- with(my_data2, ave(my_resid_rm, word), FUN= mean)
resid_s_mean_rm <- resid_mean_rm*resid_mean_rm
my_data2$resid_s_mean_rm <- resid_s_mean_rm
minres_rm <- min(resid_s_mean_rm)
maxres_rm <- max(resid_s_mean_rm)

### get predictions
my_data$pred <- predict(en_decomp)

### Get AIC
AICen_dual_rm <- 9999
pp_en_dual_rm <- NA

### build iterator along unique residuals
iter <- unique(my_data2[my_data2$only_suffix == 'en', 'resid_s_mean_rm'])

### fit model
for (i in 1:length(iter)) {
  t <- iter[i]
  my_data2 <- my_data
  my_data2$resid_s_mean_rm <- resid_s_mean_rm
  my_data2$pnlty <- ifelse(my_data2$only_suffix != 'en', 0,
                           ifelse(my_data2$resid_s_mean_rm > t, 0, 1))
  my_data2$lookup <- ifelse(my_data2$only_suffix != 'en', 1,
                            ifelse(my_data2$resid_s_mean_rm > t, 1, 0))
  my_data2$decomp <- ifelse(my_data2$only_suffix != 'en', 0,
                            ifelse(my_data2$resid_s_mean_rm > t, 0, 1))
  en_dual_rm <- lmer(two_routes, data = my_data2)
  
  # check if new prediction shorter
  my_data2$pred_new <- predict(en_dual_rm)
  my_data2$pnlty <- ifelse(my_data$only_suffix == '-', 0,
                           ifelse(my_data2$pnlty == 0, 0,
                                  ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  my_data2$lookup <- ifelse(my_data$only_suffix == '-', 1,
                            ifelse(my_data2$lookup == 1, 1,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 0, 1)))
  my_data2$decomp <- ifelse(my_data$only_suffix == '-', 0,
                            ifelse(my_data2$decomp == 0, 0,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  
  en_dual_rn <- lmer(two_routes, data = my_data2)
  
  new_aic <- AIC(en_dual_rm)

  if (new_aic < AICen_dual_rm) {
    mse_en_dual_rm <- get_mse(my_data2$logRT, predict(en_dual_rm))
    pp_en_dual_rm <- get_pp(en_dual_rm)
    AICen_dual_rm <- new_aic
    threshold <- t
    final_data_rm <- my_data2
    en_dual_rm_final <- en_dual_rm
  }
}
decomped_dual_rm <- final_data_rm[final_data_rm$pnlty == 1, 'word']
decomped_dual_rm_per <- length(decomped_dual_rm) / (length(my_data2[my_data2$only_suffix == 'en', 'word'])/100)
en_results <- rbind(en_results, c(AICen_dual_rm, mse_en_dual_rm, pp_en_dual_rm, 'dual_rm', decomped_dual_rm_per))

### Process dataframe for inspection
en_results <- unique(en_results)
en_results$AIC <- round(as.numeric(en_results$AIC), 3)
en_results$MSE <- round(as.numeric(en_results$MSE), 7)
en_results$PP <- round(as.numeric(en_results$PP))
en_results$decomped_per <- round(as.numeric(en_results$decomped_per))
print(en_results)


############################
###   Smallest-AIC-LMM   ###
############################

### make df for storing results
en_results <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("AIC", "MSE", 'PP', "type", "decomped_per"))

### make model
two_routes <- formula('logRT ~ logSTMfreq_c:decomp + penalty:pnlty + logRTprev_c + logword_duration_c + logtrial_c + (1 | word) + (1 | subject)')
my_lookup <- formula('logRT ~ logSTMorFRMfreq_c:lookup + logRTprev_c + logword_duration_c + logtrial_c + (1 | word) + (1 | subject)')

### get data
my_data <- en_words[complete.cases(en_words[,c('logRT', 'logSTMfreq_c', 'logFRMfreq_c')]),]
my_data2 <- my_data

### full decomposing model
en_decomp <- lmer(two_routes, data = my_data2)
AICen_decomp <- AIC(en_decomp)
mse_en_decomp <- get_mse(my_data2$logRT, predict(en_decomp))
pp_en_decomp <- get_pp(en_decomp)
newrow <- c(AICen_decomp, mse_en_decomp, pp_en_decomp, 'full_decomp', 1)
en_results[1, ] <- newrow

### full lookup model
en_lookup <- lmer(my_lookup, data = my_data2)
mse_en_lookup <- get_mse(my_data2$logRT, predict(en_lookup))
AICen_lookup <- AIC(en_lookup)
en_results <- rbind(en_results, c(AICen_lookup, mse_en_lookup, NA, 'full_listing', 0))

rm(AICen_dual_add, new_aic, threshold, final_data_add, en_dual_add_final, decomped_dual_add, decomped_dual_add_per)


##### Additve search algorithm #####

### get squared residuals
my_resid_add <- resid(en_lookup)
my_data2$my_resid_add <- my_resid_add
resid_mean_add <- with(my_data2, ave(my_resid_add, word), FUN= mean)
resid_s_mean_add <- resid_mean_add*resid_mean_add
my_data2$resid_s_mean_add <- resid_s_mean_add
minres_add <- min(resid_s_mean_add)
maxres_add <- max(resid_s_mean_add)

### get predictions
my_data$pred <- predict(en_lookup)

### Get AIC
AICen_dual_add <- 99999
pp_en_dual_add <- NA

### build iterator along unique residuals
iter <- unique(my_data2[my_data2$only_suffix == 'en', 'resid_s_mean_add'])

### fit model
for (i in 1:length(iter)) {
  t <- iter[i]
  my_data2 <- my_data
  my_data2$resid_s_mean_add <- resid_s_mean_add
  my_data2$pnlty <- ifelse(my_data2$only_suffix != 'en', 0,
                           ifelse(my_data2$resid_s_mean_add > t, 1, 0))
  my_data2$lookup <- ifelse(my_data2$only_suffix != 'en', 1,
                            ifelse(my_data2$resid_s_mean_add > t, 0, 1))
  my_data2$decomp <- ifelse(my_data2$only_suffix != 'en', 0,
                            ifelse(my_data2$resid_s_mean_add > t, 1, 0))
  en_dual_add <- lmer(two_routes, data = my_data2)
  
  new_aic <- AIC(en_dual_add)

  if (new_aic < AICen_dual_add) {
    mse_en_dual_add <- get_mse(my_data2$logRT, predict(en_dual_add))
    pp_en_dual_add <- get_pp(en_dual_add)
    AICen_dual_add <- new_aic
    threshold <- t
    final_data_add <- my_data2
    final_data_add$resid <- resid(en_dual_add)
    en_dual_add_final <- en_dual_add
  }
}
decomped_dual_add <- final_data_add[final_data_add$pnlty == 1, 'word']
decomped_dual_add_per <- length(decomped_dual_add) / (length(my_data2[my_data2$only_suffix == 'en', 'word'])/100)
en_results <- rbind(en_results, c(AICen_dual_add, mse_en_dual_add, pp_en_dual_add, 'dual_add', decomped_dual_add_per))


##### Subtractive search algorithm #####

rm(AICen_dual_rm, new_aic, threshold, final_data_rm, en_dual_rm_final, decomped_dual_rm, decomped_dual_rm_per)

### get squared residuals
my_resid_rm <- resid(en_lookup)
my_data2$my_resid_rm <- my_resid_rm
resid_mean_rm <- with(my_data2, ave(my_resid_rm, word), FUN= mean)
resid_s_mean_rm <- resid_mean_rm*resid_mean_rm
my_data2$resid_s_mean_rm <- resid_s_mean_rm
minres_rm <- min(resid_s_mean_rm)
maxres_rm <- max(resid_s_mean_rm)

### get predictions
my_data$pred <- predict(en_decomp)

### Get AIC
AICen_dual_rm <- 9999
pp_en_dual_rm <- NA

### build iterator along unique residuals
iter <- unique(my_data2[my_data2$only_suffix == 'en', 'resid_s_mean_rm'])

### fit model
for (i in 1:length(iter)) {
  t <- iter[i]
  my_data2 <- my_data
  my_data2$resid_s_mean_rm <- resid_s_mean_rm
  my_data2$pnlty <- ifelse(my_data2$only_suffix != 'en', 0,
                           ifelse(my_data2$resid_s_mean_rm > t, 0, 1))
  my_data2$lookup <- ifelse(my_data2$only_suffix != 'en', 1,
                            ifelse(my_data2$resid_s_mean_rm > t, 1, 0))
  my_data2$decomp <- ifelse(my_data2$only_suffix != 'en', 0,
                            ifelse(my_data2$resid_s_mean_rm > t, 0, 1))
  en_dual_rm <- lmer(two_routes, data = my_data2)
  
  new_aic <- AIC(en_dual_rm)

  if (new_aic < AICen_dual_rm) {
    mse_en_dual_rm <- get_mse(my_data2$logRT, predict(en_dual_rm))
    pp_en_dual_rm <- get_pp(en_dual_rm)
    AICen_dual_rm <- new_aic
    threshold <- t
    final_data_rm <- my_data2
    en_dual_rm_final <- en_dual_rm
  }
}
decomped_dual_rm <- final_data_rm[final_data_rm$pnlty == 1, 'word']
decomped_dual_rm_per <- length(decomped_dual_rm) / (length(my_data2[my_data2$only_suffix == 'en', 'word'])/100)
en_results <- rbind(en_results, c(AICen_dual_rm, mse_en_dual_rm, pp_en_dual_rm, 'dual_rm', decomped_dual_rm_per))

en_results <- unique(en_results)
en_results$AIC <- round(as.numeric(en_results$AIC), 3)
en_results$MSE <- round(as.numeric(en_results$MSE), 7)
en_results$PP <- round(as.numeric(en_results$PP))
en_results$decomped_per <- round(as.numeric(en_results$decomped_per))
en_results


#################################
###   Independent-token LMM   ###
#################################

### make df for storing results
en_results <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("AIC", "MSE", 'PP', "type", "decomped_per"))

### make model
two_routes <- formula('logRT ~ logSTMfreq_c:decomp + penalty:pnlty + logRTprev_c + logword_duration_c + logtrial_c + (1 | word) + (1 | subject)')
my_lookup <- formula('logRT ~ logSTMorFRMfreq_c:lookup + logRTprev_c + logword_duration_c + logtrial_c + (1 | word) + (1 | subject)')

### get data
my_data <- en_words[complete.cases(en_words[,c('logRT', 'logSTMfreq_c', 'logFRMfreq_c')]),]
my_data2 <- my_data

### full decomposing model
en_decomp_rt_trial <- lmer(two_routes, data = my_data2)
AICen_decomp <- AIC(en_decomp_rt_trial)
mse_en_decomp <- get_mse(my_data2$logRT, predict(en_decomp_rt_trial))
pp_en_decomp <- get_pp(en_decomp_rt_trial)
newrow <- c(AICen_decomp, mse_en_decomp, pp_en_decomp, 'full_decomp', 1)
en_results[1, ] <- newrow

### full lookup model
en_lookup_rt_trial <- lmer(my_lookup, data = my_data2)
mse_en_lookup <- get_mse(my_data2$logRT, predict(en_lookup_rt_trial))
AICen_lookup <- AIC(en_lookup_rt_trial)
en_results <- rbind(en_results, c(AICen_lookup, mse_en_lookup, NA, 'full_listing', 0))

rm(AICen_dual_add, new_aic, threshold, final_data_add, en_dual_add_final, decomped_dual_add, decomped_dual_add_per)


##### Additive search algorithm #####

### get squared residuals
my_resid_add <- resid(en_lookup_rt_trial)
my_data2$my_resid_add <- my_resid_add^2

### get predictions
my_data$pred <- predict(en_lookup_rt_trial)

### Get AIC
AICen_dual_add <- 99999
pp_en_dual_add <- NA

### build iterator along unique residuals
iter <- unique(my_data2[my_data2$only_suffix == 'en', 'my_resid_add'])

### fit model
for (i in 1:length(iter)) {
  t <- iter[i]
  my_data2 <- my_data
  my_data2$my_resid_add <- my_resid_add
  my_data2$pnlty <- ifelse(my_data2$only_suffix != 'en', 0,
                           ifelse(my_data2$my_resid_add > t, 1, 0))
  my_data2$lookup <- ifelse(my_data2$only_suffix != 'en', 1,
                            ifelse(my_data2$my_resid_add > t, 0, 1))
  my_data2$decomp <- ifelse(my_data2$only_suffix != 'en', 0,
                            ifelse(my_data2$my_resid_add > t, 1, 0))
  en_dual_add <- lmer(two_routes, data = my_data2)
  
  # check if new prediction shorter
  my_data2$pred_new <- predict(en_dual_add)
  my_data2$pnlty <- ifelse(my_data$only_suffix == '-', 0,
                           ifelse(my_data2$pnlty == 0, 0,
                                  ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  my_data2$lookup <- ifelse(my_data$only_suffix == '-', 1,
                            ifelse(my_data2$lookup == 1, 1,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 0, 1)))
  my_data2$decomp <- ifelse(my_data$only_suffix == '-', 0,
                            ifelse(my_data2$decomp == 0, 0,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  
  en_dual_add <- lmer(two_routes, data = my_data2)
  
  new_aic <- AIC(en_dual_add)

  if (new_aic < AICen_dual_add) {
    mse_en_dual_add <- get_mse(my_data2$logRT, predict(en_dual_add))
    pp_en_dual_add <- get_pp(en_dual_add)
    AICen_dual_add <- new_aic
    threshold <- t
    final_data_add_rt_trial <- my_data2
    en_dual_add_final_rt_trial <- en_dual_add
  }
}
decomped_dual_add <- final_data_add_rt_trial[final_data_add_rt_trial$pnlty == 1, 'word']
decomped_dual_add_per <- length(decomped_dual_add) / (length(my_data2[my_data2$only_suffix == 'en', 'word'])/100)
en_results <- rbind(en_results, c(AICen_dual_add, mse_en_dual_add, pp_en_dual_add, 'dual_add', decomped_dual_add_per))


##### Subtractive search algorithm #####

rm(AICen_dual_rm, new_aic, threshold, final_data_rm, en_dual_rm_final, decomped_dual_rm, decomped_dual_rm_per)

### get squared residuals
my_resid_rm <- resid(en_lookup_rt_trial)
my_data2$my_resid_rm <- my_resid_rm^2

### get predictions
my_data$pred <- predict(en_decomp_rt_trial)

### Get AIC
AICen_dual_rm <- 9999
pp_en_dual_rm <- NA

### build iterator along unique residuals
iter <- unique(my_data2[my_data2$only_suffix == 'en', 'my_resid_rm'])

### fit model
for (i in 1:length(iter)) {
  t <- iter[i]
  my_data2 <- my_data
  my_data2$my_resid_rm <- my_resid_rm
  my_data2$pnlty <- ifelse(my_data2$only_suffix != 'en', 0,
                           ifelse(my_data2$my_resid_rm > t, 0, 1))
  my_data2$lookup <- ifelse(my_data2$only_suffix != 'en', 1,
                            ifelse(my_data2$my_resid_rm > t, 1, 0))
  my_data2$decomp <- ifelse(my_data2$only_suffix != 'en', 0,
                            ifelse(my_data2$my_resid_rm > t, 0, 1))
  en_dual_rm <- lmer(two_routes, data = my_data2)
  
  # check if new prediction shorter
  my_data2$pred_new <- predict(en_dual_rm)
  my_data2$pnlty <- ifelse(my_data$only_suffix == '-', 0,
                           ifelse(my_data2$pnlty == 0, 0,
                                  ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  my_data2$lookup <- ifelse(my_data$only_suffix == '-', 1,
                            ifelse(my_data2$lookup == 1, 1,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 0, 1)))
  my_data2$decomp <- ifelse(my_data$only_suffix == '-', 0,
                            ifelse(my_data2$decomp == 0, 0,
                                   ifelse(my_data2$pred_new <= my_data2$pred, 1, 0)))
  
  en_dual_rn <- lmer(two_routes, data = my_data2)
  
  new_aic <- AIC(en_dual_rm)

  if (new_aic < AICen_dual_rm) {
    mse_en_dual_rm <- get_mse(my_data2$logRT, predict(en_dual_rm))
    pp_en_dual_rm <- get_pp(en_dual_rm)
    AICen_dual_rm <- new_aic
    threshold <- t
    final_data_rm_rt_trial <- my_data2
    en_dual_rm_final_rt_trial <- en_dual_rm
  }
}
decomped_dual_rm <- final_data_rm_rt_trial[final_data_rm_rt_trial$pnlty == 1, 'word']
decomped_dual_rm_per <- length(decomped_dual_rm) / (length(my_data2[my_data2$only_suffix == 'en', 'word'])/100)
en_results <- rbind(en_results, c(AICen_dual_rm, mse_en_dual_rm, pp_en_dual_rm, 'dual_rm', decomped_dual_rm_per))

en_results <- unique(en_results)
en_results$AIC <- round(as.numeric(en_results$AIC), 3)
en_results$MSE <- round(as.numeric(en_results$MSE), 7)
en_results$PP <- round(as.numeric(en_results$PP))
en_results$decomped_per <- round(as.numeric(en_results$decomped_per))
en_results
