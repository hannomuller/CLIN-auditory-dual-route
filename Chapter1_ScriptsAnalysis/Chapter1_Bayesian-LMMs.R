#! usr/bin/R
# Chapter1_Bayesian-LMMs.R
# H. Muller
# 2019-01-18

set.seed(123) # set this seed to reproduce the reported results


#############################
###     Load packages     ###
#############################

wants <- c('rstan')
has <- wants %in% rownames(installed.packages())
if(any(!has))install.packages(wants[!has])
lapply(wants, library, character.only = TRUE)

Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


##########################
###     Read files     ###
##########################

setwd("~/WORKSTATION/Dissertation/Chapter1") # this line needs to be adjusted

en_words <- read.delim("./Chapter1_DataProcessed/2022_Muller_Baldey_Processed.txt", sep = "\t")


####################################
###   Make data-list for Rstan   ###
####################################

# Transform variables for random effect estimation
en_words$item <- factor(en_words$item)
en_words$subject <- factor(en_words$subject)
en_words$dummy_number <- ifelse(en_words$number == 'sing', -1, 1)
stanDat <- list(N=nrow(my_df), y=my_df$logRT, x1=my_df$logRTprev_c, x2=my_df$logword_duration_c, 
                x3=my_df$logtrial_c, x4=my_df$logSTMfreq_c, x5=my_df$logFRMfreq_c, 
                x6=my_df$dummy_number, subj=as.integer(my_df$subject), J=nlevels(my_df$subject), 
                item=as.integer(my_df$word), K=nlevels(my_df$word))


##########################################
###   Fit 1. Bayesian dual-route LMM   ###
##########################################

bayesian_lmm1 <- stan(file="./Chapter1_ScriptsAnalysis/Bayesian-LMM1.stan", data=stanDat, iter=20000, chains=8, control=list(max_treedepth=20, adapt_alpha=0.99))

# Compute MSE
y_pred1 <- get_posterior_mean(bayesian_lmm1, pars = "y_pred")[,"mean-all chains"]
mse1 <- sum((en_words$logRT-y_pred1)^2) / nrow(en_words)

##########################################
###   Fit 2. Bayesian dual-route LMM   ###
##########################################

bayesian_lmm2 <- stan(file="./Chapter1_ScriptsAnalysis/Bayesian-LMM2.stan", data=stanDat, iter=20000, chains=8, control=list(max_treedepth=20, adapt_alpha=0.99))

# Compute MSE
y_pred2 <- get_posterior_mean(bayesian_lmm2, pars = "y_pred")[,"mean-all chains"]
mse2 <- sum((en_words$logRT-y_pred2)^2) / nrow(en_words)
