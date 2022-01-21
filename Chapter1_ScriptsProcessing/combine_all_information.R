#! usr/bin/R
# combine_all_information.R
# H. Muller
# 2019-01-18


##########################
###     Read files     ###
##########################

#setwd("~/WORKSTATION/Dissertation/Chapter1") # this line needs to be adjusted

# read the baldey file that contains the necessary frequency information and whether the word belongs to the -en plural paradigm
my_baldey_old <- read.delim("./Chapter1_ScriptsProcessing/my_baldey_frequencies.txt", sep = "\t")
my_baldey_old <- my_baldey_old[,c('word', 'CELEXform_freq', 'CELEXstem_freq', 'CELEXen_form_freq', 'en', 's')]
my_baldey_old <- unique(my_baldey_old)

# read the baldey file that contains all other information
my_baldey_new <- read.delim("./Chapter1_ScriptsProcessing/my_baldey4.txt", sep = "\t")

# merge the files
my_baldey <- merge(my_baldey_new, my_baldey_old, by='word', all.x = TRUE)
my_baldey <- my_baldey[order(my_baldey$session, my_baldey$subject, my_baldey$trial),]
row.names(my_baldey) <- 1:nrow(my_baldey)


####################################
###     Subset the right words   ###
####################################

# log transform variables
my_baldey$logCELEXform_freq <- log(my_baldey$CELEXform_freq + 1)
my_baldey$logCELEXstem_freq <- log(my_baldey$CELEXstem_freq + 1)
my_baldey$logCELEXen_form_freq <- log(my_baldey$CELEXen_form_freq + 1)
my_baldey$logRT <- log(my_baldey$RT + 1)
my_baldey$logRTprev <- log(my_baldey$RTprev + 1)
my_baldey$logword_duration <- log(my_baldey$word_duration + 1)
my_baldey$logtrial <- log(my_baldey$trial + 1)

# subset the words that belong to the -en plural paradigm
my_en <- my_baldey[complete.cases(my_baldey$en),]
my_en <- my_en[my_en$en == 1,]

# restrict subset to unambiguous nouns
my_en <- my_en[my_en$word_class == "nom",]
my_en <- my_en[my_en$number != "sing,plur",]

# remove verbs based on CELEX
verbs_celex <- c('plagen', 'lonen', 'kuilen', 'lijsten', 'schroeven', 'typen', 'bulten', 'dieven', 'snuiten', 'planten', 'polsen', 'geuren', 'buren', 'katten', 'vlaggen', 'sloten', 'vallen', 'klunzen', 'kluiven', 'loodsen', 'kruiden', 'ballen', 'tolken', 'riemen', 'treinen', 'tinten', 'peulen', 'rijen', 'buizen', 'bossen', 'dreunen', 'tonnen', 'schaduwen', 'buurten', 'vissen', 'krikken', 'vloeken', 'veren', 'middelen', 'stallen', 'stollen', 'grimassen', 'spruiten', 'hoeden')
my_en <- my_en[my_en$word %in% verbs_celex == FALSE,]

# remove words that can build their plural also with -s
my_en <- my_en[my_en$s == 0,]


#########################################
### Initiate variables for analysis   ###
#########################################

# initialize dummy variables
my_en$lookup <- ifelse(my_en$only_suffix == 'en', 0, 1)
my_en$decomp <- ifelse(my_en$only_suffix == 'en', 1, 0)
my_en$pnlty <- ifelse(my_en$only_suffix == 'en', 1, 0)

# set up penalty column
my_en$penalty <- ifelse(my_en$only_suffix == 'en', 1, 0) 

# set up frequency columns
my_en$logSTMfreq <- my_en$logCELEXstem_freq
my_en$logFRMfreq <- ifelse(my_en$only_suffix == 'en', my_en$logCELEXen_form_freq, my_en$logCELEXform_freq)
my_en$logSTMorFRMfreq <- ifelse(my_en$only_suffix == 'en', my_en$logCELEXform_freq, my_en$logCELEXstem_freq)


################################################
###   remove outliers and center variables   ###
################################################

# write function for centering variables
center_variables <- function(df, addition=NULL) {
  variables = addition
  for (variable in variables) {
    variable_c <- paste0(variable, '_c')
    df[,variable_c] <- df[,variable] - mean(df[,variable], na.rm = TRUE)
  }
  return(df)  
}

# remove outliers that are 2 SD away from mean
en_mean <- mean(my_en$RT)
en_sd <- sd(my_en$RT)
with_o <- nrow(my_en)
my_en <- my_en[my_en$RT > en_mean-2*en_sd & my_en$RT < en_mean+2*en_sd,]
without_o <- nrow(my_en)
removed <- with_o - without_o; print(removed)
removed_p <- 100-without_o / (with_o/100); print(removed_p)

# center variables
my_en <- my_en[complete.cases(my_en[,c('logRT', 'logSTMfreq', 'logFRMfreq', 'logSTMorFRMfreq')]),]
my_en <- center_variables(my_en, addition = c('logRTprev', 'logword_duration', 'logtrial', 'logSTMfreq', 'logFRMfreq', 'logSTMorFRMfreq'))

# recode en_form_freq
my_en$CELEXen_form_freq <- ifelse(my_en$number == 'sing', my_en$CELEXform_freq, my_en$CELEXen_form_freq)
my_en$logPLfreq <- log(my_en$CELEXen_form_freq + 1)

# The 1997-audio-model differentiates four groups: singular dominant singulars, plural dominant singulars, singular dominant plurals and plural dominant plurals, which are initated here
my_en$groups <- ifelse(my_en$number == 'sing',
                          ifelse(my_en$logPLfreq < my_en$logSTMfreq, 'sg_sgDom', 'sg_plDom'),
                          ifelse(my_en$logPLfreq < my_en$logSTMfreq, 'pl_sgDom', 'pl_plDom'))
my_en$groups <- as.factor(my_en$groups)


##############################
###   write data to file   ###
##############################

write.table(my_en, "Chapter1_DataProcessed/2022_Muller_Baldey_Processed.txt", sep = "\t") 
