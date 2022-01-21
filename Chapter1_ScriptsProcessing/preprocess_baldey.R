### load packages if not installed ###
# wants <- c("ggplot2", "tidyr", "lme4", "corrplot", "car", "rcompanion", "sjPlot", "languageR", "ranger", "caret")
# has <- wants %in% rownames(installed.packages())
# if(any(!has))install.packages(wants[!has])
# lapply(wants, library, character.only = TRUE)

# Load the relevant files into the environment
setwd("./")

### source function
source('./preprocessing.R')
#source('./BALDEY/cross_valid_lmer.R')

### get data
#baldey_raw <- read.delim("baldey_new.txt", sep = " ")
baldey_raw <- read.delim("./baldey_rcp.txt", sep = " ")
#rcp_compounds <- read.delim("baldey_only_rcp_1stcomp.txt", sep = " ")
rcp_compounds <- read.delim("./baldey_only_rcp_1stcomp.txt", sep = " ", header=FALSE)
colnames(rcp_compounds) <- 'rcp_short'
#crups <- read.delim("baldey_crup.txt", sep = "\t")

# Combine baldey_rcp & rcp_compounds & baldey_crups
baldey_raw <- cbind(baldey_raw, rcp_compounds)

# Preprocess the data
my_baldey <- preprocess_baldey(baldey_raw)

my_baldey$only_suffix <- ifelse(my_baldey$one_compound == 1, '-',
                                ifelse(my_baldey$one_prefix == 1, '-',
                                       ifelse(is.na(my_baldey$affix_trans) == FALSE, my_baldey$affix, '-')))
my_baldey$only_suffix_trans <- ifelse(my_baldey$one_compound == 1, '-',
                                      ifelse(my_baldey$one_prefix == 1, '-',
                                             ifelse(is.na(my_baldey$affix_trans) == FALSE, my_baldey$affix_trans, '-')))

my_baldey$only_suffix <- ifelse(my_baldey$only_suffix == 'er',
                                ifelse(my_baldey$only_suffix_trans == '@r', 'er', '-'), my_baldey$only_suffix)
my_baldey$only_suffix <- ifelse(my_baldey$only_suffix == 'ig',
                                ifelse(my_baldey$only_suffix_trans == '@x', 'ig', '-'), my_baldey$only_suffix)
my_baldey$only_suffix_trans <- ifelse(my_baldey$only_suffix == '-', '-', my_baldey$only_suffix_trans)

my_baldey$only_suffix <- ifelse(my_baldey$one_compound == 1, '-',
                                ifelse(my_baldey$one_prefix == 1, '-',
                                       ifelse(my_baldey$only_suffix != '-', my_baldey$only_suffix, 
                                              ifelse(my_baldey$inflection == '@', ifelse(my_baldey$inflection_ort == 'en', 'en', 'e'),
                                                     ifelse(my_baldey$inflection == 's', 's',
                                                            ifelse(my_baldey$inflection == 't', 't',
                                                                   ifelse(my_baldey$inflection == 't@', 'te',
                                                                          ifelse(my_baldey$inflection == 'd@', 'de', '-'))))))))
my_baldey$only_suffix_trans <- ifelse(my_baldey$one_compound == 1, '-',
                                      ifelse(my_baldey$one_prefix == 1, '-',
                                             ifelse(my_baldey$only_suffix_trans != '-', my_baldey$only_suffix_trans, 
                                                    ifelse(my_baldey$inflection == '@', ifelse(my_baldey$inflection_ort == 'en', '@', '@'),
                                                           ifelse(my_baldey$inflection == 's', 's',
                                                                  ifelse(my_baldey$inflection == 't', 't',
                                                                         ifelse(my_baldey$inflection == 't@', 't@',
                                                                                ifelse(my_baldey$inflection == 'd@', 'd@', '-'))))))))

write.table(my_baldey, "./my_baldey.txt", sep = "\t") 
