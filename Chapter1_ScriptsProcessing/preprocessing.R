preprocess_baldey <- function(baldey_rcp, opaque = TRUE, correctness = TRUE, unreasonable = TRUE, outliers = FALSE) {
  
  #### TRANSFORMATION ##################################################################
  # transform factor to character
  baldey_rcp$rcp <- as.character(baldey_rcp$rcp)
  baldey_rcp <- droplevels(baldey_rcp)
  baldey_rcp$rcp_short <- as.character(baldey_rcp$rcp_short)
  baldey_rcp$lip <- as.character(baldey_rcp$lip)
  baldey_rcp$fip <- as.character(baldey_rcp$fip)
  baldey_rcp$transcription <- as.character(baldey_rcp$transcription)
  baldey_rcp$word <- as.character(baldey_rcp$word)
  baldey_rcp$affix <- as.character(baldey_rcp$affix)
  
  # transform factor to integer
  baldey_rcp$fip.ms <- as.numeric(baldey_rcp$fip.ms)
  baldey_rcp$lip.ms <- as.numeric(baldey_rcp$lip.ms)
  baldey_rcp$RT <- as.numeric(baldey_rcp$RT)
  baldey_rcp$RTprev <- as.numeric(baldey_rcp$RTprev)
  
  
  ### CLEAN UP LIP  ##################################################################
  # there are some items in which instead of a lip a list of transcribed words is given
  # I assume that these words are the competing candidates and derive the lip from them
  # However, it would be good to check this
  baldey_rcp[baldey_rcp$word == 'hopen', 'lip'] <- 'ho' 
  baldey_rcp[baldey_rcp$word == 'zoetere', 'lip'] <- 'zu'
  baldey_rcp[baldey_rcp$word == 'lijfelijk', 'lip'] <- 'lKf'
  baldey_rcp[baldey_rcp$word == 'labieler', 'lip'] <- 'labil'
  baldey_rcp[baldey_rcp$word == 'laf', 'lip'] <- '[end]'
  baldey_rcp[baldey_rcp$word == 'forsere', 'lip'] <- 'fOrs'
  baldey_rcp[baldey_rcp$word == 'brommerig', 'lip'] <- 'brOm'
  baldey_rcp[baldey_rcp$word == 'verhaalde', 'lip'] <- 'v@rhal'
  baldey_rcp[baldey_rcp$word == 'vergrendelt', 'lip'] <- 'verG'
  baldey_rcp[baldey_rcp$word == 'lacherig', 'lip'] <- 'lAx'
  
  ### EXCLUSION  ##################################################################
  # exclude unknown words
  baldey_rcp = subset(baldey_rcp, word_unknown != "36" | word_unknown != "34" | word_unknown != "32" | word_unknown != "46")
  
  # exclude incorrect responses
  if (correctness == TRUE) {
    baldey_rcp = subset(baldey_rcp, response != "incorrect")
  }
  
  # exclude opaque words
  if (opaque == TRUE) {
    baldey_rcp <- subset(baldey_rcp, rcp != 'opaque' & word_status == 'word')
  }
  
  # exclude too fast and too responses
  if (unreasonable == TRUE) {
  baldey_rcp <- subset(baldey_rcp, RT > 150 & RTprev > 150)
  baldey_rcp <- subset(baldey_rcp, RT < 2000 & RTprev < 2000)
  }
  
  # exclude outliers
  if (outliers == TRUE) {
    baldey_rcp$logRT <- log(baldey_rcp$RT)
    grandmean = mean(baldey_rcp$logRT, na.rm = TRUE)
    SD = sd(baldey_rcp$logRT)
    upper = grandmean + 2*SD
    lower = grandmean - 2*SD
    baldey_rcp <- subset(baldey_rcp, logRT < upper & logRT > lower)   
  }
  
  
  #### RCP ##################################################################
  # improve RCP
  baldey_rcp[baldey_rcp$word == 'drukkerij', 'rcp'] <- 'dr}k'
  baldey_rcp[baldey_rcp$word == 'hijgerig', 'rcp'] <- 'hKG'
  baldey_rcp[baldey_rcp$word == 'goudachtig', 'rcp'] <- 'xMd'
  baldey_rcp[baldey_rcp$word == 'luidheid', 'rcp'] <- 'lLt'
  baldey_rcp[baldey_rcp$word == 'thuisloos', 'rcp'] <- 'tLs'
  baldey_rcp[baldey_rcp$word == 'kruidachtig', 'rcp'] <- 'krLd'
  baldey_rcp[baldey_rcp$word == 'honingwijn', 'rcp'] <- 'honINwKn'
  baldey_rcp[baldey_rcp$word == 'vispastei', 'rcp'] <- 'vIspAstK'
  baldey_rcp[baldey_rcp$word == 'moessonregen', 'rcp'] <- 'musOnreG@'
  baldey_rcp[baldey_rcp$word == 'richtinggevoel', 'rcp'] <- 'rIxtING@vul'
  baldey_rcp[baldey_rcp$word == 'visusstoornis', 'rcp_short'] <- 'viz}s'
  
  ### LIP & FIP ##################################################################
  # determine length
  baldey_rcp$rcp_length <- nchar(baldey_rcp$rcp)
  baldey_rcp$rcp_short_length <- nchar(baldey_rcp$rcp_short)
  baldey_rcp$fip_length <- ifelse(baldey_rcp$fip != '[end]', nchar(baldey_rcp$fip), nchar(baldey_rcp$transcription) + 1)
  baldey_rcp$lip_length <- ifelse(baldey_rcp$lip == '-', nchar(baldey_rcp$transcription) + 1,
                           ifelse(baldey_rcp$lip == '[end]', nchar(baldey_rcp$transcription) + 1,
                                  nchar(baldey_rcp$lip)))
  baldey_rcp$transcription_length <- nchar(baldey_rcp$transcription)
  
  # tidy up LIP & FIP
  baldey_rcp$lip[baldey_rcp$lip == '[end]'] <- baldey_rcp$transcription[baldey_rcp$lip == '[end]']
  baldey_rcp$lip[baldey_rcp$lip == '-'] <- baldey_rcp$transcription[baldey_rcp$lip == '-']
  baldey_rcp$fip[baldey_rcp$fip == '[end]'] <- baldey_rcp$transcription[baldey_rcp$fip == '[end]']
  baldey_rcp$fip[baldey_rcp$old_fip == '[end]'] <- baldey_rcp$transcription[baldey_rcp$old_fip == '[end]']
  
  
  ### RCP to UP ##################################################################
  # rcp longer, shorter or equal lip/fip?
  baldey_rcp$rcplip <- ifelse(baldey_rcp$rcp_length > baldey_rcp$lip_length, 'rcp>lip', 
                              ifelse(baldey_rcp$rcp_length < baldey_rcp$lip_length, 'rcp<lip', 
                                     'rcp=lip'))
  baldey_rcp$rcpfip <- ifelse(baldey_rcp$rcp_length > baldey_rcp$fip_length, 'rcp>fip', 
                              ifelse(baldey_rcp$rcp_length < baldey_rcp$fip_length, 'rcp<fip', 
                                     'rcp=fip'))
  baldey_rcp$rcplip_short <- ifelse(baldey_rcp$rcp_short_length > baldey_rcp$lip_length, 'rcp>lip', 
                                    ifelse(baldey_rcp$rcp_short_length < baldey_rcp$lip_length, 'rcp<lip', 
                                           'rcp=lip'))
  baldey_rcp$rcpfip_short <- ifelse(baldey_rcp$rcp_short_length > baldey_rcp$fip_length, 'rcp>fip', 
                                    ifelse(baldey_rcp$rcp_short_length < baldey_rcp$fip_length, 'rcp<fip', 
                                           'rcp=fip'))
  baldey_rcp$rcp_trans <- ifelse(baldey_rcp$rcp_length > baldey_rcp$transcription_length, 'rcp>trans', 
                                 ifelse(baldey_rcp$rcp_length < baldey_rcp$transcription_length, 'rcp<trans', 
                                        'rcp=trans'))
  baldey_rcp$rcptrans_short <- ifelse(baldey_rcp$rcp_short_length > baldey_rcp$transcription_length, 'rcp>trans', 
                                      ifelse(baldey_rcp$rcp_short_length < baldey_rcp$transcription_length, 'rcp<trans', 
                                             'rcp=trans'))
  
  baldey_rcp$rcp_by_lip <- baldey_rcp$rcp_length / baldey_rcp$lip_length
  
  # convert new condition to factor
  baldey_rcp$rcpfip <- factor(baldey_rcp$rcpfip, levels = c('rcp=fip', 'rcp<fip', 'rcp>fip'))
  baldey_rcp$rcplip <- factor(baldey_rcp$rcplip, levels = c('rcp=lip', 'rcp<lip', 'rcp>lip'))
  baldey_rcp$rcptrans <- factor(baldey_rcp$rcptrans_short, levels = c('rcp=fip', 'rcp<fip', 'rcp>fip'))
  baldey_rcp$rcpfip_short <- factor(baldey_rcp$rcpfip_short, levels = c('rcp=fip', 'rcp<fip', 'rcp>fip'))
  baldey_rcp$rcplip_short <- factor(baldey_rcp$rcplip_short, levels = c('rcp=lip', 'rcp<lip', 'rcp>lip'))
  baldey_rcp$rcptrans_short <- factor(baldey_rcp$rcptrans_short, levels = c('rcp=lip', 'rcp<lip', 'rcp>lip'))
  
  # encode distance RCP and UP as continuous variable
  baldey_rcp$rcplip_c <- baldey_rcp$lip_length - baldey_rcp$rcp_length
  baldey_rcp$rcpfip_c <- baldey_rcp$fip_length - baldey_rcp$rcp_length
  baldey_rcp$rcplip_short_c <- baldey_rcp$lip_length - baldey_rcp$rcp_short_length
  baldey_rcp$rcpfip_short_c <- baldey_rcp$fip_length - baldey_rcp$rcp_short_length
  baldey_rcp$rcptrans_c <- baldey_rcp$transcription_length - baldey_rcp$rcp_length
  baldey_rcp$rcptrans_short_c <- baldey_rcp$transcription_length - baldey_rcp$rcp_short_length
  
  # min rcp or fip/lip
  baldey_rcp$minrcplip <- transform(baldey_rcp, min = pmin(rcp_length, lip_length))[['min']]
  baldey_rcp$minrcpfip <- transform(baldey_rcp, min = pmin(rcp_length, fip_length))[['min']]
  baldey_rcp$minrcp_shortlip <- transform(baldey_rcp, min = pmin(rcp_short_length, lip_length))[['min']]
  baldey_rcp$minrcp_shortfip <- transform(baldey_rcp, min = pmin(rcp_short_length, fip_length))[['min']]
  baldey_rcp$minrcptrans <- transform(baldey_rcp, min = pmin(rcp_length, transcription_length))[['min']]
  baldey_rcp$minrcp_shorttrans <- transform(baldey_rcp, min = pmin(rcp_short_length, transcription_length))[['min']]
  
  
  ### MORPHOLOGICAL STRUCTURE   ##################################################################
  # analyze circumfixes
  circumfixed <- unique(baldey_rcp[grep('\\bge.+[td]\\b|\\bge.+en\\b', baldey_rcp$word), 'word'])
  exclude <- c('geelzucht', 'getut', 'gemakzucht', 'geiten', 'geuren', 'gehaktmolen')
  circumfixed <- circumfixed[which(!(circumfixed %in% exclude))]
  baldey_rcp$one_circumfix <- ifelse(baldey_rcp$word %in% circumfixed,1,0)
  baldey_rcp$circumfix <- ifelse(baldey_rcp$one_circumfix == 1 & endsWith(baldey_rcp$word, 't'), 'ge+t',
                             ifelse(baldey_rcp$one_circumfix == 1 & endsWith(baldey_rcp$word, 'd'), 'ge+d',
                                    ifelse(baldey_rcp$one_circumfix == 1 & endsWith(baldey_rcp$word, 'en'), 'ge+en', NA)))
  
  # analyze morphological composition
  baldey_rcp$morph_composition <- ifelse(baldey_rcp$affix == '-', 'none',
                                  ifelse(endsWith(baldey_rcp$word, baldey_rcp$affix) |
                                           endsWith(baldey_rcp$word, paste0(baldey_rcp$affix, 'e')) |
                                           endsWith(baldey_rcp$word, paste0(baldey_rcp$affix, 's')), 'suffixed',
                                  ifelse(startsWith(baldey_rcp$word, baldey_rcp$affix), 'prefixed', 'both')))
  baldey_rcp$morph_composition <- as.factor(baldey_rcp$morph_composition)
  
  # transfrom affixes into phonetic representation
  baldey_rcp$affix_trans <- ifelse(baldey_rcp$affix == "erig", '@r@x',
                            ifelse(baldey_rcp$affix == "be", 'b@',
                            ifelse(baldey_rcp$affix == "ver", 'v@r',
                            ifelse(baldey_rcp$affix == "heid", 'hKt',
                            ifelse(baldey_rcp$affix == "baar", 'bar',
                            ifelse(baldey_rcp$affix == "ig", '@x',
                            ifelse(baldey_rcp$affix == "er", '@r',
                            ifelse(baldey_rcp$affix == "schap", 'sxAp',
                            ifelse(baldey_rcp$affix == "over", 'ov@r',
                            ifelse(baldey_rcp$affix == "elijk", '@l@k',
                            ifelse(baldey_rcp$affix == "achtig", 'Axt@x',
                            ifelse(baldey_rcp$affix == "erij", '@rK',
                            ifelse(baldey_rcp$affix == "loos", 'los',
                                                                                                                          ifelse(baldey_rcp$affix == "ont", 'Ont', NA))))))))))))))
  baldey_rcp$affix_trans <- ifelse(endsWith(baldey_rcp$word, paste0(baldey_rcp$affix, 'e')) & baldey_rcp$affix == 'ig', '@G@',
                            ifelse(endsWith(baldey_rcp$word, paste0(baldey_rcp$affix, 'e')) & baldey_rcp$affix == 'er', '@r@',
                            ifelse(endsWith(baldey_rcp$word, paste0(baldey_rcp$affix, 's')), '@rs', baldey_rcp$affix_trans)))
  
  # and make the continuous variable
  baldey_rcp$prefix <- ifelse(baldey_rcp$morph_composition != 'prefixed', 0, nchar(baldey_rcp$affix_trans))
  baldey_rcp$suffix <- ifelse(baldey_rcp$morph_composition != 'suffixed', 0, nchar(baldey_rcp$affix_trans))
  
  # and get also stem_length
  baldey_rcp$stem_length <- baldey_rcp$rcp_length - baldey_rcp$prefix
  
  # determine inflectional suffixes and interfixes
  baldey_rcp$inflection <- ifelse(baldey_rcp$transcription_length == (baldey_rcp$rcp_length + baldey_rcp$suffix), "", substr(baldey_rcp$transcription, baldey_rcp$rcp_length + baldey_rcp$suffix +1, baldey_rcp$transcription_length))
  
  baldey_rcp[baldey_rcp$word == 'gevangenschap', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'zieners', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'skier', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'weddenschap', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'vleierij', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'wetenschap', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'buiig', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'zeggenschap', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'vleierig', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'vrijere', 'inflection'] <- ''
  baldey_rcp[baldey_rcp$word == 'kruier', 'inflection'] <- ''
  
  # And get the inflectional transcriptions
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  baldey_rcp$inflection_ort <- ifelse(baldey_rcp$inflection == '@' & substrRight(baldey_rcp$word, 2) == 'en', 'en',
                                  ifelse(baldey_rcp$inflection == '@r@', 'en',
                                  ifelse(baldey_rcp$inflection == 'j@', 'en',
                                  ifelse(baldey_rcp$inflection == 'j@r@', 'en',
                                  ifelse(baldey_rcp$inflection == 'w@', 'en',
                                  ifelse(baldey_rcp$inflection == 'rop@', 'en',
                                  ifelse(baldey_rcp$inflection == 'rOmp@', 'en',                                    
                                  ifelse(baldey_rcp$inflection == 'ns', 's',
                                  ifelse(baldey_rcp$inflection == 's', 's',
                                  ifelse(baldey_rcp$inflection == 't', 't',
                                  ifelse(baldey_rcp$inflection == 'mEnt', 't',
                                  ifelse(baldey_rcp$inflection == 'd@', 'de/te',
                                  ifelse(baldey_rcp$inflection == 't@', 'de/te',
                                  ifelse(baldey_rcp$inflection == '@' & substrRight(baldey_rcp$word, 1) == 'e', 'e',
                                         '-'))))))))))))))
  
  baldey_rcp$one_inflectional_suffix <- ifelse(baldey_rcp$inflection_ort == '-', 0, 1)

  baldey_rcp$inflection_length <- nchar(baldey_rcp$inflection)
  
  # Measure the length of prefix, suffix and inflection w.r.t. the length of the stem
  baldey_rcp$stemprefix <- baldey_rcp$prefix - baldey_rcp$stem_length
  baldey_rcp$stemsuffix <- baldey_rcp$suffix - baldey_rcp$stem_length
  baldey_rcp$steminffix <- baldey_rcp$inflection_length - baldey_rcp$stem_length
  
  
  ### TRANSFORMATIONS 2 ##################################################################
  # log transform RTs, word duration and frequencies
  baldey_rcp$logRT <- log(baldey_rcp$RT)
  baldey_rcp$logRTprev <- log(baldey_rcp$RTprev)
  baldey_rcp$logword_duration = log(baldey_rcp$word_duration)
  baldey_rcp$logtrial = log(baldey_rcp$trial)
  baldey_rcp$logCGN_form_freq = log(1+baldey_rcp$CGN_form_freq) # avoid 0
  baldey_rcp$logCGN_lemma_freq = log(1+baldey_rcp$CGN_lemma_freq) # avoid 0
  baldey_rcp$loginflection_length <- log(baldey_rcp$inflection_length + 1)

  # THIS IS NOW DONE USING THE FUNCTION CENTER_VARIABLES FOR EVERY SUBSET
  # # center continuous variables
  # baldey_rcp$logRT_c <- baldey_rcp$logRT - mean(baldey_rcp$logRT, na.rm = TRUE)
  # baldey_rcp$logRTprev_c <- baldey_rcp$logRTprev - mean(baldey_rcp$logRTprev, na.rm = TRUE)
  # baldey_rcp$logword_duration_c <- baldey_rcp$logword_duration - mean(baldey_rcp$logword_duration, na.rm = TRUE)
  # baldey_rcp$logCGN_form_freq_c <- baldey_rcp$logCGN_form_freq - mean(baldey_rcp$logCGN_form_freq, na.rm = TRUE)
  # baldey_rcp$trial_c <- baldey_rcp$trial - mean(baldey_rcp$trial, na.rm = TRUE)
  # baldey_rcp$rcpfip_cc <- baldey_rcp$rcpfip_c - mean(baldey_rcp$rcpfip_c, na.rm = TRUE)
  # baldey_rcp$rcplip_cc <- baldey_rcp$rcplip_c - mean(baldey_rcp$rcplip_c, na.rm = TRUE)
  # baldey_rcp$rcptrans_c <- baldey_rcp$rcptrans_c - mean(baldey_rcp$rcptrans_c, na.rm = TRUE)
  # baldey_rcp$rcpfip_short_c <- baldey_rcp$rcpfip_short_c - mean(baldey_rcp$rcpfip_short_c, na.rm = TRUE)
  # baldey_rcp$rcplip_short_c <- baldey_rcp$rcplip_short_c - mean(baldey_rcp$rcplip_short_c, na.rm = TRUE)
  # baldey_rcp$rcptrans_short_c <- baldey_rcp$rcptrans_short_c - mean(baldey_rcp$rcptrans_short_c, na.rm = TRUE)
  # baldey_rcp$fip.ms_c <- baldey_rcp$fip.ms - mean(baldey_rcp$fip.ms, na.rm = TRUE)
  # baldey_rcp$lip.ms_c <- baldey_rcp$lip.ms - mean(baldey_rcp$lip.ms, na.rm = TRUE)
  # baldey_rcp$rcp_length_c <- baldey_rcp$rcp_length - mean(baldey_rcp$rcp_length, na.rm = TRUE)
  # baldey_rcp$rcp_short_length_c <- baldey_rcp$rcp_short_length - mean(baldey_rcp$rcp_short_length, na.rm = TRUE)
  # baldey_rcp$fip_length_c <- baldey_rcp$fip_length - mean(baldey_rcp$fip_length, na.rm = TRUE)
  # baldey_rcp$lip_length_c <- baldey_rcp$lip_length - mean(baldey_rcp$lip_length, na.rm = TRUE)
  # baldey_rcp$transcription_length_c <- baldey_rcp$transcription_length - mean(baldey_rcp$transcription_length, na.rm = TRUE)
  # baldey_rcp$minrcpfip_c <- baldey_rcp$minrcpfip - mean(baldey_rcp$minrcpfip, na.rm = TRUE)
  # baldey_rcp$minrcplip_c <- baldey_rcp$minrcplip - mean(baldey_rcp$minrcplip, na.rm = TRUE)
  # baldey_rcp$minrcptrans_c <- baldey_rcp$minrcptrans - mean(baldey_rcp$minrcptrans, na.rm = TRUE)
  # baldey_rcp$minrcp_shortfip_c <- baldey_rcp$minrcp_shortfip - mean(baldey_rcp$minrcp_shortfip, na.rm = TRUE)
  # baldey_rcp$minrcp_shortlip_c <- baldey_rcp$minrcp_shortlip - mean(baldey_rcp$minrcp_shortlip, na.rm = TRUE)
  # baldey_rcp$minrcp_shorttrans_c <- baldey_rcp$minrcp_shorttrans - mean(baldey_rcp$minrcp_shorttrans, na.rm = TRUE)
  # baldey_rcp$prefix_c <- baldey_rcp$prefix - mean(baldey_rcp$prefix, na.rm = TRUE)
  # baldey_rcp$suffix_c <- baldey_rcp$suffix - mean(baldey_rcp$suffix, na.rm = TRUE)
  # baldey_rcp$inflection_length_c <- baldey_rcp$inflection_length - mean(baldey_rcp$inflection_length, na.rm = TRUE)
  # baldey_rcp$loginflection_length_c <- baldey_rcp$loginflection_length  - mean(baldey_rcp$loginflection_length , na.rm = TRUE)
  # baldey_rcp$stemprefix_c <- baldey_rcp$stemprefix - mean(baldey_rcp$stemprefix, na.rm = TRUE)
  # baldey_rcp$stemsuffix_c <- baldey_rcp$stemsuffix - mean(baldey_rcp$stemsuffix, na.rm = TRUE)
  # baldey_rcp$steminffix_c <- baldey_rcp$steminffix - mean(baldey_rcp$steminffix, na.rm = TRUE)
  # baldey_rcp$stem_length_c <- baldey_rcp$stem_length - mean(baldey_rcp$stem_length, na.rm = TRUE)
  
  # one hot encoding
  baldey_rcp$one_compound <- ifelse(baldey_rcp$morph_classification == 'compound', 1, 0)
  baldey_rcp$one_prefix <- ifelse(baldey_rcp$morph_composition == 'prefixed', 1, 0)
  baldey_rcp$one_suffix <- ifelse(baldey_rcp$morph_composition == 'suffixed', 1, 0)
  baldey_rcp$one_inflection <- ifelse(baldey_rcp$inflected == "yes", 1, 0)
  baldey_rcp$rcp_before_fip <- ifelse(baldey_rcp$rcpfip == "rcp<fip", 1, 0)
  baldey_rcp$rcp_after_fip <- ifelse(baldey_rcp$rcpfip == "rcp>fip", 1, 0)
  baldey_rcp$rcp_equal_fip <- ifelse(baldey_rcp$rcpfip == "rcp=fip", 1, 0)
  baldey_rcp$rcp_before_lip <- ifelse(baldey_rcp$rcplip == "rcp<fip", 1, 0)
  baldey_rcp$rcp_after_lip <- ifelse(baldey_rcp$rcplip == "rcp>fip", 1, 0)
  baldey_rcp$rcp_equal_lip <- ifelse(baldey_rcp$rcplip == "rcp=fip", 1, 0)
  
  return(baldey_rcp)
}

center_variables <- function(df, 
                             variables = c('rcpfip_c', 'rcplip_c', 'minrcpfip', 'minrcplip',
                                           'rcp_length', 'stem_length', 'transcription_length', 'lip_length', 'fip_length',
                                           'prefix', 'suffix', 'crup', 'rcp_by_lip', 'logsuffix_probability', 'suffix_probability', 'loglip.ms',
                                           'logfamily_size', 'logfamily_size1to42mio', 'logfamily_size10to42mio', 'logfamily_size50to42mio', 
                                           'cum_surprise_form_lip', 'cum_surprise_lemma_lip', 'suffix_entropy', 'suffix_types',
                                           'logRT', 'logRTprev', 'logword_duration', 'logCGN_form_freq', 'logCGN_lemma_freq', 'trial', 'logtrial'), 
                             addition = NULL) {
  variables = c(variables, addition)
  for (variable in variables) {
    variable_c <- paste0(variable, '_c')
    df[,variable_c] <- df[,variable] - mean(df[,variable], na.rm = TRUE)
  }
  return(df)  
}
