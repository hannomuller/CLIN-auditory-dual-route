preprocess_baldey2 <- function(my_baldey, rm_false_resp=TRUE, rm_implausible=TRUE, rm_outliers=FALSE) {
  
  # auto-correlation
  my_baldey <- my_baldey[with(my_baldey, order(subject, session, trial)), ]
  my_baldey$AR.start <- FALSE
  my_baldey$AR.start[my_baldey$trial==1] <- TRUE
  
  # included observations
  my_baldey$include <- TRUE

  # exclude incorrect response
  if (rm_false_resp == TRUE) {
    my_baldey$include[my_baldey$response == 'incorrect'] <- FALSE 
  }
  
  # exclude too fast and too late responses
  if (rm_implausible == TRUE) {
    my_baldey$include[my_baldey$RT < 150] <- FALSE
    my_baldey$include[my_baldey$RT > 2000] <- FALSE
  }
  
  # exclude outliers
  if (rm_outliers == TRUE) {
    my_baldey$logRT <- log(my_baldey$RT)
    grandmean = mean(my_baldey$logRT, na.rm = TRUE)
    SD = sd(my_baldey$logRT, na.rm=TRUE)
    upper = grandmean + 2*SD
    lower = grandmean - 2*SD
    my_baldey$include[my_baldey$logRT > upper] <- FALSE
    my_baldey$include[my_baldey$logRT < lower] <- FALSE
  }
  
  # transform some variables
  my_baldey$cob <- as.numeric(as.character(my_baldey$cob))
  my_baldey$stem_cob <- as.numeric(as.character(my_baldey$stem_cob))
  my_baldey$continuations <- as.numeric(as.character(my_baldey$continuations))
  my_baldey$morph_continuations <- as.numeric(as.character(my_baldey$morph_continuations))
  my_baldey$unrelated_conti <- my_baldey$continuations - my_baldey$morph_continuations
  my_baldey$semanticND <- as.numeric(as.character(my_baldey$semanticND))
  
  my_baldey$root_family_size <- as.numeric(as.character(my_baldey$root_family_size))
  my_baldey$root_family_size_oa <- as.numeric(as.character(my_baldey$root_family_size_oa))
  my_baldey$root_family_size_o <- as.numeric(as.character(my_baldey$root_family_size_o))
  my_baldey$root_semantic_family_size <- as.numeric(as.character(my_baldey$root_semantic_family_size))

  my_baldey$family_size_baayen <- as.numeric(as.character(my_baldey$family_size_baayen))
  my_baldey$family_size_baayen_oa <- as.numeric(as.character(my_baldey$family_size_baayen_oa))
  my_baldey$family_size_baayen_o <- as.numeric(as.character(my_baldey$family_size_baayen_o))
  my_baldey$semantic_family_size_baayen <- as.numeric(as.character(my_baldey$semantic_family_size_baayen))
  
  my_baldey$family_size_baayen_rstr <- as.numeric(as.character(my_baldey$family_size_baayen_rstr))
  my_baldey$family_size_baayen_oa_rstr <- as.numeric(as.character(my_baldey$family_size_baayen_oa_rstr))
  my_baldey$family_size_baayen_o_rstr <- as.numeric(as.character(my_baldey$family_size_baayen_o_rstr))
  my_baldey$semantic_family_size_baayen_rstr <- as.numeric(as.character(my_baldey$semantic_family_size_baayen_rstr))
  
  my_baldey$crup <- as.numeric(as.character(my_baldey$crup))
  
  # log-transform some variables
  my_baldey$logRT <- log(my_baldey$RT + 1)
  my_baldey$logmaRT <- log(my_baldey$maRT + 1)
  my_baldey$logword_duration <- log(my_baldey$word_duration + 1)
  my_baldey$logRTprev <- log(my_baldey$RTprev + 1)
  my_baldey$logtrial <- log(my_baldey$trial + 1)
  
  my_baldey$logcob <- log(my_baldey$cob + 1)
  my_baldey$logcontinuations <- log(my_baldey$continuations + 1)
  my_baldey$logmorph_continuations <- log(my_baldey$morph_continuations + 1)
  my_baldey$logunrelated_conti <- log(my_baldey$unrelated_conti + 1)
  my_baldey$logsemanticND <- log(my_baldey$semanticND + 1)
  my_baldey$logPhonND <- log(my_baldey$PhonND + 1)
  
  my_baldey$logroot_family_size <- log(my_baldey$root_family_size + 1)
  my_baldey$logroot_family_size_oa <- log(my_baldey$root_family_size_oa + 1)
  my_baldey$logroot_family_size_o <- log(my_baldey$root_family_size_o + 1)
  my_baldey$logroot_semantic_family_size <- log(my_baldey$root_semantic_family_size + 1)
  
  my_baldey$logfamily_size_baayen <- log(my_baldey$family_size_baayen + 1)
  my_baldey$logfamily_size_baayen_oa <- log(my_baldey$family_size_baayen_oa + 1)
  my_baldey$logfamily_size_baayen_o <- log(my_baldey$family_size_baayen_o + 1)
  my_baldey$logsemantic_family_size_baayen <- log(my_baldey$semantic_family_size_baayen + 1)
  
  my_baldey$logfamily_size_baayen_rstr <- log(my_baldey$family_size_baayen_rstr + 1)
  my_baldey$logfamily_size_baayen_oa_rstr <- log(my_baldey$family_size_baayen_oa_rstr + 1)
  my_baldey$logfamily_size_baayen_o_rstr <- log(my_baldey$family_size_baayen_o_rstr + 1)
  my_baldey$logsemantic_family_size_baayen_rstr <- log(my_baldey$semantic_family_size_baayen_rstr + 1)
  
  my_baldey$logcrup <- log(my_baldey$crup+1)

  return(my_baldey)
}