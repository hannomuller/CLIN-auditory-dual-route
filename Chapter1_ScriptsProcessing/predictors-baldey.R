source('./preprocessing.R')
source('./preprocess_baldey2.R')

baldey <- read.delim("baldey.txt", sep = " ")
  
# my_baldey preprocessed according to baldey_rcp3.Rmd
my_baldey <- read.delim("my_baldey_frequencies.txt", sep = "\t")
my_baldey <- subset(my_baldey, select = c(word, one_compound, one_prefix, one_suffix, one_inflection, one_circumfix, fip_length, rcp))
my_baldey <- unique(my_baldey)
baldey <- merge(baldey, my_baldey, by = 'word', all.x = TRUE)

# add information about the affixes
# add to baldey_famsize_gam.R:vvv
#my_baldey$include[my_baldey$affix=='erij' | my_baldey$affix=='erig'] <- FALSE
# ^^^

# stem frequency and family sizes --- NEWEST VERSION (IGNORE OHTER DATA BELOW)
my_baldey <- read.delim("./my_baldey_morphology.txt", sep = "\t")
#my_baldey <- subset(my_baldey, select=c(word, continuations, morph_continuations, cob, stem_cob))
my_baldey <- subset(my_baldey, select=c(word, cob, stem_cob))
my_baldey <- unique(my_baldey)
my_baldey <- merge(baldey, my_baldey, by = 'word', all.x = TRUE)

my_baldey <- my_baldey[order(my_baldey$session, my_baldey$subject, my_baldey$trial),]
row.names(my_baldey) <- 1:nrow(my_baldey)

# preprocess BALDEY
#my_baldey <- average_rt(my_baldey, 'session', 'subject')
my_baldey <- preprocess_baldey2(my_baldey)
my_baldey$include[my_baldey$word_status != 'word'] <- FALSE
nrow(my_baldey)

write.table(my_baldey, './my_baldey4.txt', sep='\t')
