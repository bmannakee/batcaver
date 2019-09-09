library(tidyverse)
library(SomaticSignatures)
library(usethis)
library(BSgenome.Hsapiens.UCSC.hg38)
# use Somatic signatures to create a table that has the fraction of the whole genome for each of the 32 possible contexts
seqs <- SomaticSignatures::kmerFrequency(BSgenome.Hsapiens.UCSC.hg38,n=1e7,k=3)
fr <- data.frame(ctxt=as.character(names(seqs)),global_prior=as.numeric(seqs))
fr <- fr %>% dplyr::mutate(context=paste0(stringr::str_sub(ctxt,1L,1L),'.',stringr::str_sub(ctxt,-1L,-1L)),
                           cref=stringr::str_sub(ctxt,2L,2L))

# There are lots of N bases in the reference
no_n_fr <- fr %>% dplyr::filter(!grepl("N",ctxt))

# Fold reference bases G becomes C, A becomes T
renamed_fr <- no_n_fr %>% mutate(cref=case_when(cref == "G" ~ "C",cref == "A" ~ "T", TRUE ~ cref))
# Now group and sum the probabilities, This is the frame that contains the global prior for the folded spectrum 32 Rows
# sums to 95% with something like 5% of the genome bases being N
global_prior_fr <- renamed_fr %>% dplyr::group_by(context,cref) %>% summarize(global_prior = sum(global_prior))

usethis::use_data(global_prior_fr, internal = TRUE)
