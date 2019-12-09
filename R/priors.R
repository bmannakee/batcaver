# internal functions to compute the tumor and site specific prior probability of mutation

# base function to compute all parts of the prior P(m,M|C)
.compute_priors <- function(vr,reference,mu,min_odds, high_conf_variant_path, profile_path=NULL,plot_path=NULL){
  # Obtain the set of variants that have MuTect-computed odds greater than the specified threshold
  # and pass all filters.
  prior_vars <- vr %>% dplyr::filter((mutect_odds >= min_odds) & pass_all)
  prior_vars %>% write_tsv(high_confidence_variant_path)
  message(crayon::green(glue::glue("There are : ",length(prior_vars)," high confidence variants")))
  message(crayon::green(glue::glue("High confidence variant minimum allele frequency is : ",max(prior_vars$freq))))
  mutation_prior = .compute_empirical_prior(prior_vars = prior_vars, reference = reference, profile_path = profile_path, plot_path = plot_path)

  context_prior <- .compute_context_prior(prior_vars) #P(C | M) This can be zero, so need to normalize, this is where the dirichlet will come in.
  rm(prior_vars)
  gc()

  global_prior <- global_prior_fr
  #P(M) = 3e-6 = mu
  # full prior= mutation_prior * context_prior * P(M)/P(C)
  # 4/23/2018: Add the real global context prior. replacing P(C) = 1/96
  joint_prior <- mutation_prior %>%
    dplyr::left_join(global_prior,by=c('context','cref')) %>%
    dplyr::left_join(context_prior,by=c('context','cref')) %>%
    dplyr::mutate(joint_prior=((mutation_prior * context_prior)/global_prior)*mu) %>%
    dplyr::select(cref,calt,context,mutation_prior,joint_prior,global_prior)
  rm(mutation_prior)
  rm(context_prior)
  rm(global_prior)
  if (!is.null(plot_path)){
    .save_empirical_signature_plot(joint_prior, plot_path)
  }
  gc()
  vr1 <- vr %>% dplyr::left_join(joint_prior,by=c('context','cref','calt'))
  rm(vr)
  rm(joint_prior)
  gc()
  vr2 <- vr1 %>% dplyr::mutate(joint_q = 1 - joint_prior)
  vr2 <- vr2 %>% dplyr::mutate(log_joint_prior = log10(joint_prior))
  vr2 <- vr2 %>% dplyr::mutate(log_joint_q = log10(joint_q))
  rm(vr1)
  gc()
  vr3 <- vr2 %>% dplyr::mutate(logit_prior=log_joint_prior - log_joint_q)
  rm(vr2)
  gc()
  vr4 <- vr3 %>% dplyr::mutate(prior_odds=10**(TLOD + logit_prior))
  vr4 <- vr4 %>% dplyr::mutate(pprob_variant = prior_odds/(1 + prior_odds))

}


.compute_empirical_prior <- function(prior_vars, reference, profile_path, plot_path){

  fr <- prior_vars %>% dplyr::select("Sample"="sampleNames","chr"="seqnames","pos"="start",ref,alt)
  # Workaround to get mut.to.sigs.input to work when using "BSgenome.Hsapiens.1000genomes.hs37d5"
  # It converts everything to work with UCSC.hg19 coordinates, and while this works with NCBI.GRCh38, it fails with 1000genomes.hs37d5
  # This solution allows deconstructSigs to use its default hg19 reference when the reference is 1000genomes.hs37d5
  if (reference@pkgname == "BSgenome.Hsapiens.1000genomes.hs37d5") {
    ds_input <- deconstructSigs::mut.to.sigs.input(as.data.frame(fr),bsg = NULL)
  } else {
    ds_input <- deconstructSigs::mut.to.sigs.input(as.data.frame(fr),bsg = reference)
  }

  rm(fr) # Clean up memory
  gc()
  ds_input <- ds_input %>% t() %>% as_tibble(rownames = "ctxt")
  ds_input <- ds_input %>% mutate(TUMOR = TUMOR + 1) # This is generatating the dirichlet distribution with concentration parameter 1.

  mut_prior <- ds_input %>% dplyr::mutate(context=paste0(stringr::str_sub(ctxt,1L,1L),'.',stringr::str_sub(ctxt,-1L,-1L)),
                                          cref=paste0(stringr::str_sub(ctxt,3L,3L)),
                                          calt=paste0(stringr::str_sub(ctxt,5L,5L))) %>%
    dplyr::select(context,cref, calt, "mutation_prior" = "TUMOR") %>%
    dplyr::mutate(mutation_prior = mutation_prior/sum(mutation_prior))
  if (!is.null(profile_path)){
    mut_prior %>% readr::write_tsv(profile_path)
  }

  if (!is.null(plot_path)){
  .save_empirical_signature_plot(mut_prior, plot_path)
  }

  rm(ds_input) # Clean up memory
  gc()
  mut_prior <- mut_prior %>%
    group_by(context,cref) %>%
    mutate(mutation_prior=mutation_prior/sum(mutation_prior)) %>%
    ungroup()
  if (!is.null(prior_values_path)){
    mut_prior %>% readr::write_tsv(prior_values_path)
  }
  mut_prior
}

.compute_mutation_prior <- function(vars) {
  vars %>% dplyr::group_by(context,cref,calt) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(mutation_prior=n/sum(n)) %>%
    dplyr::select(-n)
}

.compute_context_prior <- function(vars){
  in_fr <- vars %>% dplyr::group_by(context,cref)
  in_fr <- in_fr %>%dplyr::summarise(n=n()) %>% dplyr::ungroup()

  # left_join to global_prior_fr to ensure that all contexts are represented
  # get rid of the extra column from global prior.
  # Anything not in context prior will now be present with n value = NA
  joined_fr <- global_prior_fr %>%
    left_join(in_fr,by=c("context","cref")) %>%
    dplyr::select(-global_prior)
  joined_fr <- joined_fr %>% replace_na(list(n=1)) # Make the count 1 for missing contexts
  # Now get proportions
  out_fr <- joined_fr %>% dplyr::mutate(context_prior=n/sum(n)) %>%
    dplyr::select(-n)
  out_fr
}

.save_empirical_signature_plot <- function(contexts, file_path){
  fr <- contexts %>% dplyr::mutate(alteration = paste0(cref, ">", calt))
  p <- ggplot(fr, aes(x=context, y=mutation_prior, fill=alteration)) +
    geom_bar(stat = "identity", position = "identity") +
    facet_grid( ~ alteration) + theme_classic() +
    theme(axis.text.x = element_text(size = 6, face = "bold", angle = 90),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 12),
          strip.text.x = element_text(size = 10, face = "bold")) + xlab("Context") + ylab("Proportion") + guides(fill = F)
  ggsave(filename = file_path, plot = p, height = 4, width = 7.5, units = "in")
}

gather_context_prior_by_site <- function(fr){
  # input data frame (fr) is the output from smartcallr
  # Create the data frame needed to compute the KL divergence
  # between prior and target as each mutation is added in
  # descending order of allele frequency
  # This function is using the development version of tidyr to
  # get the pivot_wider and pivot_longer functions which make this much much easier
  # The end result is the signature prior at every new mutation evaluated, order by descinding vaf

  # Get the 96 mutation categories
  all_96_categories <- get_signature(1) %>% pull(context)

  # Get only the mutations we want
  ppt <- fr %>% dplyr::filter(pass_all) %>% # only variants that pass
    tidyr::separate(context, into=c("left","right")) %>% # make substitution_type that matches COSMIC key A.T -> left=A,right=T
    dplyr::mutate(mutation_type = paste0(left,"[",cref,">",calt,"]",right)) %>% # A[C>A]T
    dplyr::select(seqnames,start,mutation_type,TLOD) %>%
    dplyr::mutate(mutation_type = factor(mutation_type, levels = all_96_categories)) %>%
    unique() %>% # Try making this a factor before spreading to see if we get all 96 columns.
    dplyr::arrange(desc(TLOD))





  # now spread it into a one-hot matrix
  ppt <- ppt %>% mutate(seen = 1) %>%
    tidyr::pivot_wider(names_from = mutation_type,
                       values_from = seen,
                       values_fill = list(seen=0))

  # The tumor might not have a mutation in each of the 96 contexts, but we need them all.
  these_names <- colnames(ppt %>% dplyr::select(contains(">")))
  missing <- setdiff(all_96_categories,these_names)
  ppt[missing] <- 0

  # now get the cumulative sum as each mutation_type appears
  ppt <- ppt %>% dplyr::mutate_at(.funs = funs("cum" = cumsum(.) + 1),
                                  .vars = vars(contains(">"))) %>% # adds a new column with "_cum" that is the cumsum +1 for each. This is now the dirichlet prior!
    dplyr::select(seqnames,start,TLOD,contains("_cum")) %>% # get rid of the non-cumsum cols
    purrr::set_names(~str_replace_all(.,"_cum",""))# rename the columns so they are normal again

  # Now make the dirichlet prior into expected value of the probability vector for the multinomial.
  ppt <- ppt %>% dplyr::mutate_(row_total = ~Reduce(`+`, .[4:ncol(.)])) %>%
    mutate_at(.funs = funs(prop = ./row_total),
              .vars = vars(contains(">"))) %>%
    dplyr::select(seqnames,start,TLOD,contains("_")) %>%
    purrr::set_names(~str_replace_all(.,"_prop","")) %>%
    dplyr::select(-row_total)

  # Here you go!
  ppt
}

get_signature <- function(...){
  dots <- list(...)
  dots
  fr <- suppressMessages(read_tsv('./data/signatures_probabilities.txt') %>% dplyr::select(-contains("X")))

  names(fr) %<>% stringr::str_replace_all("\\s","_") %>% tolower
  props <- fr %>% dplyr::select(paste0("signature_",dots)) %>%
    dplyr::mutate_(row_total = ~Reduce(`+`,.)) %>%
    dplyr::mutate(row_prop = row_total/sum(row_total)) %>% pull(row_prop)
  out <- tibble(context = fr$somatic_mutation_type,prop = props)
  out
}
