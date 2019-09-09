# internal functions to compute the tumor and site specific prior probability of mutation

# base function to compute all parts of the prior P(m,M|C)
.compute_priors <- function(vr,reference,mu,prior_lod,profile_path=NULL,plot_path=NULL){
  # Obtain the set of variants that have MuTect-computed odds greater than the specified threshold
  # and pass all filters.
  prior_vars <- vr %>% dplyr::filter((mutect_odds >= prior_lod) & pass_all)

  mutation_prior = .compute_empirical_prior(prior_vars, reference,profile_path)

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
  vr4

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
