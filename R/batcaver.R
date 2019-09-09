#' Make a function definition


run_batcave <- function(vcf, reference, file_type = "vcf", seq_type="wgs",
                        sample_name='TUMOR', min_vaf = 0.05, min_odds = 10.0,
                        plot_path = NULL, signature_path = NULL){
  if (!file.exists(vcf)){
    stop('The VCF path is not valid.')
  }

  tictoc::tic("Run time")
  vars <- load_variants(vcf,reference,sample_name, file_type)
  print("computing priors")  # compute per-site mutation probability
  # compute mu
  if (seq_type=="wes"){ N <- 3e7}
  else {N <- 3e9}
  if (is.null(mu)) {
    Nalpha <- vars %>% dplyr::filter(pass_all & ((freq >= alpha) & (freq <= 0.25))) %>% nrow()
    mu <- Nalpha/((1/alpha - 1/.25)*N)
  }

  vars <- .compute_priors(vars,reference = reference,mu = mu,min_odds = min_odds,
                          sig_plot_path = sig_plot_path, empirical = empirical,
                          prior_values_path = prior_values_path, signature_path = signature_path)

  tictoc::toc()
  vars
}
