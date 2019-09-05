#' Make a function definition


run_batcave <- function(vcf, reference, file_type = "vcf", seq_type="wgs", sample_name='TUMOR', alpha = 0.05, prior_lod=5.0, plot_path=NULL, signature_path=NULL){
  if (!file.exists(vcf)){
    stop('The VCF path is not valid.')
  }

  # load the vcf as vranges
  # Rules applied are filters for SNV and sample name
  tictoc::tic("total")
  print("Beginning processing")
  vars <- load_variants(vcf,reference,sample_name, file_type)
  print("computing priors")  # compute per-site mutation probability
  # compute mu
  if (seq_type=="wes"){ N <- 3e7}
  else {N <- 3e9}
  if (is.null(mu)) {
    Nalpha <- vars %>% dplyr::filter(pass_all & ((freq >= alpha) & (freq <= 0.25))) %>% nrow()
    mu <- Nalpha/((1/alpha - 1/.25)*N)
  }
  tictoc::tic("compute priors")
  vars <- .compute_priors(vars,reference,mu,prior_lod,sig_plot_path=sig_plot_path, empirical=empirical, prior_values_path=prior_values_path, signature_path=signature_path)
  # classify. This should be separate function in the refactored package
  # need a new way to figure out if we have a wgs or exome. matters for N*mu
  # First, p1 = mu/f and p0=1-p1. mu is now a parameter for the caller

  # recalibrate for probability null given mutation rate and allele frequency.
  en <- vars %>% dplyr::filter(pass_all | tlod_only) %>% nrow() # expect they are all in the things we observe.
  print(paste0("en = ",en))
  print(en)
  vars <- vars %>% dplyr::mutate(p1 = (N*mu)/(en*(freq**2)),p0=1-p1)
  # Here hard coding the acceptable fpr as .1, so minimum threshold is 9 SEe efron and comment this better elsewhere.
  vars <- vars %>% dplyr::mutate(prior_pass = p1*prior_odds/p0)
  tictoc::toc()
  tictoc::toc()
  vars
}
