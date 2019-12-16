#' Run the BATCAVE algorithm
#'
#'
#' @param vcf Output of MuTect 1.1.7 (call stats) or 2.0 (VCF 4.0).
#' @param reference The genome reference used for alignment and variant calling (BSgenome Object)
#' @param file_type The file type of "vcf" (default = "vcf")
#' @param seq_type The type of sequencing experiment. Whole Genome "wgs" or Whole Exome "wes" (default = "wgs")
#' @param sample_name The name of the tumor sample (default = "TUMOR")
#' @param tumor_fraction Tumor purity. Fraction of sample derived from the tumor (default = 1.0). Used to adjust allele frequencies
#' @param min_vaf The minimum allele frequency to use to compute mutation rate per base (default = "0.05)
#' @param min_odds The odds ratio cutoff for high-confidence mutations used to compute prior probability of mutation (default = "10" which corresponds to MuTect TLOD = 7.3)
#' @param plot_path The file name to use for a plot of the empirical mutation profile (default = "NULL" for no plot)
#' @param profile_path The file name to use for a tab separated file containing the empirical mutation profile (default = "NULL" for no file)
#' @return a data frame with 11 columns: \code{chrom} the chromosome the variant is on,
#'  \code{start} and \code{end} the position of the variant, \code{ref} and \code{alt} the reference and alternate alleles,
#'  \code{context} the tri-nucleotide context of the variant, \code{TLOD} the score reported by MuTect, \code{freq} the variant allele frequency,
#'  \code{pass_all} a boolean column that is TRUE is the variant passed all MuTect filters and FALSE otherwise,
#'  \code{tlod_only} a boolean column that is TRUE is the variant passed all MuTect filters except the TLOD filter and false otherwise,
#'  \code{pprob_variant} the posterior probability for the variant computed using the BATCAVE algorithm
#'
#'  @examples
#'  \dontrun{
#'  ## input files and varianbles
#'  library(BSgenome.Hsapiens.UCSC.hg38) # the BSgenome reference
#'  vcf <- /path/to/MuTect/vcf
#'  reference <- BSgenome.Hsapiens.UCSC.hg38
#'  file_type <- "vcf"
#'  seq_type <- "wes"
#'  sample_name <- "TUMOR"
#'  tumor_fraction <- 1.0
#'  min_vaf <- .1
#'  min_odds <- 10
#'  plot_path <- /path/to/empirical/profile/plot.pdf
#'  profile_path <- /path/to/empirical/profile/profile.tsv
#'  fr <- run_batcave(vcf = vcf, reference = reference, file_type = file_type,
#'                    seq_type = seq_type, sample_name = sample_name,
#'                    tumor_fraction = tumor_fraction, min_vaf = min_vaf, min_odds = min_odds,
#'                    plot_path = plot_path, profile_path = profile_path)
#'  }


run_batcave <- function(vcf, reference, file_type = "vcf", seq_type="wgs",
                        sample_name='TUMOR', tumor_fraction = 1.0,
                        min_vaf = 0.05, min_odds = 10.0,
                        contamination_fraction = "0.0",
                        high_conf_variant_path = "./high_condidence_variants.tsv",
                        plot_path = NULL, profile_path = NULL){
  if (!file.exists(vcf)){
    stop('The VCF path is not valid.')
  }

  tictoc::tic("Total running time")
  vars <- load_variants(vcf = vcf, reference = reference, sample_name = sample_name, file_type = file_type)
  # Adjust allele frequencies if tumor_purity of given
  stopifnot(tumor_fraction <= 1.0)
  stopifnot(tumor_fraction > 0.0)
  vars <- vars %>% dplyr::mutate(adjusted_af = freq * tumor_fraction)

  # compute per-site mutation probability
  # compute mu
  if (seq_type=="wes")
    { N <- 3e7}
  else
  {N <- 3e9}

  Nalpha <- vars %>% dplyr::filter(pass_all & ((adjusted_af >= alpha) & (adjusted_af <= 0.25))) %>% nrow()
  mu <- Nalpha/((1/alpha - 1/.25)*N)

  message(crayon::green(glue::glue("Estimated mutation rate: ",mu)))

  vars <- .compute_priors(vars,reference = reference, mu = mu, min_odds = min_odds,
                          plot_path = plot_path)

  tictoc::toc()
  vars <- vars %>% dplyr::select("chrom" = "seqnames", start, end, ref, alt, context, TLOD, freq, adjusted_af, pass_all, tlod_only, pprob_variant)

  if (!is.null(profile_path)){
    gather_context_prior_by_site(vars) %>% readr::write_tsv(profile_path)
  }

  vars
}
