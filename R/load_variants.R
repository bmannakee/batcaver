#' Load a variant call file.
#' Allowed types are either call_stats (MuTect1 format) or MuTect2 (VCF 3 and 4) format
#'
#' The function is generally internal, but is exported to help the user in debugging input
#'
#' @param variant_file /path/to/variant/file
#' @param file_type vcf or call_stats default is vcf
#' @param reference /path/to/reference Path to the reference genome used to call variants. Must be a BSGenome object
#' @param sample_name The sample name of the tumor in the vcf file
#' @return Data frame with 5 columns: chr, start, end, ref_allele, alt_allele
#'
#' @examples
#' # Load variant calls from MuTect2
#' \dontrun{
#'
#'
#' }
load_variants <- function(vcf,file_type="vcf",reference,sample_name,){

  if (! file.exists(vcf)){
    stop('Variant file does not exist - Please check the path')
  }
  if (! (file_type %in% c("vcf","call_stats"))){
    stop("file_type must be either vcf or call_stats")
  }
  if (class(reference) != "BSgenome"){
    stop("The reference must be provided as a Bioconductor BSgenome object")
  }
  if (file_type == "vcf"){
    .load_and_filter_vcf(vcf,reference,sample_name)
  }
  else {
    .load_and_filter_cs(vcf,reference,sample_name)
  }
}

.load_and_filter_vcf <- function(vcf,reference,sample_name){
  # Load vcf as ranges, filter.
  # SNVs that pass plus SNVs that only fail tlod
  hard_filters <- S4Vectors::FilterRules(list(snv=VariantAnnotation::isSNV,
                                              sample=function(x) {VariantAnnotation::sampleNames(x)==sample_name}))
  vr <- VariantAnnotation::readVcfAsVRanges(vcf)
  vr <- S4Vectors::subsetByFilter(vr,filter=hard_filters)

  # Get rid of the GL sequences, leave just the chromosomes that appear in the variant file
  GenomeInfoDb::seqlevels(vr) <- GenomicAlignments::seqlevelsInUse(vr)
  GenomeInfoDb::genome(vr) <- GenomeInfoDb::genome(reference)[1:length(GenomeInfoDb::genome(vr))] # These need to have exactly the same name

  # Get contexts from SomaticSignatures
  vr <- SomaticSignatures::mutationContext(vr,reference)
  mcols(vr)$TLOD <- as.numeric(mcols(vr)$TLOD)

  # softFilterMatrix values are NA in many cases, this happens because they put a "." in where all filters pass in GATK4
  # This fixes it
  VariantAnnotation::softFilterMatrix(vr)[is.na(VariantAnnotation::softFilterMatrix(vr))] <- TRUE

  # Add 0-1 columns tlod_only and pass_all to identify variants that pass all heuristic filters
  # and variants that fail only the tlod filter.
  # GATK 4 changes the filter t_lod_fstar to t_lod. need different functions for mutect2 in gatk3-4
  if (sum(stringr::str_detect(colnames(VariantAnnotation::softFilterMatrix(vr)),'t_lod_fstar'))>0){
    mcols(vr)$tlod_only <- rowSums(VariantAnnotation::softFilterMatrix(vr))==ncol(VariantAnnotation::softFilterMatrix(vr))-1 & !VariantAnnotation::softFilterMatrix(vr)[,"t_lod_fstar"]
  }else{
    mcols(vr)$tlod_only <- rowSums(VariantAnnotation::softFilterMatrix(vr))==ncol(VariantAnnotation::softFilterMatrix(vr))-1 & !VariantAnnotation::softFilterMatrix(vr)[,"t_lod"]
  }
  mcols(vr)$pass_all <- rowSums(VariantAnnotation::softFilterMatrix(vr))==ncol(VariantAnnotation::softFilterMatrix(vr))

  # GATK 4 has lists of probabilities in mcols. We don't need them and they make the tibble() call fail
  mcols(vr)$SA_MAP_AF <- ""
  mcols(vr)$SA_POST_PROB <- ""

  # Format as a tibble and get cref, calt, and mutect_odds columns
  vars <- tibble::as_tibble(vr)
  vars <- vars %>% dplyr::mutate(cref=stringr::str_sub(alteration,1L,1L),
                                 calt=stringr::str_sub(alteration,2L,2L),
                                 mutect_odds=10**(TLOD - 6.0))
  vars
}


.load_and_filter_cs <- function(variant_file, reference, sample_name){
  # load the data from the call_stats file
  col_spec <- readr::cols_only(contig='c',position='i',ref_allele='c',alt_allele='c',judgement='c',tumor_f='d','t_lod_fstar'='d',failure_reasons = 'c')

  fr <- readr::read_tsv(variant_file,comment='#',col_types=col_spec) %>%dplyr::mutate(start=position,end=position) %>%
    dplyr::mutate(pass_all = judgement == "KEEP",
                  n = lengths(failure_reasons),
                  tlod_included = map(failure_reasons, str_detect, pattern = "fstar_tumor_lod"),
                  tlod_included = map_lgl(tlod_included, any),
                  tlod_only = (n == 1) & tlod_included) %>%
    dplyr::select("seqnames"=contig,position,"ref"=ref_allele,"alt"=alt_allele,'freq'=tumor_f,'TLOD'=t_lod_fstar, pass_all, tlod_only)

  # Generate VRanges object for SomaticSignatures
  vr <- VariantAnnotation::VRanges(seqnames = fr$seqnames,
                                   ranges = IRanges(start = fr$position, width = rep(1,nrow(fr))),
                                   ref = fr$ref,
                                   alt = fr$alt,
                                   TLOD = fr$TLOD,
                                   freq = fr$freq,
                                   sampleNames = rep(sample_name,nrow(fr)),
                                   pass_all = fr$pass_all,
                                   tlod_only = fr$tlod_only)
  GenomeInfoDb::seqlevels(vr) <- GenomicAlignments::seqlevelsInUse(vr)
  GenomeInfoDb::genome(vr) <- GenomeInfoDb::genome(reference)[1:length(GenomeInfoDb::genome(vr))]
  vr <- SomaticSignatures::mutationContext(vr,reference)
  mcols(vr)$TLOD <- as.numeric(mcols(vr)$TLOD)
  vars <- tibble::as_tibble(vr)
  vars <- vars %>% dplyr::mutate(cref=stringr::str_sub(alteration,1L,1L),
                                 calt=stringr::str_sub(alteration,2L,2L),
                                 mutect_odds=10**(TLOD - 6.0))

  vars
}
