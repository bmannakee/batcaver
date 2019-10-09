# Title
**batcaver** -  BATCAVE: Calling somatic mutations with a tumor-and-site-specific prior
## Description
This package is an implementation of the **BATCAVE** algorithm described in `Mannakee et al. ...`.
The algorithm is a post-variant-calling filter for variants that incorporates a tumor-and-site-specific prior probability of mutation.

## Dependencies and Installation
**batcaver** depends on several bioconductor packages for bam file manipulation as well as functions included in the *tidyverse*.
The dependencies can be installed prior to installation of **batcaver** using the following commands. 
```
### Not run
## install dependencies
depencies <- c("dplyr","tidyr","readr","magrittr","deconstructSigs","ggplot2","tictoc","crayon","glue","devtools")
install.packages(dependcies)
source("https://bioconductor.org/biocLite.R")
biocLite(c("Rsamtools","SomaticSignatures","GenomicAlignments","VariantAnnotation"))

## install batcaver
devtools::install_github("bmannakee/batcaver",build_vignettes=TRUE)
```


## Usage
The algorithm takes as input a MuTect output file in either VCF 4.0 or call_stats format.
The package exports one main function `run_batcave()` (See function documentation) as well as the `load_variants` (see function documentation) function which is exported to facilitate trouble-shooting of the parsing of variant file formats.
Example code for `run_batcave()`
```
### not run

## input files and variables
#'  library(BSgenome.Hsapiens.UCSC.hg38) # the BSgenome reference
#'  vcf <- /path/to/MuTect/vcf
#'  reference <- BSgenome.Hsapiens.UCSC.hg38
#'  file_type <- "vcf"
#'  seq_type <- "wes"
#'  sample_name <- "TUMOR"
#'  min_vaf <- .1
#'  min_odds <- 10
#'  plot_path <- /path/to/empirical/profile/plot.pdf
#'  profile_path <- /path/to/empirical/profile/profile.tsv
#'  fr <- run_batcave(vcf = vcf, reference = reference, file_type = file_type,
#'                    seq_type = seq_type, sample_name = sample_name, min_vaf = min_vaf,
#'                    min_odds = min_odds, plot_path = plot_path, profile_path = profile_path)
#'  }
```

