# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.b = import('ebits/base')
.io = import('ebits/io')
.df = import('data_frame')
.p = import('../../path')
util = import('./util')

#' Regular expression for CNA files
archive_regex = "[A-Z]+\\.CopyNumber_Gistic2.Level_4.*\\.tar\\.gz$"

#' Read a single clinical file and return results
#'
#' @param fname  File name to load
#' @param quiet  Print file name currently processing
#' @return       An expression matrix with genes x samples
file2cna = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = .io$read_table(fname, sep="\t", header=TRUE) %>%
        tidyr::gather(key="barcode", value="gistic", -`Gene Symbol`, -`Locus ID`, -Cytoband) %>%
        filter(gistic != 0)
}

#' Process all GISTIC CNA files
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       Clinical data.frame if save is NULL
cna_thresholded = function(regex=archive_regex, dir=util$analyses_dir, save=NULL) {
    clist = util$list_files(regex, dir) %>%
        util$unpack() %>%
        util$select("all_thresholded.by_genes\\.txt")

    cna = lapply(clist, file2cna) %>%
        setNames(.b$grep("gdac.broadinstitute.org_([A-Z]+)", clist)) %>%
        .df$add_name_col("cohort", bind=TRUE) %>%
        transmute(cohort = cohort,
                  barcode = barcode,
                  hgnc = sub("\\|.*$", "", `Gene Symbol`),
                  locus_id = `Locus ID`,
                  cytoband = Cytoband,
                  gistic = gistic)

    if (is.null(save))
        cna
    else
        save(cna, file=.p$file("tcga", save))
}