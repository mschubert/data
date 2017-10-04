# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
b = import('ebits/base')
io = import('ebits/io')
df = import('data_frame')
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

    re = io$read_table(fname, sep="\t", header=TRUE) %>%
        tidyr::gather(key="barcode", value="gistic",
                      -`Gene Symbol`, -`Locus ID`, -Cytoband)
}

#' Process all GISTIC CNA files
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @return       Clinical data.frame if save is NULL
process = function(fnames) {
    cna = lapply(fnames, file2cna) %>%
        setNames(b$grep("gdac.broadinstitute.org_([A-Z]+)", .)) %>%
        df$add_name_col("cohort", bind=TRUE) %>%
        transmute(cohort = cohort,
                  barcode = barcode,
                  hgnc = sub("\\|.*$", "", `Gene Symbol`),
                  locus_id = `Locus ID`,
                  cytoband = Cytoband,
                  gistic = gistic)
}

if (is.null(module_name())) {
    files = util$list_files(util$analyses_dir, archive_regex) %>%
        util$unpack()

    # process focal amplifications
    focal = util$list_files(files, "focal_data_by_genes\\.txt") %>%
        process() %>%
        narray::construct(gistic ~ hgnc + barcode,
                          fun.aggregate = function(x) x[1])
    fname = file.path(module_file(), "../cache", "cna_focal.gctx")
    io$save(focal, file=fname)

    # process all thresholded gistic
    focal = util$list_files(files, "all_thresholded.by_genes\\.txt") %>%
        process() %>%
        narray::construct(gistic ~ hgnc + barcode,
                          fun.aggregate = function(x) x[1])
    fname = file.path(module_file(), "../cache", "cna_gistic.gctx")
    io$save(focal, file=fname)
}
