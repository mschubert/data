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
#' @return       An gene copy number matrix with genes x samples
file2cna = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = io$read_table(fname, sep="\t", header=TRUE, check.names=FALSE)
    mat = data.matrix(re[,-c(1:3)])
    rownames(mat) = re$`Gene Symbol`
    mat
}

#' Process all GISTIC CNA files
#'
#' @param fnames  List of all unpacked files
#' @param regex   Regex filter for unpacked files
#' @return        Big matrix of HGNC x TCGA barcode
process = function(fnames, regex) {
    cna = fnames %>%
        util$list_files(regex) %>%
        lapply(file2cna) %>%
        setNames(b$grep("gdac.broadinstitute.org_([A-Z]+)", fnames)) %>%
        narray::stack(along=2)
}

if (is.null(module_name())) {
    files = util$list_files(util$analyses_dir, archive_regex) %>%
        util$unpack()

    # raw gistic scores
    cna = process(files, "all_data_by_genes\\.txt")
    fname = file.path(module_file(), "../cache", "cna_gistic.gctx")
    io$save(cna, file=fname)

    # thresholded
    cna = process(files, "all_thresholded.by_genes\\.txt")
    cna[] = as.integer(cna)
    fname = file.path(module_file(), "../cache", "cna_gistic_thresholded.gctx")
    io$save(cna, file=fname)
}
