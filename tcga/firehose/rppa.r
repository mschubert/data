# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
#b = import('ebits/base')
io = import('ebits/io')
ar = import('ebits/array')
util = import('./util')

#' Regular expression for RPPA files
archive_regex = "RPPA_AnnotateWithGene.Level_3.*\\.tar(\\.gz)?$"

#' Read a single RPPA file and return results
#'
#' @param fname  File name to load
#' @param quiet  Print file name currently processing 
#' @return       An expression matrix with genes x samples
file2rppa = function(fname, quiet=FALSE) {
    message(fname)

    if (!quiet)
        message(fname)

    re = io$read_table(fname, header=TRUE)
    mat = data.matrix(re[,-1])
    rownames(mat) = re[[1]]
    mat
}

#' Process all RPPA files
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @return       Protein matrix if save is NULL
rppa = function(regex=archive_regex, dir=util$data_dir) {
    elist = util$list_files(dir, regex) %>%
        util$unpack() %>%
        util$list_files("\\.rppa\\.txt")

    rppa = elist %>%
        lapply(file2rppa) %>%
        ar$stack( along=2)
}

if (is.null(module_name())) {
    rppa = rppa()
    fname = file.path(module_file(), "../cache", "rppa.RData")
    io$save(rppa, file=fname)
}
