# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.b = import('ebits/base')
.io = import('ebits/io')
.ar = import('ebits/array')
.p = import('../path')
util = import('./util')

#' Regular expression for RPPA files
archive_regex = "RPPA_AnnotateWithGene.Level_3.*\\.tar(\\.gz)?$"

#' Read a single RPPA file and return results
#'
#' @param fname  File name to load
#' @param quiet  Print file name currently processing 
#' @return       An expression matrix with genes x samples
file2rppa = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = .io$read_table(fname, header=TRUE)
    mat = data.matrix(re[,-1])
    rownames(mat) = re[,1]
    mat
}

#' Process all RPPA files
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       Protein matrix if save is NULL
rppa = function(regex=archive_regex, dir=util$data_dir, save=NULL) {
    elist = util$list_files(regex, util$data_dir) %>%
        util$unpack() %>%
        util$select("\\.rppa\\.txt")

    rppa = elist %>%
        lapply(file2rppa) %>%
#        setNames(.b$grep("gdac.broadinstitute.org_([A-Z]+)", elist)) %>%
        .ar$stack( along=2)

    if (is.null(save))
        rppa
    else
        base::save(rppa, file=.p$file("tcga", save))
}
