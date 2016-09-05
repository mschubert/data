# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
b = import('ebits/base')
io = import('ebits/io')
util = import('./util')

#' Regular expression for RPPA files
archive_regex = "[A-Z]+\\.Merge_Clinical.Level_1.*\\.tar\\.gz$"

#' Read a single clinical file and return results
#'
#' @param fname  File name to load
#' @param quiet  Print file name currently processing
#' @return       An expression matrix with genes x samples
file2clinical = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    ff = t(io$read_table(fname, sep="\t", quote=""))
    colnames(ff) = ff[1,]
    ff = ff[-1,]
    ff = as.data.frame(ff)

    if (!"admin.disease_code" %in% colnames(ff))
        stop("disease code not provided")

    ff
}

#' Process all clinical files
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @return       Clinical data.frame if save is NULL
clinical = function(regex=archive_regex, dir=util$data_dir) {
    clist = util$list_files(dir, regex) %>%
        util$unpack() %>%
        util$list_files("clin\\.merged\\.txt")

    cdata = lapply(clist, function(...) try(file2clinical(...)))
    cdata = cdata[sapply(cdata, class) != "try-error"]

    cnames = do.call(b$intersect, lapply(cdata, colnames))
    cc = do.call(rbind, lapply(cdata, function(x) x[,cnames]))
    rownames(cc) = 1:nrow(cc)
    cc
}
