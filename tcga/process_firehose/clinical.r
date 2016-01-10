# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.p = import('../path')
.b = import('../../base')
.io = import('../../io')
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

    ff = t(.io$read_table(fname, sep="\t", quote=""))
    colnames(ff) = ff[1,]
    ff = ff[-1,]
    ff = as.data.frame(ff)
    ff$study = .b$grep("gdac.broadinstitute.org_([A-Z]+)", fname)
    ff
}

#' Process all clinical files
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       Clinical data.frame if save is NULL
clinical = function(regex=archive_regex, dir=util$data_dir, save=NULL) {
    clist = util$list_files(regex, dir) %>%
        util$unpack() %>%
        util$select("clin\\.merged\\.txt")

#    clist = clist[-24] #FIXME: remove OV, invalid number of fields
    cdata = lapply(clist, file2clinical)

    cnames = do.call(.b$intersect, lapply(cdata, colnames))
    cc = do.call(rbind, lapply(cdata, function(x) x[,cnames]))
    rownames(cc) = 1:nrow(cc)

    if (is.null(save))
        cc
    else
        save(cc, file=.p$file("tcga", save))
}
