# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.b = import('ebits/base')
.io = import('ebits/io')
.ar = import('ebits/array')
.p = import('../../path')
util = import('./util')

#' Regular expression for PARADIGM files
archive_regex = "Pathway_Paradigm_mRNA\\.Level_4.*\\.tar(\\.gz)?$"
#archive_regex = "Pathway_Paradigm_RNASeq_And_Copy_Number\\.Level_4.*\\.tar(\\.gz)?$"
#archive_regex = "Pathway_Paradigm_RNASeq\\.Level_4.*\\.tar(\\.gz)?$"

#' Read a single PARADIGM file and return results
#'
#' @param fname  File name to load
#' @param quiet  Print file names currently processing
#' @return       An expression matrix with genes x samples
file2paradigm = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = .io$read_table(fname, header=TRUE, sep="\t")
    mat = data.matrix(re[,-1])
    rownames(mat) = re$pid_entity
    mat
}

#' Process all PARADIGM files
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       PARADIGM score matrix if save is NULL
paradigm = function(regex=archive_regex, dir=util$analyses_dir, save=NULL) {
    elist = util$list_files(regex, dir) %>%
        util$unpack() %>%
        util$select("inferredPathwayLevels\\.tab")

    elist = elist[-4] # remove GBM, duplicate w/ GBMLGG
    elist = elist[-5] # remove KIPAN, duplicate w/ other kidney
    elist = elist[-3] # remove COADREAD, overlaps w/ COAD

    paradigm = elist %>%
        lapply(file2paradigm) %>%
#        setNames(.b$grep("gdac.broadinstitute.org_([A-Z]+)", elist)) %>%
        .ar$stack( along=2)

    if (is.null(save))
        paradigm
    else
        base::save(paradigm, file=.p$file("tcga", file))
}
