# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.io = import('ebits/io')
.ar = import('ebits/array')
.path = import('../../path')
util = import('./util')

#' Regular expression for RNA seq v1 files
archive_regex = "mRNAseq_Preprocess\\.Level_3.*\\.tar(\\.gz)?$"

#' Read a single RNA seq v1 file and return results
#'
#' @param fname  File name to load
#' @param ids    Gene identifiers to map to: "hgnc", "entrez"
#' @param quiet  Print file name currently processing
#' @return       An expression matrix with genes x samples
file2expr = function(fname, ids="hgnc", quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = .io$read_table(fname, header=TRUE, check.names=FALSE)
    mat = data.matrix(re[,-1])
    rownames(mat) = re[,1]
    util$voom_transform(mat, ids=ids)
}

#' Process all RNA seq v1 files with voom
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       Expression matrix if save is NULL
rna_seq_voom = function(regex=archive_regex, dir=util$data_dir, save=NULL) {
    elist = util$list_files(regex, util$data_dir) %>%
        util$unpack() %>%
        util$select("[^2]\\.mRNAseq_raw_counts\\.txt") %>%
        lapply(file2expr) %>%
        .ar$stack(along=2)

    if (is.null(save))
        elist
    else
        .io$h5save(elist, file=.path$file("tcga", save))
}