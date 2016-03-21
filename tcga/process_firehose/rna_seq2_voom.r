# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
.b = import('ebits/base')
.io = import('ebits/io')
.ar = import('ebits/array')
.p = import('../../path')
util = import('./util')

#' Regular expression for RNA seq v2 files
archive_regex = "RSEM_genes__data.Level_3.*\\.tar(\\.gz)?$"

#' Read a single RNA seq v2 file and return results
#'
#' @param fname  File name to load
#' @param ids    Gene identifiers to map to: "hgnc"[, "entrez"]
#' @param quiet  Print file name currently processing
#' @return       An expression matrix with genes x samples
file2expr = function(fname, ids="hgnc", quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = .io$read_table(fname, header=TRUE, check.names=FALSE)
    re = re[-1, re[1,] %in% c("gene_id", "raw_count")]
    mat = data.matrix(re[,-1])
    rownames(mat) = re[[1]]
    util$voom_transform(mat, ids)
}

#' Process all RNA seq v2 files with voom
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @param save   File name to save results to (NULL: return)
#' @return       Expression matrix if save is NULL
rna_seq2_voom = function(regex=archive_regex, dir=util$data_dir, save=NULL) {
    elist = util$list_files(regex, util$data_dir) %>%
        util$unpack() %>%
        util$select("rnaseqv2")
    
    expr = elist %>%
        lapply(file2expr) #%>%
#        .ar$stack(along=2, fun.aggregate=mean) #TODO: see which sample + how to handle

    names = .b$grep("gdac.broadinstitute.org_([A-Z]+)", elist)

    if (is.null(save))
        setNames(expr, names)
    else
        .io$h5save(elist, file=.p$file("tcga", save))
}
