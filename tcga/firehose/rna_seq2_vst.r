# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
b = import('ebits/base')
io = import('ebits/io')
ar = import('ebits/array')
util = import('./util')
rnaseq = import('ebits/process/rna-seq')

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

    re = io$read_table(fname, header=TRUE, check.names=FALSE)
    re = re[-1, re[1,] %in% c("gene_id", "raw_count")]
    mat = round(data.matrix(re[,-1]))
    rownames(mat) = re[[1]]
    rnaseq$vst(mat)
}

#' Process all RNA seq v2 files with voom
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @return       Expression matrix if save is NULL
rna_seq2_vst = function(regex=archive_regex, dir=util$data_dir) {
    elist = util$list_files(dir, regex) %>%
        util$unpack() %>%
        util$list_files("rnaseqv2")

    expr = elist %>%
        lapply(file2expr)

    names = b$grep("gdac.broadinstitute.org_([A-Z]+)", elist)

    setNames(expr, names)
}

if (is.null(module_name())) {
    exprs = rna_seq2_vst()
    # need to allow overwrite here because some TCGA cohorts have same patient
    bigmat = ar$stack(exprs, along=2, allow_overwrite=TRUE)
    fname = file.path(module_file(), "../cache", "rna_seq2_vst.gctx")
    io$save(bigmat, file=fname)
}
