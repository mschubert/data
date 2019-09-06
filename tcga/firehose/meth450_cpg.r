library(dplyr)
b = import('ebits/base')
io = import('ebits/io')
util = import('./util')

#' Regular expression for RNA seq v2 files
archive_regex = "humanmethylation450.*\\.tar(\\.gz)?$"

#' Read a single file and return results
#'
#' @param fname  File name to load
#' @param quiet  Print file name currently processing
#' @return       An expression matrix with genes x samples
file2expr = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = io$read_table(fname, header=TRUE, check.names=FALSE)
    re = re[-1, re[1,] %in% c("Composite Element REF", "Beta_value")]
    mat = data.matrix(re[,-1])
    rownames(mat) = re[[1]]
    mat
}

#' Process all files with
#'
#' @param regex  Regular expression for archive files
#' @param dir    Directory for archive dirs
#' @return       Expression matrix if save is NULL
meth450_cpg = function(regex=archive_regex, dir=util$data_dir) {
    elist = util$list_files(dir, regex) %>%
        util$unpack() %>%
        util$list_files("data\\.txt")

    meth = elist %>%
        lapply(file2expr)

    names = b$grep("gdac.broadinstitute.org_([A-Z]+)", elist)

    setNames(meth, names)
}

if (is.null(module_name())) {
    exprs = meth450_cpg()
    # need to allow overwrite here because some TCGA cohorts have same patient
    bigmat = narray::stack(exprs, along=2, allow_overwrite=TRUE)
    fname = file.path(module_file(), "../cache", "meth450_cpg.gctx")
    io$save(bigmat, file=fname)
}
