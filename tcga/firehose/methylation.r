# read raw data from .txt.gz files
# save into R objects for quicker loading
library(dplyr)
b = import('ebits/base')
io = import('ebits/io')
util = import('./util')

#' Regular expression for Methylation files
archive_regex = "Methylation_Preprocess.Level_3.*\\.tar(\\.gz)?$"

#' Read a single RNA seq v2 file and return results
#'
#' @param fname  File name to load
#' @param ids    Gene identifiers to map to: "hgnc"[, "entrez"]
#' @param quiet  Print file name currently processing
#' @return       An expression matrix with genes x samples
file2meth = function(fname, ids="hgnc", quiet=FALSE) {
    if (!quiet)
        message(fname)

    re = io$read_table(fname, sep="\t", header=TRUE, check.names=FALSE)
    mat = data.matrix(re[-1,-1])
    rownames(mat) = re$`Hybridization REF`[-1]
    mat
}

#' Process all GISTIC CNA files
#'
#' @param fnames  List of all unpacked files
#' @param regex   Regex filter for unpacked files
#' @return        Big matrix of HGNC x TCGA barcode
process = function(fnames, regex) {
    meth = fnames %>%
        util$list_files(regex) %>%
        lapply(file2meth) %>%
        setNames(b$grep("gdac.broadinstitute.org_([A-Z-]+)", fnames)) %>%
        narray::stack(along=2)
}

if (is.null(module_name())) {
    # we exclude FFPE samples here; that drops:
    # 7 BRCA (885 non-FFPE)
    # 3 BLCA (434)
    # 3 KIRC (480)
    # 4 PRAD (549)
    # 4 UCEC (478)
    files = util$list_files(util$data_dir, archive_regex)
    files = files[!grepl("-FFPE\\.", files)]
    files = util$unpack(files)

    # maximum standard deviation
    meth = process(files, "meth\\.by_max_stddev\\.data\\.txt")
    colnames(meth) = paste0(colnames(meth), "A") # fake portion
    fname = file.path(module_file(), "../cache", "meth_stdev.gctx")
    io$save(meth, file=fname)

    # mean of beta-values
    meth = process(files, "meth\\.by_mean\\.data\\.txt")
    colnames(meth) = paste0(colnames(meth), "A") # fake portion
    fname = file.path(module_file(), "../cache", "meth_avg.gctx")
    io$save(meth, file=fname)
}
