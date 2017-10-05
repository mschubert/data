library(dplyr)
b = import('ebits/base')
io = import('ebits/io')
util = import('./util')

#' Regular expression for miRNA archive files
archive_regex = "miRseq_Mature_Preprocess.Level_3.*\\.tar(\\.gz)?$"

#' Read a single miRNA file
#'
#' @param fname  File name to load
#' @param quiet  Print file name currently processing
#' @return       An expression matrix with genes x samples
file2mir = function(fname, quiet=FALSE) {
    if (!quiet)
        message(fname)

    # raw reads
    re = io$read_table(fname, sep="\t", header=TRUE, check.names=FALSE)
    mat = data.matrix(re[,-1])
    rownames(mat) = sub("\\|.*$", "", re$miRseqMature)
    mat
}

#' Process all miRNA files
#'
#' @param fnames  List of all unpacked files
#' @param regex   Regex filter for unpacked files
#' @return        Big matrix of HGNC x TCGA barcode
process = function(fnames, regex) {
    mir = fnames %>%
        util$list_files(regex) %>%
        lapply(file2mir) %>%
        narray::stack(along=2)
}

if (is.null(module_name())) {
    files = util$list_files(util$data_dir, archive_regex)
    files = files[!grepl("-FFPE\\.", files)]
    files = util$unpack(files)

    mir = process(files, "miRseq_mature_RPM_log2\\.txt")
    fname = file.path(module_file(), "../cache", "mir_seq_mature.RData")
    io$save(mir, file=fname)
}
