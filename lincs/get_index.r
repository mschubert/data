io = import('ebits/io')

#' Returns the LINCS metadata
#'
#' @return  A data.frame containing the experimental metadata
get_index = function() {
    fname = module_file("data", "inst.info", mustWork=TRUE)
    io$read_table(fname, quote="", header=TRUE, sep="\t")
}
