io = import('ebits/io')
config = import('./config')

#' Return the index of all specimen
#'
#' @return  A data.frame of details of all specimen
index = function() {
    fname = file.path(config$raw_data, "Summary", "specimen.all_projects.tsv.gz")
    io$read_table(fname, header=TRUE)
}
