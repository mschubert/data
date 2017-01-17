# read raw data from .txt.gz files
# save into R objects for quicker loading
data_dir = module_file("..", "data", "stddata__2016_01_28", mustWork=TRUE)
analyses_dir = module_file("..", "data", "analyses__2016_01_28", mustWork=TRUE)

#' Finds files that match a regular expression
#'
#' @param dir    The directory where to look for files
#' @param regex  A regular expression to match files
#' @return       A character vector of matched files
list_files = function(dir, regex) {
    list.files(dir, pattern=regex, full.names=TRUE, recursive=TRUE)
}

#' Unpacks a Firehose archive and returns the directory
#'
#' @param files       A character vector of archives to unpack
#' @param unpack_dir  The directory in which to unpack the files (default: tmp)
#' @return            The file path to the directory
unpack = function(files, unpack_dir=tempdir()) {
    #TODO: only unpack if file is not there
    for (file in files)
        untar(file, exdir = unpack_dir)
    file.path(unpack_dir, sub(".tar", "", sub(".gz", "", basename(files))))
}
