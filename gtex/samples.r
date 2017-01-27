io = import('io')

#' Return the same information
#'
#' @return  A data.frame
samples = function() {
	fpath = module_file("cache", "samples.RData", mustWork=TRUE)
	io$load(fpath)
}
