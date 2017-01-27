io = import('io')

#' Return the same information
#'
#' @return  A data.frame
samples = function() {
	fpath = file.path(module_file("cache", mustWork=TRUE), "samples.RData")
    if (!file.exists(fpath)) {
        p = import('./process_data')
        p$samples()
    }
	io$load(fpath)
}
