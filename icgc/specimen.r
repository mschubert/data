io = import('ebits/io')
config = import('./config')

#' Return the index of all specimen
#'
#' @return  A data.frame of details of all specimen
index = function() {
    io$load(file.path(config$cached_data, "specimen.RData"))
}
