config = import('./config')

#' Return the index of all specimen
#'
#' @return  A data.frame of details of all specimen
studies= function() {
    pdir = file.path(config$raw_data, "Projects")
    list.files(pdir, include.dirs=TRUE)
}
