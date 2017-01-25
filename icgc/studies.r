config = import('./config')

#' Return the index of all specimen
#'
#' @param tcga  Only (TRUE) or no (FALSE) TCGA cohorts
#' @return  A data.frame of details of all specimen
studies= function(tcga=NULL) {
    pdir = file.path(config$raw_data, "Projects")
    studies = grep("[A-Z]+-[A-Z]+", list.files(pdir, include.dirs=TRUE), value=TRUE)

    if (is.null(tcga))
        studies
    else if (tcga)
        studies[grepl("-US$", studies)]
    else
        studies[!grepl("-US$", studies)]
}
