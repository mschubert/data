util = import('./util')
index = import('./index')$index

#' Returns a vector of tissues, with COSMIC IDs as names
#'
#' @param tissue        Character vector of tissues to filter for
#' @param unknown       Data type to encode unknown tissue; default: NA
#' @param drop_unknown  Remove cell lines where tissue is unknown
#' @param TCGA          Use TCGA tissue descriptors
#' @param minN          Minimum number of cell lines to include per tissue
tissues = function(tissue=NULL, unknown=NA, drop_unknown=TRUE, TCGA=TRUE, minN=2) {
    stopifnot(!drop_unknown || is.na(unknown)) # if drop_unknown, unknown needs to be NA

    tissueVec = setNames(index$`Site Primary`, index$`Expression arrays`)

    if (!is.null(tissue))
        tissueVec = tissueVec[tissueVec %in% tissue]

    n = sapply(tissueVec, function(t) sum(t==tissueVec, na.rm=TRUE))
    tissueVec[n>=minN]
}
