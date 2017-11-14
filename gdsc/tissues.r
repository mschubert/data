cosmic = import('./cosmic')
MASTER_LIST = cosmic$MASTER_LIST

#' Returns a vector of tissues, with COSMIC IDs as names
#'
#' @param tissue        Character vector of tissues to filter for
#' @param unknown       Data type to encode unknown tissue; default: NA
#' @param drop_unknown  Remove cell lines where tissue is unknown
#' @param TCGA          Use TCGA tissue descriptors
#' @param minN          Minimum number of cell lines to include per tissue
tissues = function(tissue=NULL, unknown=NA, drop_unknown=TRUE, TCGA=TRUE, minN=2) {
    stopifnot(!drop_unknown || is.na(unknown)) # if drop_unknown, unknown needs to be NA

    if (TCGA)
        tissueVec = as.character(MASTER_LIST$Study.Abbreviation) # v16 TCGA
    else
        tissueVec = as.character(MASTER_LIST$GDSC.description_1) # v16 tissues
    names(tissueVec) = MASTER_LIST$COSMIC.ID # v16 cosmic id

    tissueVec[is.na(tissueVec)] = unknown
    tissueVec[tissueVec %in% c("unknown", "", "UNABLE TO CLASSIFY")] = unknown
    tissueVec = na.omit(tissueVec)

    if (!is.null(tissue))
        tissueVec = tissueVec[tissueVec %in% tissue]

    n = sapply(tissueVec, function(t) sum(t==tissueVec, na.rm=T))
    tissueVec[n>=minN]
}
