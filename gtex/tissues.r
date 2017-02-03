.samples = import('./samples')$samples()

#' Return tissue identifiers with sample IDs as names
#'
#' @param tissue  A character vector of which tissues to return
#' @param detail  Use more detailed tissue descriptors
tissues = function(tissue=NULL, detail=FALSE) {
    field = "SMTS"
    if (detail)
        field = paste0(field, "D")

    re = setNames(samples[[field]], .samples$SAMPID)
    if (is.null(tissue))
        re
    else
        re[re %in% tissue]
}
