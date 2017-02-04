.samples = import('./samples')$samples()
.omit = import('ebits/base/omit')

#' Return tissue identifiers with sample IDs as names
#'
#' @param tissue  A character vector of which tissues to return
#' @param detail  Use more detailed tissue descriptors
tissues = function(tissue=NULL, detail=FALSE, na_rm=TRUE) {
    field = "SMTS"
    if (detail)
        field = paste0(field, "D")

    re = setNames(.samples[[field]], .samples$SAMPID)
    if (is.null(tissue))
        .omit$na(re, omit=na_rm)
    else
        .ommit$na(re[re %in% tissue], omit=na_rm)
}
