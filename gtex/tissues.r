.samples = import('./samples')$samples()

#' Return tissue identifiers with sample IDs as names
#'
#' @param tissue  Return only IDs for a specific tissue
#' @param detail  Use more detailed tissue descriptors
tissues = function(tissue, detail=FALSE) {
    field = "SMTS"
    if (detail)
        field = paste0(field, "D")

    setNames(.samples$SAMPID, .samples[[field]])
}
