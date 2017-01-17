.io = import('ebits/io')

.load = function(...) .io$load(module_file(..., mustWork=TRUE))

#' Get a data.frame listing all clinical data available
#'
#' @param tissue   Limit selection to a single/set of tissues
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the clinical data
clinical = function(tissue=NULL, id_type=NULL) {
    re = .load("cache", "clinical.RData")
    if (!is.null(tissue))
        re = dplyr::filter(re, study %in% tissue)
    if (!is.null(id_type))
        rownames(re) = substr(re$Tumor_Sample_Barcode, 1,
                              setNames(.id_lengths, id_types)[id_type])
    re
}
