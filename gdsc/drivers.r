.file = import('./file')

#' Returns a list of drivers for each tissue
#'
#' @param tissue  A vector of tissues to filter for
#' @return        A list of drivers per tissue
drivers = function(tissue=NULL) {
    ig = .file$get('INTOGEN_DRIVERS')

    if (!is.null(tissue)) {
        if (tissue %in% c("COAD", "READ", "COADREAD"))
            tissue = "COREAD"

        ig = dplyr::filter(ig, Tumor_Type %in% tissue)
    }

    dplyr::transmute(ig, HGNC=ActingDriver_Symbol, tissue=Tumor_Type)
}
