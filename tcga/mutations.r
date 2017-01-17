`%>%` = magrittr::`%>%`
.io = import('ebits/io')
.map_id = import('./map_id')$map_id

.load = function(...) .io$load(module_file(..., mustWork=TRUE))

#' Get a data.frame listing all mutations and types
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the simple mutations
mutations = function(id_type="specimen", ...) {
    # 19k/22k in READ and all 20k OV have no portion in barcode, assuming 'A'
    .load("cache", "mutations.RData") %>%
        dplyr::mutate(Tumor_Sample_Barcode =
                      ifelse(nchar(Tumor_Sample_Barcode) == 15,
                             paste0(Tumor_Sample_Barcode, "A"),
                             Tumor_Sample_Barcode)) %>%
        .map_id(id_type=id_type, along="Tumor_Sample_Barcode", ...)
}
