`%>%` = magrittr::`%>%`
.io = import('ebits/io')
.map_id = import('./map_id')$map_id

.load = function(...) .io$load(module_file(..., mustWork=TRUE))

#' Get a matrix for all RPPA measurements
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A matrix with antibodies x TCGA samples
rppa = function(id_type="specimen", ...) {
    .load("cache", "rppa.RData") %>%
        .map_id(id_type=id_type, along=1, ...)
}
