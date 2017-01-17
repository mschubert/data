.io = import('ebits/io')
.map_id = import('./map_id')$map_id

.load = function(...) .io$load(module_file(..., mustWork=TRUE))

#' Get a data.frame listing all GISTIC scores for CNAs
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the simple mutations
cna = function(id_type="specimen", ...) {
    .load("cache", "cna.RData") %>%
        .map_id(id_type=id_type, ...)
}
