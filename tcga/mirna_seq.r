`%>%` = magrittr::`%>%`
.io = import('ebits/io')
.filter = import('./filter')$filter
.map_id = import('./map_id')$map_id

load = function(...) .io$load(module_file(..., mustWork=TRUE))

#' Get a matrix with miRNA expression (variance-stabilization transformed)
#'
#' @param tissue   Character vector of studies to include
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param ...      Passed to map_id()
#' @return         A matrix [genes x samples]
mirna_seq = function(tissue, id_type="specimen", ...) {
    .load("cache", "mir_seq.RData") %>%
        .filter(tissue=tissue) %>%
        .map_id(id_type=id_type, ...)
}

#' Get a matrix with mature miRNA expression (log2 RPM)
#'
#' @param tissue   Character vector of studies to include
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param ...      Passed to map_id()
#' @return         A matrix [genes x samples]
mirna_seq_mature = function(tissue, id_type="specimen", ...) {
    .load("cache", "mir_seq_mature.RData") %>%
        .filter(tissue=tissue) %>%
        .map_id(id_type=id_type, ...)
}
