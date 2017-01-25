.io = import('ebits/io')
.bc = import('./barcode')
.map_id = import('./map_id')$map_id

.load = function(...) .io$load(module_file(..., mustWork=TRUE))

#' Get a matrix for all RNA-seq measurements
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A matrix with HGNC symbols x TCGA samples
rna_seq = function(tissue, id_type="specimen", ...) {
    library(methods) # required; otherwise h5 error
#    .load("cache", paste0(tissue, "_voom.RData")) %>%
#        .map_id(id_type=id_type, ...)
    file = h5::h5file(module_file("cache", "rna_seq2_vst.h5"), mode="r")

    barcodes = file["row"][]
    studies = .bc$barcode2study(barcodes)
    keep = studies %in% tissue

    data = file["data"][which(keep),]
    rownames(data) = barcodes[keep]
    colnames(data) = file["col"][]

    h5::h5close(file)
    .map_id(t(data), id_type=id_type, ...)
}
