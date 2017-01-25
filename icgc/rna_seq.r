idx = import('./indexing')
pu = import('./process_util')

#' Get a matrix for all RNA-seq measurements
#'
#' @param tissue   The tissue(s) to get expression for
#' @return         A matrix with HGNC symbols x TCGA samples
rna_seq = function(tissue) {
    library(methods) # required; otherwise h5 error
    file = h5::h5file(module_file("cache", "rna_seq2_vst.gctx"), mode="r")

    barcodes = file["/0/META/ROW/id"][]
    studies = .bc$barcode2study(barcodes)
    keep = studies %in% tissue

    data = file["/0/DATA/0/matrix"][which(keep),]
    rownames(data) = barcodes[keep]
    colnames(data) = file["/0/META/COL/id"][]

    h5::h5close(file)
    t(data)
}
