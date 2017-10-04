.bc = import('./barcode')
.map_id = import('./map_id')$map_id

#' Get a matrix for all RNA-seq measurements
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param genes    Genes to retrieve (default: all)
#' @return         A matrix with HGNC symbols x TCGA samples
rna_seq = function(tissue, id_type="specimen", genes=TRUE, ...) {
    library(methods) # required; otherwise h5 error
    file = h5::h5file(module_file("cache", "rna_seq2_vst.gctx"), mode="r")

    barcodes = file["/0/META/COL/id"][]
    studies = .bc$barcode2study(barcodes)
    sample_idx = which(studies %in% tissue)

    all_genes = file["/0/META/ROW/id"][]
    if (identical(genes, TRUE))
        gene_idx = seq_along(all_genes)
    else
        gene_idx = match(genes, all_genes)

    data = file["/0/DATA/0/matrix"][sample_idx, gene_idx]
    rownames(data) = barcodes[sample_idx]
    colnames(data) = all_genes[gene_idx]

    h5::h5close(file)
    .map_id(t(data), id_type=id_type, ...)
}

if (is.null(module_name())) {
    library(testthat)

    import('./rna_seq', attach=TRUE)
    expr = rna_seq("ACC", genes=c("TP53", "KRAS"))

    expect_true(all(.bc$is_barcode(colnames(expr))))
    expect_equal(rownames(expr), c("TP53", "KRAS"))
}
