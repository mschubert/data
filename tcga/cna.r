`%>%` = magrittr::`%>%`
.io = import('ebits/io')
.map_id = import('./map_id')$map_id

.load = function(...) .io$load(module_file(..., mustWork=TRUE))

#' Get a data.frame listing all GISTIC scores for CNAs
#'
#' @param tissue   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A data.frame with data for all the simple mutations
cna_gistic = function(tissue, id_type="specimen", thresh=FALSE, ...) {
    if (thresh)
        fname = "cna_gistic_thresholded.gctx"
    else
        fname = "cna_gistic.gctx"

    library(methods) # required; otherwise h5 error
    file = h5::h5file(module_file("cache", fname), mode="r")

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

#' Get ABSOLUTE copy numbers from synapse
#'
#' this is from synapse: https://www.synapse.org/#!Synapse:syn1710464
#'
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return  Data.frame with 'Sample' ID and genomic regions
cna_absolute = function(id_type="specimen", ...) {
    fpath = module_file("data/pancan12_absolute.segtab.txt")
    stopifnot(file.exists(fpath))
    cna = readr::read_tsv(fpath) %>%
        .map_id(id_type=id_type, along="Sample", ...)
}
