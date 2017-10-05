.bc = import('./barcode')
.map_id = import('./map_id')$map_id

#' Get a matrix of CpG island methylation per gene
#'
#' @param tissue   The tissue(s) to get expression for
#' @param cpg      Either max standard deviation ("stdev") or average of islands ("avg")
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param mvalues  Transform beta to m-values (default: FALSE)
#' @param genes    Genes to retrieve (default: all)
#' @return         A data.frame with data for all the simple mutations
methylation = function(tissue, id_type="specimen", cpg=c("stdev", "avg"),
                       mvalues=FALSE, genes=TRUE, ...) {
    fname = sprintf("meth_%s.gctx", match.arg(cpg))
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
    re = .map_id(t(data), id_type=id_type, ...)

    if (mvalues)
        re[] = log2(re/(1-re))
    re
}
