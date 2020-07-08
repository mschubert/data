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
    file = rhdf5::H5Fopen(module_file("cache", fname))

    barcodes = file$"/0/META/COL/id"
    studies = .bc$barcode2study(barcodes)
    sample_idx = which(studies %in% tissue)

    all_genes = file$"/0/META/ROW/id"
    if (identical(genes, TRUE)) {
        gene_idx = seq_along(all_genes)
    } else {
        gene_idx = match(genes, all_genes)
        no_match = is.na(gene_idx)
        if (any(no_match)) {
            warning("No data for genes: ", paste(genes[no_match], collapse=", "),
                    immediate.=TRUE)
            gene_idx = gene_idx[!no_match]
        }
    }

    data = file&"/0/DATA/0/matrix"
    # drop=F is unsupported: https://github.com/grimbough/rhdf5/issues/68
    data = matrix(data[gene_idx, sample_idx], nrow=length(gene_idx), ncol=length(sample_idx))
    colnames(data) = barcodes[sample_idx]
    rownames(data) = all_genes[gene_idx]

    rhdf5::H5Fclose(file)
    re = .map_id(data, id_type=id_type, ...)

    if (mvalues)
        re[] = log2(re/(1-re))
    re
}
