`%>%` = magrittr::`%>%`
.io = import('ebits/io')
.filter = import('./filter')$filter
.map_id = import('./map_id')$map_id

#' Get a matrix with miRNA expression (variance-stabilization transformed)
#'
#' @param tissue   Character vector of studies to include
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param trans    "raw" or "rpm"
#' @param ...      Paramters passed to tcga$filter
#' @return         A matrix [genes x samples]
mirna_seq = function(tissue, id_type="specimen", trans="raw", ...) {
    fpath = module_file("TCGAbiolinks-downloader/mirna_seq", mustWork=TRUE)
    expr = .io$load(file.path(fpath, paste0("TCGA-", tissue, ".RData")))

    if (trans == "raw")
        keep = "^read_count_"
    else if (trans == "rpm")
        keep = "^reads_per_million_"
    else
        stop("Transformation ", sQuote(trans), " mot recognized")

    dset = data.matrix(expr[,grepl(keep, colnames(expr))])
    colnames(dset) = sub(keep, "", colnames(dset))
    rownames(dset) = expr$miRNA_ID

    dset %>%
        .filter(...) %>%
        .map_id(id_type=id_type)
}
