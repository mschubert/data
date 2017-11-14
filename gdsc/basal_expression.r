.io = import('ebits/io')
.file = import('./file')
cosmic = import('./cosmic')

#'
#' @param probes      Whether to return an expression values for each probe [FALSE]
#' @param duplicated  Allow two arrays for the same cell line [defunct; TRUE]
#' @return            A matrix with (genes x cell lines)
basal_expression = function(probes=FALSE, duplicated=FALSE) {
    if (probes) {
        expr_set = .io$load(module_file("cache", "expr_probes.RData", mustWork=TRUE))
        re = Biobase::exprs(expr_set)
        colnames(re) = cosmic$name2id(Biobase::pData(expr_set)$Characteristics.cell.line., warn=FALSE)
        re[,!duplicated(colnames(re))]
    } else {
        obj = .file$get('BASAL_EXPRESSION')
        rownames(obj$DATA) = obj$GENE_SYMBOLS
        obj$DATA[rownames(obj$DATA) != "",]
    }
}
