b = import('ebits/base')
io = import('ebits/io')
util = import('./util')
index = import('./index')$index

#' Returns a gene expression matrix (genes x COSMIC IDs)
#'
#' @param index_type  The column of `index` to be used for column names
#' @return            The expression matrix
basal_expression = function(map_cosmic=FALSE) {
    fname = "CCLE_Expression.Arrays_2013-03-18"
    if (!util$exists("cache", fname, ext=".RData")) {
        warning("data file does not exist, creating from raw .CEL files")
        ma = import('ebits/process/microarray')
        fnames = list.files(fname, full.names=TRUE)
        expr = oligo::read.celfiles(fnames) %>%
            ma$normalize() %>%
            ma$annotate("hgnc_symbol")
        save(expr, file=util$file("cache", fname, ext=".RData"))
    } else
        expr = io$load(util$file("cache", fname, ext=".RData"))

    expr = Biobase::exprs(expr)
    colnames(expr) = sub("\\.CEL$", "", colnames(expr))

    if (map_cosmic) {
        cn = b$match(colnames(expr),
                     from = index$`Expression arrays`,
                     to = index$COSMIC)
    
        colnames(expr) = unname(cn)
        nas = is.na(cn)
        if (any(nas)) {
            warning("dropping ", sum(nas), " arrays for COSMIC mapping")
            expr = expr[,!nas]
            cn = cn[!nas]
        }
    }

    expr
}
