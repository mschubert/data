.b = import('base')
.tissues = import('./tissues')
.donors = import('./donors')

#' Filters IDs based on characteristics
#'
#' @param x       Object to work on
#' @param tissue  Only return from tissue
#' @param hardy   Which kinds of deaths to include (hardy scale)
#' @return        Object with filtered IDs or logical vector with names
filter = function(x, tissue=NULL, hardy=1:4) {
    UseMethod("filter")
}

filter.character = function(x, tissue=NULL, hardy=1:4) {
    keep = rep(TRUE, length(x))

    # filter by tissue
    if (!is.null(tissue))
        keep = keep & x %in% names(.tissues$tissues(tissue))

# do not filter donors because list in data dir incomplete
    # filter by hardy death scale
#    donors = .donors$donors(hardy=hardy)
#    keep = keep & substr(x, 1, 10) %in% donors$SUBJID

    if (sum(keep) == 0)
        warning("No entries left after filtering TCGA barcodes\n", immediate.=TRUE)

    setNames(keep, x)
}

filter.matrix = function(x, tissue=NULL, hardy=1:4, along=2) {
    args = as.list(.b$match_call_defaults())[-1]

    if (along == 1) {
        args$x = rownames(x)
        x[do.call(filter, args),]
    } else if (along == 2) {
        args$x = colnames(x)
        x[,do.call(filter, args)]
    } else
        stop("only row/col implemented")
}
