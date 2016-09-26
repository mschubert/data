import('ebits/base/operators')

#' Returns a data.frame with ARACNe regulons
#'
#' This respects the principles of clean data, unlike the package
#'
#' @param tissue      TCGA cancer type, eg. "BRCA" (default: all)
#' @param data_frame  Whether to convert regulons to a data.frame
#' @param hgnc        Map Ensembl genes to HGNC symbols (default: FALSE)
#' @return            A list (or data_frame) or regulons
regulons = function(tissue=NULL, data_frame=TRUE, hgnc=FALSE) {
    data = new.env()
    avail = data(package="aracne.networks")$results[, "Item"]
    if (!is.null(tissue))
        avail = intersect(paste0("regulon", tolower(tissue)))
    for (reg in avail)
        data(list=reg, package="aracne.networks", envir=data)
    data = as.list(data)
    names(data) = toupper(sub("^regulon", "", avail))

    if (!data_frame)
        return(data)

    regulon2df = function(r) {
        tf2tab = function(x) tibble::rownames_to_column(as.data.frame(x), "target")
        regs = sapply(r, tf2tab, simplify=FALSE, USE.NAMES=TRUE) %>%
            data.table::rbindlist(idcol="regulator")
    }
    reg = sapply(data, regulon2df, simplify=FALSE, USE.NAMES=TRUE) %>%
        data.table::rbindlist(idcol="tissue")

    if (hgnc) {
        stop("not implemented")
    } else {
        reg
    }
}
