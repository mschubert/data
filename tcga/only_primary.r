.barcode = import('./barcode')

#' Filters a vector of identifiers (or array) to only primary tumors
#'
#' @param obj      Character vector or array to subset
#' @param id_type  Where to subset returned IDs to (default: patient)
#' @param along    If `id` is an array, along which axes to subset
only_primary = function(obj, id_type, ...) UseMethod("only_primary")

only_primary.character = function(obj, id_type="patient") {
    if (any(sapply(obj, function(x) nchar(x) < 16)))
        stop("IDs too short for at least one case")

    keep = substr(obj, 14, 16) == "01A"
    .barcode$map(obj, id_type=id_type)
    obj[keep]
}

only_primary.array = function(obj, id_type="patient", along=2) {
    if (along == 1)
        obj[only_primary(rownames(obj), id_type=id_type),]
    else
        obj[,only_primary(colnames(obj), id_type=id_type)]
}

only_primary.list = function(obj_list, ..., simplify=FALSE, USE.NAMES=TRUE) {
    sapply(obj_list, function(x) only_primary(x, ...), simplify=simplify, USE.NAMES=USE.NAMES)
}
