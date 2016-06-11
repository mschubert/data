.barcode = import('./barcode')

#' List all available ID types
id_types = c("patient", "specimen", "full")

#' Lengths of the respective ID types
id_lengths = setNames(c(12, 16, 100), id_types)

#' Types of data subsets
subset_types = c("primary")

#' IDs of data subsets (part of the barcode)
subset_ids = setNames(c("01A"), subset_types)

#' Maps barcodes to a specific identifier type
map_id = function(obj, id_type, subset, drop, ...) UseMethod("map_id")

map_id.character = function(obj, id_type=NULL, subset=NULL, drop=TRUE, ...) {
    if (!is.null(subset) && ! subset %in% subset_types)
        stop("Invalid subset. Choose one of: ", paste(subset_types, sep=", "))

    if (any(sapply(obj, function(x) nchar(x) < 16)))
        stop("IDs too short for at least one case")

    if (!is.null(subset))
        keep = substr(obj, 14, 16) == subset_ids[subset]
    else
        keep = TRUE

    if (!is.null(id_type))
        obj = substr(obj, 1, id_lengths[id_type])

    if (drop)
        obj[keep]
    else {
        obj[!keep] = NA
        obj
    }
}

map_id.matrix = function(obj, id_type=NULL, along=2, subset=NULL, drop=TRUE) {
    ids = map_id(dimnames(obj)[[along]], id_type=id_type, subset=subset, drop=FALSE)
    keep = !is.na(ids)

    if (along == 1) {
        rownames(obj) = ids
        if (drop)
            obj[keep,]
        else
            obj
    } else {
        colnames(obj) = ids
        if (drop)
            obj[,keep]
        else
            obj
    }
}

map_id.data.frame = function(obj, along, id_type=NULL, subset=NULL, drop=TRUE) {
    if (is.numeric(along))
        map_id.matrix(obj, id_type=id_type, subset=subset, along=along, drop=drop)
    else {
        ids = map_id(obj[[along]], id_type=id_type, subset=subset, drop=FALSE)
        keep = !is.na(ids)

        obj[[along]] = ids
        if (drop)
            obj[keep,]
        else
            obj
    }
}

map_id.list = function(obj_list, ..., simplify=FALSE, USE.NAMES=TRUE) {
    sapply(obj_list, function(x) map_id(x, ...), simplify=simplify, USE.NAMES=USE.NAMES)
}

if (is.null(module_name())) {
    library(testthat)

    ids = c("TCGA-4Z-AA83-01A", "TCGA-4Z-AA83-01C")
    expect_equal(map_id(ids, subset="primary"), ids[1])
    expect_equal(map_id(ids, subset="primary", drop=FALSE), c(ids[1], NA))
    expect_equal(map_id(ids, id_type="patient"), substr(ids, 1,12))

    idf = data.frame(idx=1:2, ids=ids)
    expect_equal(map_id(idf, along="ids", subset="primary"), idf[1,,drop=FALSE])
    expect_equal(map_id(idf, along="ids", subset="primary", drop=FALSE)$ids, c(ids[1], NA))
    expect_equal(map_id(idf, along="ids", id_type="patient")$ids, substr(ids, 1,12))
}
