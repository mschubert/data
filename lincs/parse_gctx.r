#' Returns a subset of a .gctx object as named R matrix
#'
#' @param fname  File name
#' @param rid    Vector of probe IDs
#' @param cid    Vector of experiment IDs
parse_gctx = function(fname, rid=NULL, cid=NULL) {
    rows = gsub("\\ ", "", rhdf5::h5read(fname, "/0/META/ROW/id"))
    cols = gsub("\\ ", "", rhdf5::h5read(fname, "/0/META/COL/id"))

    if (is.null(cid))
        stop("need to specify column (experiment) ids")
    else
        col_i = match(cid, cols)

    if (is.null(rid))
        row_i = 1:length(rows)
    else
        row_i = match(rid, rows)

    structure(rhdf5::h5read(fname, "/0/DATA/0/matrix", index=list(row_i, col_i)),
              .Dimnames = list(rows[row_i], cols[col_i]))
}
