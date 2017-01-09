ar = import('array')
pryr = import_package('pryr')
bc = import('./barcode')

#' Function to intersect barcodes
#'
#' This is an array/intersect with additional shortest barcode
#'
#' @param ...    Objects to intersect
#' @param along  Which dimension to intersect
#' @param envir  The environment in which to perform
#' @return       Objects with the same identifiers, rest dropped
intersect = function(..., along=1, envir=parent.frame()) {
    dots = pryr$named_dots(...)
	objs = setNames(lapply(dots, function(x) eval(x, envir=envir)), names(dots))
    all_bc = unlist(lapply(objs, function(x) {
        if (is.vector(x) && is.null(names(x)))
            x
        else
            ar$dimnames(x, along=along)
    }))
    bc$must_barcode(all_bc)
    min_len = min(nchar(all_bc))

    # use shortest barcodes
    for (i in seq_along(dots)) {
        if (is.call(dots[[i]])) {
            if (along == 1 && is.data.frame(eval(dots[[i]][[2]], envir=envir))) {
                obj = as.data.frame(eval(dots[[i]][[2]], envir=envir))
                field = as.character(dots[[i]][[3]])
                obj[[field]] = substr(obj[[field]], 1, min_len)
                assign(as.character(dots[[i]][[2]]), obj, envir=envir)
            } else
                stop("calls can only reference `data.frame` fields with along=1")
        } else {
            obj = as.array(objs[[i]])
            dimnames(obj)[[along]] = substr(dimnames(obj)[[along]], 1, min_len)
            assign(names(dots)[i], obj, envir=envir)
        }
    }

    # intersect the arrays
    ar$intersect(..., along=along, envir=envir)
}

if (is.null(module_name())) {
    library(testthat)

    dn = list(c('TCGA-OR-A5J1-01A','TCGA-OR-A5J2-01A'),
              c('TCGA-OR-A5J2','TCGA-OR-A5J1-01A-11D'))
    a = setNames(1:2, rev(dn[[1]]))
    A = t(matrix(1:4, nrow=2, ncol=2, dimnames=dn))
    DF = structure(list(y=3:4, z=c(6,5), x=1:2, A=dn[[1]]),
            .Names=c("y","z","x","A"), row.names=1:2, class="data.frame")

	a1 = a
    DF2 = DF
	intersect(a1, A, along=1)
	expect_equal(names(a1), rownames(A))
	intersect(a1, DF2$A)
	expect_equal(names(a1), DF2$A)

    # check min_len also takes into account DF$A
    DF3 = DF2
    a2 = setNames(1:2, dn[[1]])
	intersect(a2, DF3$A)
	expect_equal(names(a2), DF3$A)
}
