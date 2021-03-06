.b = import('ebits/base')
.bc = import('./barcode')

#' Function to filter TCGA identifiers or arrays
#'
#' Will work on character vectors ()
#'
#' @param x        An array or data.frame
#' @param along    The dimension index, name, or field name to filter along
#' @param vial     Limit to a portion (e.g. 'A', 'B', etc., 1 if lowest letter or NULL for all)
#' @param primary  Limit to primary samples (solid tumor, blood cancer, bone marrow)
#' @param normal   Limit to normals of the same tissue
#' @param blood_normal  Limi to normals derived from blood
#' @param tumor    Limit to tumors
#' @param tissue   Character vector of TCGA identifiers
#' @param matched  Return only samples for which there is a matched cancer/normal available
#' @return         The object subset with barcodes that match filters
filter = function(x, along=NULL, vial=NULL, primary=NULL, normal=NULL,
                  blood_normal=NULL, cancer=NULL, tissue=NULL, matched=NULL,
                  include_cell_lines=FALSE, include_xenografts=FALSE) {

    UseMethod("filter")
}

filter.data.frame = function(x, along, vial=NULL, primary=NULL, normal=NULL,
                             blood_normal=NULL, cancer=NULL, tissue=NULL, matched=NULL,
                             include_cell_lines=FALSE, include_xenografts=FALSE) {

    args = as.list(.b$match_call_defaults())[-1]
    args$x = x[[args[["along"]]]]
    args[[2]] = NULL
    x[do.call(filter, args),]
}

# matrices and arrays
filter.matrix = function(x, along=2, vial=NULL, primary=NULL, normal=NULL,
                         blood_normal=NULL, cancer=NULL, tissue=NULL, matched=NULL,
                         include_cell_lines=FALSE, include_xenografts=FALSE) {

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

# vectors and lists
filter.default = function(x, along=NULL, vial=NULL, primary=NULL, normal=NULL,
                          blood_normal=NULL, cancer=NULL, tissue=NULL, matched=NULL,
                          include_cell_lines=FALSE, include_xenografts=FALSE) {

    .bc$must_barcode(x)
    fields = .bc$barcode2index(x)
    keep = rep(TRUE, length(x))

#    eq = function(x, cond) ifelse(x, cond, !cond)
    eq = function(x, cond) {
        if (x)
            cond
        else
            !cond
    }

    if (!is.null(vial)) {
        if (is.character(vial))
            keep = keep & fields$Vial == vial
        else if (is.numeric(vial))
            stop("not implemented")
        else
            stop("Invalid 'vial': ", vial)
    }

    if (!is.null(primary))
        keep = keep & eq(primary, grepl("Primary|Normal", fields$Sample.Definition))

    if (!is.null(normal))
        keep = keep & eq(normal, grepl("Normal", fields$Sample.Definition))

    #TODO: should be included in 'normal' if blood cancer?
    if (!is.null(blood_normal))
        keep = keep & eq(blood_normal, grepl("Blood", fields$Sample.Definition))

    if (!is.null(cancer))
        keep = keep &  eq(cancer, grepl("Tumor|Cancer|Metastatic", fields$Sample.Definition))

    if (!is.null(tissue))
        keep = keep & fields$Study.Abbreviation %in% tissue

    if (!is.null(matched)) {
		matched = fields %>%
            dplyr::group_by(TSS.Code, Participant.ID) %>%
            dplyr::mutate(keep = "TP" %in% Short.Letter.Code &
                                 "NT" %in% Short.Letter.Code)
        keep = keep & matched$keep
	}

    if (sum(keep) == 0)
        warning("No entries left after filtering TCGA barcodes\n", immediate.=TRUE)

    setNames(keep, x)
}

if (is.null(module_file())) {
    library(testthat)
    tumors = c("TCGA-DA-A3F3-06A", "TCGA-BR-4257-01B", "TCGA-E2-A1LE-01C",
        "TCGA-AA-3984-01A", "TCGA-D9-A6EC-06A", "TCGA-BR-8591-01A", "TCGA-FW-A3R5-06A",
        "TCGA-D1-A103-01A", "TCGA-91-6828-01A", "TCGA-L5-A4OI-01A")

    # without specifying filter, return all TRUE and barcodes as names
    expect_equal(all <- filter(tumors), setNames(rep(TRUE, length(tumors)), tumors))
    expect_equal(tumors, names(all))
    expect_equal(filter(tumors, cancer=TRUE), all)

    # tissue filtering should return same as barcode info
    expect_equal(filter(tumors, tissue="SKCM"),
                 .bc$barcode2study(tumors) == "SKCM")

    # primary tumors only
    expect_equal(unname(filter(tumors, primary=TRUE)),
                 .bc$barcode2index(tumors)$Sample.Definition != "Metastatic")

    # first vial available
#    expect_equal(filter(tumors, vial=1), all)

    # if we want normals this should be all false
    expect_warning(empty <- filter(tumors, normal=TRUE))
    expect_equal(empty, setNames(rep(FALSE, length(tumors)), tumors))

    # matrix, data.frame functions
    tmat = matrix(1:30, nrow=10, ncol=3, dimnames=list(tumors, NULL))
    tdf = data.frame(barcode = tumors, a="a", b=5)
    expect_equal(rownames(filter(tmat, tissue="SKCM")),
                 filter(tdf, "barcode", tissue="SKCM")$barcode,
                 names(filter(tumors, tissue="SKCM")))

    normals = c("TCGA-AA-3984-11A", "TCGA-55-7726-10A-01D", "TCGA-JY-A6FA-10A-01D",
        "TCGA-26-5134-10A", "TCGA-D1-A0ZS-10A", "TCGA-FW-A3R5-11A-11D",
        "TCGA-BS-A0UV-10A", "TCGA-78-7155-10A-01D")

    expect_equal(filter(normals),
                 filter(normals, normal=TRUE),
                 filter(normals, normal=TRUE, vial="A"))

    matched = c(tumors, normals)
    is_true = function(x) names(x)[x]
    pos = c("TCGA-AA-3984-01A", "TCGA-AA-3984-11A")
    expect_equal(is_true(filter(matched, matched=TRUE)), pos)
    expect_equal(is_true(filter(matched, matched=TRUE, cancer=TRUE)), pos[1])
    expect_equal(is_true(filter(matched, matched=TRUE, normal=TRUE)), pos[2])
}
