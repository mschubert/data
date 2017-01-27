#' Subset IDs or matrices by study, type, etc.
#'
#' @param study    Character vector of studies to include (default: all)
#' @param primary  IDs must be primary (either cancer or normal)
#' @param cancer   IDs must be cancer
#' @param normal   IDs must be normals
#' @return         A filtered vector or matrix
filter = function(x, study=NULL, primary=NULL, cancer=NULL, normal=NULL) {
    sp = import('./specimen')$specimen()

    if (!is.null(study))
        sp = dplyr::filter(sp, project_code %in% study)

    if (!is.null(primary)) {
        matched_normal = c(
            "Primary tumour - blood derived (bone marrow)" = "Normal - bone marrow",
            "Primary tumour - blood derived (peripheral blood)" = "Normal - blood derived",
            "Primary tumour - lymph node" = "Normal - lymph node",
            "Primary tumour - solid tissue" = "Normal - tissue adjacent to primary",
            "Primary tumour - solid tissue" = "Normal - solid tissue"
        )
        valid = c(matched_normal, names(matched_normal))

        pri_and_matched = function(x) {
            pri = unique(grep("^Primary", x, value=TRUE))
            x %in% pri | x %in% matched_normal[pri]
        }

        sp = dplyr::filter(sp, specimen_type %in% valid) %>%
            dplyr::group_by(project_code) %>%
            dplyr::filter(pri_and_matched(specimen_type)) %>%
            dplyr::ungroup()
    }

    if (!is.null(cancer)) {
        valid = "(^Primary)|(^Recurrent)|(^Metastatic)"
        sp = dplyr::filter(sp, grepl(valid, specimen_type))
    }

    if (!is.null(normal))
        sp = dplyr::filter(sp, grepl("^Normal -", specimen_type))

    intersect(x, sp$icgc_specimen_id)
}

if (is.null(module_name())) {
    library(testthat)

    mysp = c(
        "SP41621", # KIRC-US, adjacent normal
        "SP41770", # KIRC-US, solid primary
        "SP131416", # ALL-US, blood primary
        "SP822", # ALL-US, normal bone marrow
        "SP120560" # SKCM-US, distant metastasis
    )

    expect_equal(filter(mysp, study="KIRC-US"), mysp[c(1,2)])
    expect_equal(filter(mysp, study="ALL-US"), mysp[c(3,4)])
    expect_equal(filter(mysp, study=c("SKCM-US", "ALL-US")), mysp[c(3,4,5)])

    expect_equal(filter(mysp, primary=TRUE), mysp[c(1:4)])
    expect_equal(filter(mysp, cancer=TRUE), mysp[c(2,3,5)])
    expect_equal(filter(mysp, normal=TRUE), mysp[c(1,4)])

    expect_equal(filter(mysp, study="KIRC-US", primary=TRUE), mysp[c(1,2)])
    expect_equal(filter(mysp, study="KIRC-US", cancer=TRUE), mysp[2])
    expect_equal(filter(mysp, study="KIRC-US", normal=TRUE), mysp[1])

    expect_equal(filter(mysp, study="SKCM-US", primary=TRUE), character())
    expect_equal(filter(mysp, study="SKCM-US", cancer=TRUE), mysp[5])
}
