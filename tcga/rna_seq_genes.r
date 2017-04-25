#' Get all valid gene symbols
#'
#' @return  A character vector of HGNC symbols
rna_seq_genes = function() {
    library(methods) # required; otherwise h5 error
    file = h5::h5file(module_file("cache", "rna_seq2_vst.gctx"), mode="r")
    genes = file["/0/META/ROW/id"][]
    h5::h5close(file)

    genes
}

if (is.null(module_name())) {
    library(testthat)

    hgnc = rna_seq_genes()

    expect_true(is.character(hgnc))
    expect_true(all(c("A1BG", "A1CF", "ZZEF1") %in% hgnc))
}
