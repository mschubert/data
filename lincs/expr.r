b = import('ebits/base')
ar = import('ebits/array')
probes = import('./probes')
parse_gctx = import('./parse_gctx')$parse_gctx

#' Returns expression values for a subset of experiments
#'
#' Inference using deep neural networks, as in this articicle:
#' https://doi.org/10.1093/bioinformatics/btw074
#' Data downloaded from:
#' https://cbcl.ics.uci.edu/public_data/D-GEX/l1000_n1328098x22268.gctx 
#'
#' @param cid        Vector of experiment IDs to subset
#' @param rid        Vector of probe IDs to subset
#' @param map.genes  BioMart identifier of IDs to map to, or FALSE. Supported
#'                   identifiers are: hgnc_symbol, hgnc_id, entrezgene,
#'                   ensembl_gene_id, ensembl_transcript_id
expr = function(cid, rid=probes$landmarks, map_genes=FALSE) {
    fname = module_file("data", "q2norm_n1328098x22268.gctx", mustWork=TRUE)
    re = parse_gctx(fname=fname, cid=cid, rid=rid)

    if (!identical(map_genes, FALSE)) {
        annot = import('./probe_annotations')$probe_annotations()

        if (! map_genes %in% colnames(annot))
            stop(map_genes, " not found in annotation mapping (available: ",
                 paste(colnames(annot), collapse=", "), ")")

        mapped = b$match(rownames(re),
                         from = annot$affy_hg_u133_plus_2,
                         to = annot[[map_genes]])
        rownames(re) = unname(mapped)
        re = re[!is.na(rownames(re)),,drop=FALSE]
        limma::avereps(re)
    } else
        re
}

if (is.null(module_name())) {
    library(testthat)

    cid = "ASG001_MCF7_24H_X1_B7_DUO52HI53LO:G13"
    z = expr(cid)
    expect_equal(colnames(z), cid)
#    expect_equal(rownames(z), probes$landmarks) # 978 vs 977?

    z2 = expr(cid, map_genes="hgnc_symbol")
    expect_equal(colnames(z2), cid)
    expect_lt(nrow(z2), nrow(z))

    z3 = expr(cid, rid=probes$projected, map_genes="hgnc_symbol")
}
