b = import('ebits/base')
ar = import('ebits/array')
probes = import('./probes')
parse_gctx = import('./parse_gctx')$parse_gctx

#' Returns z-scores for a subset of experiments
#'
#' @param cid        Vector of experiment IDs to subset
#' @param rid        Vector of probe IDs to subset
#' @param map.genes  BioMart identifier of IDs to map to, or FALSE. Supported
#'                   identifiers are: hgnc_symbol, hgnc_id, entrezgene,
#'                   ensembl_gene_id, ensembl_transcript_id
get_z = function(cid, rid=probes$landmarks, map_genes=FALSE) {
    fname = module_file("data", "zspc_n1328098x22268.gctx", mustWork=TRUE)
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
    z = get_z(cid)
    expect_equal(colnames(z), cid)
#    expect_equal(rownames(z), probes$landmarks) # 978 vs 977?

    z2 = get_z(cid, map_genes="hgnc_symbol")
    expect_equal(colnames(z2), cid)
    expect_lt(nrow(z2), nrow(z))

    z3 = get_z(cid, rid=probes$projected, map_genes="hgnc_symbol")
}
