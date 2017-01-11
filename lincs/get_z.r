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

    if (is.character(map_genes)) {
        annot = import('./probe_annotations')$probe_annotations()
        ar$summarize(re, along=1, from="affy_hg_u133_plus_2", to=map_genes,
                     data=annot, FUN=function(x) mean(x, na.rm=TRUE))
    } else
        re
}

if (is.null(module_name())) {
    library(testthat)

    id = "ASG001_MCF7_24H_X1_B7_DUO52HI53LO:G13"
    z = get_z(id)
    expect_equal(colnames(z), id)
#    expect_equal(rownames(z), probes$landmarks) # 978 vs 977?

    z2 = get_z(id, map_genes="hgnc_symbol")
    expect_equal(colnames(z2), id)
    expect_lt(nrow(z2), nrow(z))
}
