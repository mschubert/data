util = import('./process_util')
config = import('./config')

#' Process rna_seq data
#'
#' @param force  Overwrite existing files instead of skipping
rna_seq = function(force=FALSE) {
    util$mat("expr_seq_raw.h5", '^exp_seq',
        raw_read_count ~ gene_id + icgc_sample_id, map.hgnc=FALSE, force=force)

#    util$mat("expr_seq_norm.h5", '^exp_seq',
#        normalized_read_count ~ gene_id + icgc_sample_id, map.hgnc=TRUE, force=force)
    voomfile = file.path(config$cached_data, "expr_seq_voom.h5")
    if (identical(force, TRUE) || !file.exists(voomfile)) {
        expr = getRNASeq(raw.counts=TRUE) %>% na.omit()
        expr = limma::voom(expr)$E
        h5store::h5save(t(expr), file=voomfile)
    }
}

clinical = function(force=FALSE) {
    tfun = function(x) mutate(x, tissue = .b$grep("^(\\w+)", project_code))
    util$df("clinical.RData", "clinical\\.", transform=tfun, force=force)
    util$df("clinicalsample.RData", "clinicalsample\\.", force=force)
}

mutations = function(force=FALSE) {
    mut_aggr = function(x) {
        any(x != 0)
    }
    util$mat('mutations.h5', '^simple_somatic',
             consequence_type ~ gene_affected + icgc_sample_id,
             fun.aggregate = mut_aggr, force=force, map.hgnc=TRUE)
}

cnv = function(force=FALSE) {
    util$mat("cnv.h5", '^copy_number_somatic_mutation',
             segment_median ~ gene_affected + icgc_sample_id,
             fun.aggregate = mean, map.hgnc=T, force=force)
}

rppa = function(force=FALSE) {
    util$mat("protein.h5", '^protein_expression',
        normalized_expression_level ~ antibody_id + icgc_sample_id,
        map.hgnc=FALSE, force=force)
}
