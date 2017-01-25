util = import('./process_util')
config = import('./config')
io = import('ebits/io')
ar = import('ebits/array')
rnaseq = import('ebits/process/rna-seq')

#' Process rna_seq data
#'
#' @param force  Overwrite existing files instead of skipping
rna_seq = function() {
    files = util$list_files("^exp_seq")
    files = files[!grepl("LIRI-JP", files)] # 32 bit integer overflow
    exprs = util$get_matrix(files,
                            raw_read_count ~ gene_id + icgc_specimen_id,
                            map.hgnc=TRUE)

    io$save(ar$stack(exprs, along=2),
            file=file.path(config$cached_data, "rna_seq_raw.gctx"))

    transformed = rnaseq$voom(exprs)
    io$save(ar$stack(transformed, along=2),
            file=file.path(config$cached_data, "rna_seq_voom.gctx"))

    transformed = rnaseq$vst(exprs)
    io$save(ar$stack(transformed, along=2),
            file=file.path(config$cached_data, "rna_seq_vst.gctx"))
}

clinical = function() {
    files = util$list_files("^donor\\.")
    clinical = util$read_files(files) %>%
        data.table::rbindlist() %>%
        as_data_frame()

    io$save(clinical, file=file.path(config$cached_data, "clinical.RData"))
}

specimen = function() {
    fname = file.path(config$raw_data, "Summary/specimen.all_projects.tsv.gz")
    specimen = util$read_files(fname)
    io$save(specimen, file=file.path(config$cached_data, "specimen.RData"))
}

mutations = function() {
    mut_aggr = function(x) {
        any(x != 0)
    }
    util$mat('mutations.h5', '^simple_somatic',
             consequence_type ~ gene_affected + icgc_sample_id,
             fun.aggregate = mut_aggr, force=force, map.hgnc=TRUE)
}

cnv = function() {
    util$mat("cnv.h5", '^copy_number_somatic_mutation',
             segment_median ~ gene_affected + icgc_sample_id,
             fun.aggregate = mean, map.hgnc=T, force=force)
}

rppa = function() {
    util$mat("protein.h5", '^protein_expression',
        normalized_expression_level ~ antibody_id + icgc_sample_id,
        map.hgnc=FALSE, force=force)
}
