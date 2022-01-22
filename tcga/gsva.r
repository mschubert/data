.gset = import('genesets')
.idmap = import('process/idmap')
.rnaseq = import('./rna_seq')$rna_seq
.map_id = import('../tcga/map_id')$map_id

gsva = function(cohort, setname, id_type="specimen") {
    fdir = file.path(module_file(), "cache", "gsva", cohort)
    fpath = file.path(fdir, sprintf("%s.rds", setname))
    is_std_coll = is.character(setname) && length(setname == 1)

    if (!file.exists(fpath)) {
        if (is_std_coll) {
            dir.create(fdir, showWarnings=FALSE, recursive=TRUE)
            sets = .gset$get_human(setname)
        } else
            sets = setname

        if (cohort == "pan") {
            expr = .rnaseq(cohort, id_type="full", trans="raw")
            eset = DESeq2::DESeqDataSetFromMatrix(expr, colData=data.frame(id=colnames(expr)), design=~1)
            expr = SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(eset))
        } else
            expr = .rnaseq(cohort, id_type="full", trans="vst")

        rownames(expr) = .idmap$gene(rownames(expr), to="hgnc_symbol")
        scores = GSVA::gsva(expr, sets, parallel.sz=10)
        if (is_std_coll)
            saveRDS(scores, file=fpath)
    } else
        scores = readRDS(fpath)

    colnames(scores) = .map_id(colnames(scores), id_type=id_type)
    scores = scores[,!duplicated(colnames(scores))]
    scores
}

if (is.null(module_name())) {
    library(testthat)
    scores = gsva("ACC", "MSigDB_Hallmark_2020")
    expect_is(scores, "array")
}
