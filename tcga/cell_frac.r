idmap = import('process/idmap')
.rec3 = import('../recount3/tcga')

#' Get cell fraction estimates using different methods
#'
#' @param cohort  The tissue(s) to get expression for
#' @param method  Method to use for estimation (e.g. 'epic', 'cibersort', 'xcell')
#' @return        A matrix of TCGA barcodes x estimated cell fractions
cell_frac = function(cohort, method) {
    fname = sprintf("cf_%s_%s.rds", cohort, method)
    fpath = file.path(module_file("cache"), fname)

    if (!file.exists(fpath)) {
        tpm = .rec3$tcga(cohort, trans="tpm")
        gnames = idmap$gene(rownames(tpm), to="hgnc_symbol")
        tpm = narray::map(tpm, along=1, sum, subsets=gnames)

        re = immunedeconv::deconvolute(tpm, method)
        mat = t(data.matrix(re[-1]))
        colnames(mat) = re$cell_type

        saveRDS(mat, file=fpath)
    } else
        mat = readRDS(fpath)

    mat
}
