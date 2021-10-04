.idmap = import('process/idmap')
.rec3 = import('../recount3/tcga')
.map_id = import('../tcga/map_id')$map_id

#' Get cell fraction estimates using different methods
#'
#' @param cohort   The tissue(s) to get expression for
#' @param method   Method to use for estimation (e.g. 'epic', 'cibersort', 'xcell')
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @return         A matrix of TCGA barcodes x estimated cell fractions
cell_frac = function(cohort, method, id_type="specimen") {
    fname = sprintf("cf_%s_%s.rds", cohort, method)
    fpath = file.path(module_file("cache"), fname)

    if (!file.exists(fpath)) {
        tpm = .rec3$tcga(cohort, trans="tpm", id_type="full")
        gnames = .idmap$gene(rownames(tpm), to="hgnc_symbol")
        tpm = narray::map(tpm, along=1, sum, subsets=gnames) #fixme @narray: too slow
        tpm = tpm[,!duplicated(colnames(tpm))]

        library(immunedeconv) # bug?
        re = immunedeconv::deconvolute(tpm, method)
        mat = t(data.matrix(re[-1]))
        colnames(mat) = re$cell_type

        saveRDS(mat, file=fpath)
    } else
        mat = readRDS(fpath)

    rownames(mat) = .map_id(rownames(mat), id_type=id_type)
    mat = mat[,!duplicated(colnames(mat))]
    mat
}
