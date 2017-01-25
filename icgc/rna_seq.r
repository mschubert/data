config = import('./config')
subs = import('./subset')

#' Get a matrix for all RNA-seq measurements
#'
#' @param tissue   The tissue(s) to get expression for
#' @return         A matrix with HGNC symbols x TCGA samples
rna_seq = function(study) {
    library(methods) # required; otherwise h5 error
    fname = file.path(config$cached_data, "rna_seq_vst.gctx")
    file = h5::h5file(fname, mode="r")

    ids = file["/0/META/COL/id"][]
    keep = ids %in% subs$subset(ids, study=study)

    data = file["/0/DATA/0/matrix"][which(keep),]
    rownames(data) = ids[keep]
    colnames(data) = file["/0/META/ROW/id"][]

    # remove genes with all NA
    data = data[,colSums(!is.na(data)) != 0]

    h5::h5close(file)
    t(data)
}
