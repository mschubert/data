io = import('io')
.filter = import('./filter')$filter

#' Returns the RNA-seq expression matrix
#'
#' @param tissue  Only return from tissue
#' @param sample  Only return sample IDs
#' @param gene    Only return HGNC symbols
#' @return        A gene expression matrix (genes x samples)
rna_seq = function(tissue=NULL, sample=NULL, gene=NULL) {
	fpath = module_file("cache", "rna_seq_vst.gctx", mustWork=TRUE)

    file = h5::h5file(fpath, mode="r")

    ids = file["/0/META/COL/id"][]

    if (!is.null(tissue))
        keep = ids %in% .filter(tissue=tissue)

    if (!is.null(sample))
        keep = ids & ids %in% sample

    if (is.null(gene))
        data = file["/0/DATA/0/matrix"][which(keep),]
    else
        data = file["/0/DATA/0/matrix"][which(keep),gene]

    rownames(data) = ids[keep]
    colnames(data) = file["/0/META/ROW/id"][]

    # remove genes with all NA
    data = data[,colSums(!is.na(data)) != 0]

    h5::h5close(file)
    t(data)
}
