.file = import('./file')
tissues = import('./tissues')$tissues

#' Returns a matrix of frequency- and intogen-filtered mutated genes
mutated_genes = function(frequency=0, intogen=FALSE, tissue=NULL, drop=FALSE) {
    mut = t(.file$get('NGS_BEM')$logical)

    if (!is.null(tissue))
        mut = mut[rownames(mut) %in% names(tissues(tissue)),]

    if (intogen) {
        drivers = drivers(tissue=tissue)
        mut = mut[,intersect(unique(drivers$HGNC), colnames(mut))]
    }

    if (frequency > 0)
        mut = mut[,colSums(mut)/nrow(mut) > frequency]

    if (drop)
        mut = mut[,colSums(mut)>0]

    mut
}
