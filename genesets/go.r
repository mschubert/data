import_package('dplyr', attach=TRUE)
io = import('ebits/io')

#' Retrieves GO categories for HGNC symbols
#'
#' @param names  Either 'id', 'name', or 'both'
#' @param dset   Dataset; e.g. '{hsapiens,mmusculus}_gene_ensembl'
#' @param genes  Identifier type for genes ('hgnc_symbol', 'entrezgene',
#'               'ensembl_gene_id', etc.)
#' @return  
go = function(dset="hsapiens_gene_ensembl", genes="hgnc_symbol", as_list=FALSE) {
    fname = file.path(module_file("cache", mustWork=TRUE),
                      paste0(paste("go", dset, genes, sep="-"), ".RData"))

    if (file.exists(fname)) {
        sets = io$load(fname)
    } else {
        warning("Creating cache from Biomart, this may take a while",
                immediate.=TRUE)

        mart = biomaRt::useMart(biomart="ensembl", dataset=dset)
        mapGO = biomaRt::getBM(attributes=c(genes, "go_id", "name_1006"), mart=mart)

        sets = mapGO %>%
            filter(go_id != "") %>%
            rename(id = go_id,
                   name = name_1006) 
        sets = sets[sets[[genes]] != "",]

        save(sets, file=fname)
    }

    if (as_list) {
        sets$both = paste(sets$id, sets$name)
        unstack(sets[c(genes, "both")])
    } else
        sets
}

if (is.null(module_name())) {
    cache = go()
}
