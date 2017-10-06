import_package('dplyr', attach=TRUE)
io = import('ebits/io')

#' Retrieves GO categories for HGNC symbols
#'
#' @param names  Either 'id', 'name', or 'both'
#' @param genes  Identifier type for genes ('hgnc_symbol', 'entrezgene',
#'               'ensembl_gene_id', etc.)
#' @return  
go = function(names="both", genes="hgnc_symbol") {
    fname = file.path(module_file("cache", mustWork=TRUE),
                      paste0("go-", genes, ".RData"))

    if (file.exists(fname)) {
        sets = io$load(fname)
    } else {
        warning("Creating cache from Biomart, this may take a while",
                immediate.=TRUE)

        mart = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        mapGO = biomaRt::getBM(attributes=c(genes, "go_id", "name_1006"), mart=mart)

        sets = mapGO %>%
            filter(go_id != "") %>%
            rename(id = go_id,
                   name = name_1006) 
        sets = sets[sets[[genes]] != "",]

        save(sets, file=fname)
    }

    if (names == "both")
        sets$both = paste(sets$id, sets$name)

    unstack(sets[c(genes, names)])
}

if (is.null(module_name())) {
    cache = go()
}
