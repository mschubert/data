ar = import('ebits/array')
ma = import('ebits/process/microarray')

process_data = function() {
    #TODO:
    # this does not actually work with error:
    # Error in oligo::read.celfiles(filenames = file.path(path, unique(files))) :
    #  The annotation package, pd.u133aaofav2, could not be loaded.
    # (package does not exist on bioc, cf.
    # https://support.bioconductor.org/p/55779/)
    esets = ArrayExpress::ArrayExpress("E-GEOD-5258") %>%
#       ma$qc() %>% # they all pass
        ma$normalize() %>%
        ma$annotate(summarize="hgnc_symbol") %>%
        lapply(exprs)

    expr = ar$stack(esets, along=2)
    save(expr, file=module_file("cache", "expr.RData"))
}
