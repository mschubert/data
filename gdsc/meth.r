.cosmic = import('./cosmic')

meth = function() {
    fname = file.path(module_file("cache"), "meth.rds")
    if (!file.exists(fname)) {
        raw = read.csv(file.path(module_file("data"), "GSE68379_Matrix.processed.txt.gz"), sep="\t")
        meth = data.matrix(raw[-1])
        rownames(meth) = raw$Row.names
        colnames(meth) = unname(.cosmic$name2id(sub("_AVG.Beta", "", sub("^X", "", colnames(raw)[-1]))))

        saveRDS(meth, file=fname)
    } else {
        meth = readRDS(fname)
    }

    meth
}
