config = import('./config')

mat = function(fname, regex, formula, map.hgnc=FALSE, force=FALSE, fun.aggregate=sum) {
    if (!force && file.exists(file.path(config$cached_data, fname)))
        return()

    idmap = import('ebits/process/idmap')
    efiles = list.files(config$raw_data, regex, recursive=T, full.names=T)

    file2matrix = function(fname) {
        message(fname)
        mat = .ar$construct(read.table(fname, header=T, sep="\t"),
                           formula = formula,
                           fun.aggregate=fun.aggregate)
        mat = mat[!rownames(mat) %in% c('?',''),] # aggr fun might handle his?
        if (map.hgnc)
            mat = idmap$gene(mat, to='hgnc_symbol', summarize=fun.aggregate)
        mat
    }
    obj = lapply(efiles, file2matrix)

    if (length(obj) == 0)
        stop("all file reads failed, something is wrong")

    mat = t(.ar$stack(obj, along=2, fill=0))
    if (grepl("\\.h5$", fname))
        h5store::h5save(mat, file=file.path(config$cached_data, fname))
    else
        save(mat, file=file.path(config$cached_data, fname))
}

df = function(fname, regex, transform=function(x) x, force=FALSE) {
    if (!force && file.exists(file.path(config$cached_data, fname)))
        return()

    files = list.files(config$raw_data, regex, recursive=TRUE, full.names=TRUE)
    mat = do.call(rbind, lapply(files, function(f) cbind(
        study = b$grep("/([^/]+)/[^/]+$", f),
        read.table(f, header=TRUE, sep="\t") #FIXME: io$read should work here
    ))) %>% transform()

    if (grepl("\\.h5$", fname))
        h5store::h5save(mat, file=file.path(config$cached_data, fname))
    else
        save(mat, file=file.path(config$cached_data, fname))
}
