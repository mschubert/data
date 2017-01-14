config = import('./config')
ar = import('ebits/array')

mat = function(regex, formula, map.hgnc=FALSE, force=FALSE, fun.aggregate=sum) {
    if (!force && file.exists(file.path(config$cached_data, fname)))
        return()

    efiles = list.files(config$raw_data, regex, recursive=TRUE, full.names=TRUE)
    efiles = efiles[file.size(efiles) > 0]
    idmap = import('ebits/process/idmap')

    file2matrix = function(fname) {
        message(fname)
        mat = ar$construct(data.table::fread(paste('zcat', fname), header=TRUE, sep="\t"),
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

    mat = ar$stack(obj, along=2) #TODO: check if saving the right t()
}

df = function(fname, regex, transform=function(x) x, force=FALSE) {
    if (!force && file.exists(file.path(config$cached_data, fname)))
        return()

    files = list.files(config$raw_data, regex, recursive=TRUE, full.names=TRUE)
    mat = do.call(rbind, lapply(files, function(f) cbind(
        study = b$grep("/([^/]+)/[^/]+$", f),
        data.table::fread(paste('zcat', f), header=TRUE, sep="\t")
    ))) %>% transform()
}
