config = import('./config')

list_files = function(regex) {
    efiles = list.files(config$raw_data, regex, recursive=TRUE, full.names=TRUE)
    efiles = efiles[file.size(efiles) > 0]
    efiles
}

read_files = function(fnames) {
    read_file = function(fname)
        data.table::fread(paste('zcat', fname), header=TRUE, sep="\t",
                          integer64='double')

    if (length(fnames) == 1)
        read_file(fnames)
    else
        lapply(fnames, read_file)
}

get_matrix = function(fnames, formula, map.hgnc=FALSE, fun.aggregate=sum) {
    ar = import('ebits/array')
    idmap = import('ebits/process/idmap')

    file2matrix = function(fname) {
        message(fname)
        mat = ar$construct(read_files(fname),
                           formula = formula,
                           fun.aggregate=fun.aggregate)
        mat = mat[!rownames(mat) %in% c('?',''),] # aggr fun might handle his?
        if (map.hgnc)
            mat = idmap$gene(mat, to='hgnc_symbol', summarize=fun.aggregate)
        mat
    }
    objs = lapply(fnames, file2matrix)

    if (length(objs) == 0)
        stop("all file reads failed, something is wrong")

    objs
}
