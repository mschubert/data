io = import('ebits/io')

sh_genes = function() {
    fpath = module_file("./data/Achilles_v3.3.8.Gs.gct", mustWork=TRUE)
    re = io$read_gct(fpath)@mat
    rownames(re) = sub("_.*$", "", rownames(re))
    re
}
