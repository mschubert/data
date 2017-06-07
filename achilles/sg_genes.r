io = import('ebits/io')

sg_genes = function() {
    fpath = module_file("./data/Achilles_QC_v2.4.3.rnai.Gs.gct", mustWork=TRUE)
    re = io$read_gct(fpath)@mat
    rownames(re) = sub("_.*$", "", rownames(re))
    re
}
