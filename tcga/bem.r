b = import('ebits/base')
io = import('ebits/io')

bem = function() {
#    bem = paste0("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/", c(
#        mut = "BEMs/PrimaryTumours/PrimTum_CG_BEMs.zip", # link broken
#        cna = "suppData/TableS2F.xlsx",
#        meth = "BEMs/PrimaryTumours/PrimTum_METH_BEMs.zip"
#    ))

    # only present when all data avail -> e.g. only 400 BRCA samples
    io$load(module_file("cache", "BEM.RData"))
}

process = function(from, to) {
    INFILE = commandArgs(TRUE)[1] %or% "data/PrimTum_MoBEM_PanCan.tsv"
    OUTFILE = commandArgs(TRUE)[1] %or% "cache/bem.RData"

    bem = t(io$read_table(INFILE, header=TRUE))

    save(bem, file=OUTFILE)
}
