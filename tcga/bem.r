b = import('ebits/base')
io = import('ebits/io')

bem = function() {
    io$load(module_file("cache", "BEM.RData"))
}

process = function(from, to) {
    INFILE = commandArgs(TRUE)[1] %or% "data/PrimTum_MoBEM_PanCan.tsv"
    OUTFILE = commandArgs(TRUE)[1] %or% "cache/bem.RData"

    bem = t(io$read_table(INFILE, header=TRUE))

    save(bem, file=OUTFILE)
}
