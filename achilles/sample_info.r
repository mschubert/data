io = import('ebits/io')

sample_info = function() {
    fpath = module_file("./data/Achilles_v2.4_SampleInfo_small.txt", mustWork=TRUE)
    io$read_table(fpath, header=TRUE, sep="\t", quote='"')
}
