.io = import('ebits/io')

#' Returns a binary event matrix (BEM) for mutated genes
bem = function() {
    t(.io$load(module_file("cache", "BEMs", 'PANCAN_simple_MOBEM.rdata', mustWork=TRUE)))
}
