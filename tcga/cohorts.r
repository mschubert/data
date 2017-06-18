.io = import('ebits/io')
.bc = import('./barcode')

.file = function(...) .io$file_path(module_file(), ...)

#' List all available tissues
cohorts = function() {
    file = h5::h5file(module_file("cache", "rna_seq2_vst.gctx"), mode="r")

    #TODO: check right transpose in h5 file
    barcodes = file["/0/META/COL/id"][]
    studies = .bc$barcode2study(barcodes)
    re = unique(sort(studies))
    h5::h5close(file)
    re
}
