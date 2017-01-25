.io = import('ebits/io')

.file = function(...) .io$file_path(module_file(), ...)

#' List all available tissues
tissues = function(id_type="specimen") {
    file = h5::h5file(module_file("cache", "rna_seq2_vst.gctx"), mode="r")

    barcodes = file["/0/META/ROW/id"][]
    studies = .bc$barcode2study(barcodes)
    re = unique(sort(studies))
    h5::h5close(file)
    re
}
