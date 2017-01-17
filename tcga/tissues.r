.io = import('ebits/io')

.file = function(...) .io$file_path(module_file(), ...)

#' List all available tissues
tissues = function(id_type="specimen") {
    tt = list.files(.file("cache"), pattern="_voom\\.RData")
    sub("_voom.RData", "", tt)
}
