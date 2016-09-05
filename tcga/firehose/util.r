# read raw data from .txt.gz files
# save into R objects for quicker loading
data_dir = module_file("data", "stddata__2016_01_28")
analyses_dir = module_file("data", "analyses__2016_01_28")

list_files = function(file_regex, dir=data_dir) {
    list.files(dir, pattern=file_regex, full.names=TRUE, recursive=TRUE)
}

unpack = function(files, unpack_dir=tempdir()) {
    #TODO: only unpack if file is not there
    for (file in files)
        untar(file, exdir = unpack_dir)
    file.path(unpack_dir, sub(".tar", "", sub(".gz", "", basename(files))))
}

select = function(newdir, tar_regex) {
    list.files(newdir, pattern=tar_regex, full.names=TRUE, recursive=TRUE)
}

voom_transform = function(mat, ids="hgnc") {
    if (ids != "hgnc")
        stop("not implemented")

    # apply voom transformation
    mat = limma::voom(mat)$E
    rownames(mat) = sub("\\|[0-9]+$", "", rownames(mat))
    mat = mat[rownames(mat) != "?",]
    limma::avereps(mat)
}
