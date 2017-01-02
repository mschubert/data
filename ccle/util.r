io = import('ebits/io')

#' Function to check whether file exists
exists = function(...) {
    !is.na(module_file(...))
}

#' Function to get a file path relative to the module
file = function(...) {
    io$file_path(module_file(), ...)
}

#' Function to read a table from the module
read = function(..., header=TRUE) {
    io$read_table(module_file(..., mustWork=TRUE), header=header)
}
