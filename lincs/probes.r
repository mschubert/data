.load_fun = function(fname)
    as.character(read.table(module_file('cache', fname, mustWork=TRUE))$V1)

#' Landmark probes
landmarks = .load_fun('rid_landmarks.txt')

#' Projected probes
projected = .load_fun('rid_projected.txt')

#' Best inferred probes
bing = .load_fun('rid_bing.txt')
