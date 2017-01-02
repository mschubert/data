gdsc = import('../gdsc')
util = import('./util')

index = util$read("data", "CCLE_sample_info_file_2012-10-18.txt", header=TRUE)
index$COSMIC = gdsc$cosmic$name2id(index$`Cell line primary name`, warn=FALSE)
