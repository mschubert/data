.io = import('io')

drug_props = function() {
    re = .io$read_table(module_file('drug_props.txt'), header=TRUE)
    re
}

combos = function() {
    re = .io$read_table(module_file('combo_resp.txt'), header=TRUE)
    re
}
