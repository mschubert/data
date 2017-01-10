library(dplyr)
b = import('base')
io = import('io')

process_index = function() {
    index = io$read_table("data/E-GEOD-5258.sdrf.txt", header=TRUE)
    index = index[,!duplicated(colnames(index))]
    index = index %>%
        transmute(id = b$grep("(GSM[0-9]+)", `Source Name`),
                  organism = `Characteristics [organism]`,
                  cell_line = `Characteristics [cell line]`,
                  solvent = `Factor Value [vehicle]`,
                  compound = `Factor Value [compound]`,
                  hours = `Factor Value [time]`,
                  dose = `Factor Value [dose]`)

    io$write_table(index, file.path(module_file("cache"), "index.txt"))
}
