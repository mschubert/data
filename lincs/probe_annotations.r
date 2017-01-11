io = import('ebits/io')

probe_annotations = function() {
    fname = module_file("cache", 'annot.RData', mustWork=TRUE) #FIXME: get path also when the file doesn't exist
    if (!file.exists(fname)) {
        mart = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        annot = biomaRt::getBM(attributes = c("hgnc_symbol", "hgnc_id", "entrezgene",
            "ensembl_gene_id", "ensembl_transcript_id", "affy_hg_u133_plus_2"),
            filter = "affy_hg_u133_plus_2", values = projected, mart=mart)
        annot[annot==''] = NA

        save(annot, file=fname)
        annot
    } else
        io$load(fname)
}
