io = import('ebits/io')
probes = import('./probes')
mart = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

affy_annots = c(
    'affy_hg_u133_plus_2', # 20052 of 22268
    "affy_hg_focus",
    "affy_hg_u133_plus_2",
    "affy_hg_u133a_2",
    "affy_hg_u133a",
    "affy_hg_u133b",
    "affy_hg_u95av2",
    "affy_hg_u95b",
    "affy_hg_u95c",
    "affy_hg_u95e",
    "affy_hg_u95d",
    "affy_hg_u95a",
    "affy_hugenefl",
#   "affy_hta_2_0",
#    "affy_huex_1_0_st_v2",
#    "affy_hugene_1_0_st_v1",
    "affy_primeview",
#    "affy_hugene_2_0_st_v1",
    "affy_u133_x3p"
)

get_annot = function(probe_array) {
    annot = biomaRt::getBM(attributes = c("hgnc_symbol", "hgnc_id", "entrezgene",
        "ensembl_gene_id", "ensembl_transcript_id", probe_array),
        filter = probe_array, values = probes$projected, mart=mart)
    annot[annot==''] = NA
    annot
}

probe_annotations = function() {
    fname = module_file("cache", 'annot.RData')
    if (!file.exists(fname)) {
        annot = get_annot('affy_hg_u133_plus_2')
        save(annot, file=fname)
        annot
    } else
        io$load(fname)
}
