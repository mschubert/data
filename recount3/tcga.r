.map_id = import('../tcga/map_id')$map_id

#' Get TCGA data processed by recount3 project
#'
#' @param cohort   The tissue(s) to get expression for
#' @param id_type  Where to cut the barcode, either "patient", "specimen", or "full"
#' @param trans    Count transformation (only cur avail: "tpm")
#' @param annot    defunct
#' @return         A matrix with Ensembl gene ID x TCGA sample barcode
tcga = function(cohort, id_type="specimen", trans="tpm", annot=FALSE) {
    library(recount3) # https://github.com/LieberInstitute/recount3/issues/7

    re = recount3::create_rse_manual(
        project = cohort,
        project_home = "data_sources/tcga",
        organism = "human",
        annotation = "gencode_v29",
        type = "gene"
    )

    SummarizedExperiment::assays(re)$counts = recount3::transform_counts(re) # library size norm
    mat = recount::getTPM(re)

    rownames(mat) = sub("\\.[0-9]$", "", rownames(mat)) # strip Ensembl versions
    colnames(mat) = SummarizedExperiment::colData(re)$tcga.tcga_barcode
    colnames(mat) = .map_id(colnames(mat), id_type=id_type)

    mat
}
