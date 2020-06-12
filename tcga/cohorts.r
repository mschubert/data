.io = import('ebits/io')
.bc = import('./barcode')

.file = function(...) .io$file_path(module_file(), ...)

#' List all available tissues
cohorts = function() {
    clin = .io$load(module_file("cache", "clinical.RData"))

    unique(.bc$barcode2study(toupper(clin$patient.bcr_patient_barcode)))
}
