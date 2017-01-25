#' Returns a list containing row- and column names for clinical data
#TODO: give index=, infer id type from it
clinical = function(specimen=NULL, donors=NULL) {
    re = idx$getRData('clinical.RData')

    if (!is.null(donors))
        re[match(donors, re$icgc_donor_id),]
    else if (!is.null(specimen))
        re[match(specimen, re$icgc_specimen_id),]
    else
        re
}

#' Returns a list containing row- and column names for clinical sample data
clinical_sample = function() idx$getRData('clinicalsample.RData')
