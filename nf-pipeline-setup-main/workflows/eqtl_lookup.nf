/*
    EQTL_LOOKUP_WORKFLOW

    Deze workflow roept het EQTL_LOOKUP process aan dat
    een R-script uitvoert om overlappende varianten tussen
    eQTL- en GWAS-bestanden te vinden.
*/

include { EQTL_LOOKUP } from '../modules/eqtl_lookup/main.nf'

workflow EQTL_LOOKUP_WORKFLOW {

    take:
    eqtl_file
    gwas_file

    main:
    EQTL_LOOKUP(
        eqtl_file: eqtl_file,
        gwas_file: gwas_file
    )

    emit:
    result = EQTL_LOOKUP.out
}
