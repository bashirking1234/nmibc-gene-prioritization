process EQTL_LOOKUP {

    tag "eqtl-lookup"

    input:
    path eqtl_file
    path gwas_file

    output:
    path "eqtl_lookup_result.txt"

    script:
    """
    Rscript scripts/eqtl_lookup.R \
        ${eqtl_file} \
        ${gwas_file} \
        eqtl_lookup_result.txt
    """
}
