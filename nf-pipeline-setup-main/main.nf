nextflow.enable.dsl = 2

process EQTL_LOOKUP {
    tag "eqtl-lookup"

    input:
    tuple path(eqtl_file), path(gwas_file), path(script)

    output:
    path "eqtl_lookup_result.txt", emit: eqtl_result
    publishDir "results", mode: 'copy'


    script:
    """
    Rscript ${script} \\
        ${eqtl_file} \\
        ${gwas_file} \\
        eqtl_lookup_result.txt
    """
}

workflow {
    input_ch = Channel.of([
        file(params.eqtl_file),
        file(params.gwas_file),
        file("scripts/eqtl_lookup.R")
    ])

    EQTL_LOOKUP(input_ch)
}
