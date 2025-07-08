process EQTL_LOOKUP {
    cpus 1
    tag "eqtl-lookup"
      
    publishDir "results/eqtl_lookup", mode: 'copy'

    input:
    tuple path(eqtl_file), path(gwas_file), path(script)

    output:
    path "eqtl_lookup_result.txt", emit: eqtl_result
    publishDir "results", mode: 'copy'

    script:
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    Rscript ${script} \\
        ${eqtl_file} \\
        ${gwas_file} \\
        eqtl_lookup_result.txt
    """
}
