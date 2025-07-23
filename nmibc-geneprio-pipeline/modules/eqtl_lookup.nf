/**
 * Process: EQTL_LOOKUP
 *
 * Description:
 *   Performs an eQTL lookup by running a specified R script with provided eQTL and GWAS files.
 *   The result is saved as 'eqtl_lookup_result.txt' and published to the results directory.
 *
 * Inputs:
 *   tuple path(eqtl_file)  - Path to the eQTL file to be used in the lookup.
 *   tuple path(gwas_file)  - Path to the GWAS file to be used in the lookup.
 *   tuple path(script)     - Path to the R script that performs the lookup.
 *
 * Outputs:
 *   path "eqtl_lookup_result.txt" (emit: eqtl_result) - The output file containing the lookup results.
 *
 * Resources:
 *   cpus: 1
 *
 * Tags:
 *   "eqtl-lookup"
 *
 * PublishDir:
 *   "results/eqtl_lookup" (mode: 'copy') - Publishes process output to this directory.
 *   "results" (mode: 'copy')             - Also publishes output to this directory.
 *
 * Script:
 *   - Sets environment variables to limit BLAS and OpenMP threads to 1.
 *   - Executes the R script with the provided input files and output filename.
 */
process EQTL_LOOKUP {
    cpus 1
    tag "eqtl-lookup"
      
    publishDir "results/eqtl_lookup", mode: 'copy'

    input:
    tuple path(eqtl_file), path(gwas_file), path(script)

    output:
    path "eqtl_lookup_result.txt", emit: eqtl_result

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
