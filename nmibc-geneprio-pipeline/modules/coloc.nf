/*
 * COLOC_ANALYSIS process
 *
 * This Nextflow process performs colocalization analysis between eQTL and GWAS datasets using an R script.
 *
 * Inputs:
 *   - eqtl_file: Path to the eQTL summary statistics file.
 *   - gwas_file: Path to the GWAS summary statistics file.
 *   - script: Path to the R script that performs the colocalization analysis.
 *
 * Outputs:
 *   - coloc.log: Log file capturing the run details.
 *   - coloc_hypotheses_summary.txt: Locus-level summary of posterior probabilities for hypotheses H0–H4.
 *   - coloc_snp_level.txt: SNP-level posterior probability table.
 *   - coloc_candidate_genes.txt: List of candidate genes identified by the analysis.
 *
 * Additional Information:
 *   - The process sets environment variables to limit BLAS and OpenMP threads to 1 for reproducibility.
 *   - Output files are published to the '/home/hbashir1/nmibc-geneprio/results/coloc' directory.
 *   - The R script is executed with the following arguments:
 *       1. eQTL file path
 *       2. GWAS file path
 *       3. Output directory (current directory)
 *       4. p-value threshold for eQTLs (1e-6). Can be changed to a different value if needed.
 *       5. p-value threshold for GWAS (1e-6). Can be changed to a different value if needed.
 *       6. Maximum number of SNPs to consider (500). Can be changed to a different value if needed.
 */
process COLOC_ANALYSIS {
    tag "coloc-analysis"
    cpus 1

    input:
    tuple path(eqtl_file), path(gwas_file), path(script)


    output:
    // your run‐log
    path 'coloc.log',                       emit: coloc_log
    // locus‐level summary of H0–H4
    path 'coloc_hypotheses_summary.txt',    emit: coloc_hypotheses_summary
    // full SNP‐level posterior table
    path 'coloc_snp_level.txt',             emit: coloc_snp_level
    // final candidate‐gene list
    path 'coloc_candidate_genes.txt',       emit: coloc_candidates

    publishDir "/home/hbashir1/nmibc-geneprio/results/coloc", mode: 'copy'

    script:
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript ${script} \\
        ${eqtl_file} \\
        ${gwas_file} \\
        .            \\
        1e-6         \\
        1e-6         \\
        500
    """
}
