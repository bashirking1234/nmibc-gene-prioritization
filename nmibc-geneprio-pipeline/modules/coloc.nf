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
