/**
 * Nextflow process: FUSION_TWAS
 *
 * Runs FUSION TWAS (Transcriptome-Wide Association Study) analysis for a given chromosome.
 *
 * Inputs:
 *   - tuple:
 *       - path(gwas_file):      Path to GWAS summary statistics file.
 *       - path(weights_list):   Path to FUSION weights list file.
 *       - path(weights_dir):    Directory containing FUSION weights.
 *       - path(ref_ld_chr_prefix): Prefix directory for reference LD files (per chromosome).
 *       - val(chr):             Chromosome number to analyze.
 *       - path(script):         Path to the R script for FUSION TWAS.
 *
 * Outputs:
 *   - path "fusion_twas_output" (emit: fusion_results): Output directory containing TWAS results.
 *
 * PublishDir:
 *   - Results are copied to "results/fusion/chr${chr}".
 *
 * Notes:
 *   - Sets BLAS/OpenMP environment variables to limit thread usage to allocated CPUs.
 *   - Lists available LD files for debugging.
 *   - Runs the provided R script with appropriate arguments for FUSION TWAS.
 */
process FUSION_TWAS {
    cpus 2
    tag "fusion-twas"

    input:
    tuple(
        path(gwas_file),
        path(weights_list),
        path(weights_dir),
        path(ref_ld_chr_prefix),   // Nextflow will stage this dir
        val(chr),
        path(script)
    )

    output:
    path "fusion_twas_output", emit: fusion_results

    publishDir "results/fusion/chr${chr}", mode: 'copy'

    script:
    """
    mkdir -p fusion_twas_output
    # Force BLAS/OpenMP to use exactly the number of cores allocated
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    # Debug: show LD files at the staged prefix
    echo "LD files in container at: ${ref_ld_chr_prefix}/"
    ls -lh ${ref_ld_chr_prefix}/chr${chr}.* || true

    Rscript ${script} \
      --sumstats    ${gwas_file} \
      --weights     ${weights_list} \
      --weights_dir ${weights_dir} \
      --ref_ld_chr  ${ref_ld_chr_prefix}/chr \
      --chr         ${chr} \
      --out         fusion_twas_output/chr${chr}_results.txt
    """
}
