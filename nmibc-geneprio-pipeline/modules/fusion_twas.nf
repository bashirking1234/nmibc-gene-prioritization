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
