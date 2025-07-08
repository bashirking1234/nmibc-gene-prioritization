#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EQTL_LOOKUP        } from './modules/eqtl_lookup.nf'
include { COLOC_ANALYSIS      } from './modules/coloc.nf'
include { FUSION_TWAS         } from './modules/fusion_twas.nf'
include { GENE_PRIORITIZATION } from './modules/gene_prioritization.nf'

process FILTER_GWAS {
    publishDir "results/gwas", mode: 'copy'

    input:
      tuple path(gwas_file), val(threshold)

    output:
      tuple path("filtered_gwas.tbl"), path("snp_count.log")

    script:
    """
    awk 'NR==1 || \$10 < ${threshold}' ${gwas_file} > filtered_gwas.tbl
    echo "Aantal overgebleven SNPs: \$(wc -l < filtered_gwas.tbl)" > snp_count.log
    """
}

workflow {

    println "P-value threshold = ${params.pvalue_threshold}"

    // Scripts
    def eqtl_lookup_script = file("scripts/eqtl_lookup.R")
    def coloc_script = file("${baseDir}/scripts/coloc.R")


    // 1) GWAS filtering
    Channel
        .of(file(params.gwas_file))
        .map { gwas -> tuple(gwas, params.pvalue_threshold) }
        .set { gwas_with_threshold }

    filtered_output = FILTER_GWAS(gwas_with_threshold)

    filtered_gwas_tbl_ch = filtered_output.map { gwas, log -> gwas }
    snp_log_ch           = filtered_output.map { gwas, log -> log }

    // 2) eQTL lookup
    eqtl_lookup_input = filtered_gwas_tbl_ch.map { gwas ->
        tuple(file(params.eqtl_file), gwas, eqtl_lookup_script)
    }

    

    // 3) Coloc analysis
    coloc_input = Channel
        .of(file(params.gwas_file))
        .map { raw -> tuple(file(params.eqtl_file), raw, coloc_script) }

    // 4) Fusion TWAS
    fusion_input = Channel
        .from(1..22)
        .map { chr ->
            tuple(
                file(params.gwas_file),
                file(params.weights_list),
                file(params.weights_dir),
                file(params.ld_dir),
                chr,
                file("scripts/fusion_twas.R")
            )
        }

    // 5) Conditional runs
    if ( params.run_eqtl ) {
        EQTL_LOOKUP(eqtl_lookup_input)
    }

    if ( params.run_coloc ) {
        COLOC_ANALYSIS(coloc_input)
    }

    if ( params.run_fusion ) {
        FUSION_TWAS(fusion_input)
    }

        // 6) Gene prioritization (after analyses)
    if ( params.run_prioritization ) {
        prio_in = Channel.of( tuple(
            file(params.eqtl_result_file),
            file(params.coloc_candidates_file),
            file(params.fusion_results_dir),
            file("scripts/prioritize_genes.R")
        ) )

        // capture the module invocation in a variable
        def prio_ch = GENE_PRIORITIZATION(prio_in)

        // now view the output channel
        prio_ch.prioritized_genes.view()
    }

    // nf-set up, email notification, and other configurations

    
}
// End of file: nmibc-geneprio/main.nf 
