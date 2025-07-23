#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
 * Main Nextflow workflow for NMIBC gene prioritization pipeline.
 *
 * This workflow orchestrates the following steps:
 * 1. GWAS Filtering: Filters GWAS summary statistics based on a user-defined p-value threshold.
 * 2. eQTL Lookup: Optionally performs eQTL lookup using filtered GWAS results and an eQTL file.
 * 3. Colocalization Analysis: Optionally runs colocalization analysis between GWAS and eQTL data.
 * 4. Fusion TWAS: Optionally runs FUSION TWAS analysis across chromosomes 1-22.
 * 5. Gene Prioritization: Optionally prioritizes genes based on results from previous analyses.
 *
 * Modules included:
 *   - EQTL_LOOKUP: eQTL lookup process (scripts/eqtl_lookup.R)
 *   - COLOC_ANALYSIS: Colocalization analysis process (scripts/coloc.R)
 *   - FUSION_TWAS: FUSION TWAS process (scripts/fusion_twas.R)
 *   - GENE_PRIORITIZATION: Gene prioritization process (scripts/prioritize_genes.R)
 *
 * Parameters (to be provided via Nextflow config or command line):
 *   - params.gwas_file: Path to GWAS summary statistics file.
 *   - params.pvalue_threshold: P-value threshold for GWAS filtering.
 *   - params.eqtl_file: Path to eQTL file.
 *   - params.weights_list: Path to FUSION weights list.
 *   - params.weights_dir: Directory containing FUSION weights.
 *   - params.ld_dir: Directory containing LD reference panels.
 *   - params.run_eqtl: Boolean, whether to run eQTL lookup.
 *   - params.run_coloc: Boolean, whether to run colocalization analysis.
 *   - params.run_fusion: Boolean, whether to run FUSION TWAS.
 *   - params.run_prioritization: Boolean, whether to run gene prioritization.
 *   - params.eqtl_result_file: Path to eQTL lookup results (for prioritization).
 *   - params.coloc_candidates_file: Path to colocalization candidates (for prioritization).
 *   - params.fusion_results_dir: Directory with FUSION TWAS results (for prioritization).
 *
 * Output:
 *   - Filtered GWAS file and SNP count log in results/gwas/
 *   - Results from each analysis module, as configured.
 *
 * Usage:
 *   nextflow run main.nf --gwas_file <file> --pvalue_threshold <value> [other params...]
 */


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

    
    // 3) Coloc analysis (using raw GWAS)
    coloc_input = Channel
        .of(file(params.gwas_file))
        .map { rawGwas -> tuple(
            file(params.eqtl_file),
            rawGwas,
            file("scripts/coloc.R")
        ) }

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
