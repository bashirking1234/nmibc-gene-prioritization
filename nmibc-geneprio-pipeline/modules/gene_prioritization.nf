#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**
 * Module: GENE_PRIORITIZATION
 *
 * Count & rank genes by how often they appear in:
 *   1) eQTL lookup results (gene_name)
 *   2) Coloc candidate genes (semicolon-separated)
 *   3) Fusion TWAS results (ID prefix before '_')
 *
 * Inputs:
 *   A single channel emitting tuples of four paths:
 *     - eqtl_file (eQTL lookup TSV)
 *     - coloc_file (Coloc candidate genes TSV)
 *     - fusion_dir (directory of per-chr TWAS outputs)
 *     - script (R script to run prioritization)
 *
 * Outputs:
 *   gene_priority.tsv (ranked gene counts)
 */

process GENE_PRIORITIZATION {
    tag "gene-prioritization"
    cpus 1

    publishDir "results/gene_prioritization", mode: 'copy'

    input:
    tuple path(eqtl_file), path(coloc_file), path(fusion_dir), path(script)

    output:
    path "gene_priority.tsv", emit: prioritized_genes

    script:
    """
    Rscript ${script} \
        ${eqtl_file} \
        ${coloc_file} \
        ${fusion_dir} \
        gene_priority.tsv
    """
}
