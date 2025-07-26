# NMIBC Gene Prioritization Pipeline

A Nextflow-driven workflow for integrating eQTL-lookup, colocalization, and Transcriptome-Wide-Association Study analyses to prioritize candidate genes in non–muscle invasive bladder cancer (NMIBC). Designed for execution on an HPC cluster (e.g., Snellius) using Singularity containers.

---

## Table of Contents

1. [Project Structure](#project-structure)
2. [Features](#features)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Configuration](#configuration)
6. [Workflow Components](#workflow-components)

   * [scripts/](#scripts)
   * [modules/](#modules)
   * [main.nf](#mainnf)
   * [nextflow.config](#nextflowconfig)
7. [Execution/Usage](#execution)
8. [Output Files](#output-files)
9. [Log Files](#log-files)
10. [License](#license)

---

## Project Structure

```
├── containers/          # Singularity images
├── Data/                # Custom input data
├── scripts/             # R scripts for analyses
├── modules/             # Nextflow DSL2 modules per analysis step
├── main.nf              # Main workflow definition
├── nextflow.config      # Parameters, container, and executor settings
├── README.md            # This file
├── results/             # Workflow outputs (copied)
└── work/                # Nextflow work directory (intermediate)
```

---

## Features

* **GWAS filtering** by p-value (customizable threshold)
* **eQTL lookup**: harmonization and matching between GWAS and eQTL variants
* **COLOC analysis** for colocalization of GWAS and eQTL signals
* **FUSION TWAS** for transcriptome-wide association study
* **Gene prioritization** by integrating results and ranking them across methods
* **Reproducibility** via Singularity containers and Nextflow DSL2

---

## Requirements

* [Nextflow](https://nextflow.io) (v21+)
* Java (OpenJDK 11 or higher)
* Singularity/Apptainer
* R (>=4.0) with packages: `data.table`, `coloc`, `optparse`, `plink2R`

---

## Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/<username>/nmibc-geneprio.git
   cd nmibc-geneprio
   ```
2. **Ensure Nextflow and Singularity are available** on your HPC
3. **Adjust `nextflow.config`** to reflect your data and container paths

---

## Configuration

All input paths, parameters, and container settings are defined in `nextflow.config`:

```groovy
params {
  gwas_file           = 'path/to/file'
  eqtl_file           = 'path/to/file'
  weights_dir         = 'path/to/file'
  weights_list        = 'path/to/file'
  ld_dir              = 'path/to/file'
  pvalue_threshold    = 0.1
}

singularity {
  enabled    = true
  engine     = 'apptainer'   // force Apptainer if both are installed
  autoMounts = true
}

executor {
  queueSize = 2
}
```

> Update the paths above to match your local HPC directories.

---

## Workflow Components

### scripts/

R scripts called by each module:

* **eqtl\_lookup.R**: Harmonizes and matches GWAS and eQTL variants; outputs `eqtl_lookup_result.txt`.
* **coloc.R**: Runs `coloc.abf` per locus; generates `coloc_hypotheses_summary.txt`, `coloc_snp_level.txt`, and `coloc_candidate_genes.txt`.
* **fusion\_twas.R**: Executes FUSION TWAS per gene; outputs per-chromosome results in `fusion_twas_output/chr{chr}_results.txt`.
* **prioritize\_genes.R**: Counts and ranks genes across eQTL, COLOC, and TWAS; writes `gene_priority.tsv`.

### modules/

Nextflow DSL2 modules wrapping each analysis step:

* **eqtl\_lookup.nf** → `EQTL_LOOKUP`
* **coloc.nf**       → `COLOC_ANALYSIS`
* **fusion\_twas.nf** → `FUSION_TWAS`
* **gene\_prioritization.nf** → `GENE_PRIORITIZATION`

Each module defines inputs (data files and script), executes the R script, and publishes results to `results/`.

### main.nf

Orchestrates the workflow steps:

1. **FILTER\_GWAS**: Filters GWAS by p-value threshold and writes `filtered_gwas.tbl`.
2. **EQTL\_LOOKUP**: Performs eQTL lookup on the filtered GWAS results.
3. **COLOC\_ANALYSIS**: Runs colocalization between GWAS and eQTL.
4. **FUSION\_TWAS**: Conducts transcriptome-wide association per chromosome.
5. **GENE\_PRIORITIZATION**: Integrates outputs to rank candidate genes.

Conditional execution of each step is controlled via `params.run_eqtl`, `params.run_coloc`, etc.

### nextflow\.config

See [Configuration](#configuration).

---

## Execution/Usage

Run the full pipeline:

```bash
 nextflow run main.nf --profile singularity
```

Run one module of the pipeline

```bash
 nextflow run main.nf --profile singularity --run_eqtl true
```

Run the pipeline with selective modules (true == run or false == don't run) and custom thresholds by passing flags to `nextflow run`:

```bash
# Only TWAS & prioritization
nextflow run main.nf \
  --run_eqtl false \
  --run_coloc false \
  --run_fusion true \
  --run_prioritization true \
  --pvalue_threshold 5e-8

# Full pipeline + prioritization
nextflow run main.nf \
  --run_eqtl true \
  --run_coloc true \
  --run_fusion true \
  --run_prioritization true \
  --pvalue_threshold 1e-6

# Only gene prioritization (using existing results)
nextflow run main.nf \
  --run_eqtl false \
  --run_coloc false \
  --run_fusion false \
  --run_prioritization true

# Dry-run without prioritization
nextflow run main.nf \
  --run_prioritization false
```

Only specific steps can be toggled in `nextflow.config`:

```groovy
params {
  run_eqtl           = true
  run_coloc          = true
  run_fusion         = true
  run_prioritization = true
}
```

---

## Output Files

Workflow outputs are organized under `results/`:

* **GWAS filtering**: `results/gwas/filtered_gwas.tbl`, `snp_count.log`
* **eQTL lookup**: `results/eqtl_lookup/eqtl_lookup_result.txt`
* **COLOC**: `results/coloc/coloc_hypotheses_summary.txt`, `coloc_snp_level.txt`, `coloc_candidate_genes.txt`, `coloc.log`
* **FUSION TWAS**: `results/fusion/chr{chr}/fusion_twas_output/chr{chr}_results.txt`
* **Gene Prioritization**: `results/gene_prioritization/gene_priority.tsv`

---

## Log Files

* **Nextflow work directory**: `work/` (intermediate files)
* **Console logs**: standard output and error

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
