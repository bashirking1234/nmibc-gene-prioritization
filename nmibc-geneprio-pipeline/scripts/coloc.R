#' coloc.R - Colocalization Analysis Script
#'
#' This script performs colocalization analysis between eQTL and GWAS summary statistics using the `coloc` R package.
#' It identifies loci with evidence for shared genetic signals between gene expression and disease/trait association.
#'
#' ## Usage
#' ```
#' eqtl[gz], gwas, outdir, eqtl_p, gwas_p, window_kb
#' ```
#' - `eqtl[gz]`      : Path to eQTL summary statistics file (TSV or GZ).
#' - `gwas`          : Path to GWAS summary statistics file.
#' - `outdir`        : Output directory for results.
#' - `eqtl_p`        : eQTL p-value threshold for lead SNP selection.
#' - `gwas_p`        : GWAS p-value threshold for lead SNP selection.
#' - `window_kb`     : Window size (in kilobases) around lead SNPs for locus definition.
#'
#' ## Workflow
#' 1. **Input Parsing**: Reads command-line arguments and checks input files.
#' 2. **Data Loading**: Loads eQTL and GWAS summary statistics using `data.table::fread`.
#' 3. **Preprocessing**:
#'    - Parses variant IDs to extract chromosome, position, reference, and alternate alleles.
#'    - Computes effect sizes, variances, and minor allele frequencies.
#'    - Estimates standard deviation of the eQTL trait.
#' 4. **Lead SNP Identification**: Selects lead SNPs based on p-value thresholds in both datasets.
#' 5. **Locus Definition**: For each lead SNP, defines a locus window and extracts overlapping variants from both datasets.
#' 6. **Colocalization Analysis**: Runs `coloc.abf` for each locus to compute posterior probabilities for five hypotheses (H0-H4).
#' 7. **Results Output**:
#'    - Writes locus-level hypothesis summary (`coloc_hypotheses_summary.txt`).
#'    - Writes SNP-level colocalization results (`coloc_snp_level.txt`).
#' 8. **Candidate Gene Extraction**:
#'    - Identifies loci with strong colocalization evidence (PP.H4 ≥ 0.8).
#'    - Extracts 95% credible sets and top candidate genes per locus.
#'    - Writes candidate gene list (`coloc_candidate_genes.txt`).
#'
#' ## Output Files
#' - `coloc_hypotheses_summary.txt` : Locus-level summary of posterior probabilities for each hypothesis.
#' - `coloc_snp_level.txt`          : SNP-level colocalization results.
#' - `coloc_candidate_genes.txt`    : Candidate genes and credible sets for loci with strong colocalization.
#' - `coloc.log`                    : Log file with progress and messages.
#'
#' ## Dependencies
#' - R packages: `data.table`, `coloc`
#'
#' ## Notes
#' - Assumes eQTL and GWAS files contain variant IDs in the format `chr:pos:ref:alt`.
#' - The eQTL sample size (`N_eqtl`) is hardcoded to 134 (bladder tissue).
#' - Designed for use with GTEx eQTL and GWAS summary statistics.
#'
#'
#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)  # fast data handling
  library(coloc)       # colocalization
})

# === Parse command-line arguments ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: coloc_gtex.R eqtl[gz] gwas outdir eqtl_p gwas_p window_kb")
}
# Assign args
eqtl_path        <- args[1]
gwas_path        <- args[2]
output_dir       <- args[3]
eqtl_pval_thresh <- as.numeric(args[4])
gwas_pval_thresh <- as.numeric(args[5])
window_kb        <- as.numeric(args[6])

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# === Logging ===
logfile <- file.path(output_dir, "coloc.log")
con     <- file(logfile, open = "wt")
sink(con, type = "output"); sink(con, type = "message")
message("=== Starting coloc_gtex.R ==="); flush.console()

# === Read & prep eQTL data ===
# GTEx cis-eQTL input (TSV or gzipped TSV)
eqtl <- fread(eqtl_path)
setnames(eqtl, "variant_id", "snp")  # unify column name
# Parse chromosome, position, alleles
eqtl[, c("chr", "pos", "ref", "alt") := tstrsplit(snp, ":", fixed = TRUE)]
eqtl[, pos := as.integer(pos)]

# Compute eQTL summary stats
N_eqtl <- 80  # fixed sample size for Bladder tissue # 134 in GTEx v8 and 80 in BCG dataset
eqtl[, `:=`(
  beta_eqtl    = as.numeric(slope),      # effect size
  varbeta_eqtl = (slope_se)^2,            # variance of effect
  MAF_eqtl     = as.numeric(af)           # minor allele frequency
)]

# corrected: explicitly invoke the unexported function
eqtl_sdY <- coloc:::sdY.est(
  vbeta = eqtl$varbeta_eqtl,
  maf   = eqtl$MAF_eqtl,
  n     = N_eqtl
)


# === Read & prep GWAS data ===
# Survival trait meta-GWAS input
gwas <- fread(gwas_path)
setnames(gwas, "MarkerName", "snp")
# Parse columns
gwas[, c("chr", "pos", "ref", "alt") := tstrsplit(snp, ":", fixed = TRUE)]
gwas[, pos := as.integer(pos)]

# Compute GWAS summary stats
gwas[, `:=`(
  beta_gwas    = as.numeric(Effect),      # log-hazard ratio treated as log-OR
  varbeta_gwas = (StdErr)^2,              # variance
  MAF_gwas     = as.numeric(MinFreq),      # minor allele freq
  N_gwas       = as.numeric(TotalSampleSize),
  cases_gwas   = as.numeric(TotalEvents)   # number of events
)]

# Compute controls per SNP
gwas[, controls_gwas := N_gwas - cases_gwas]

# === Identify lead SNPs based on p-value thresholds ===
eqtl_leads <- eqtl[pval_nominal < eqtl_pval_thresh, unique(snp)]
gwas_leads <- gwas[`P-value` < gwas_pval_thresh, unique(snp)]
lead_snps  <- unique(c(eqtl_leads, gwas_leads))

# Initialize containers for results
locus_list <- vector("list", length(lead_snps))
snp_list   <- vector("list", length(lead_snps))

# === Main loop: per lead SNP, run coloc ===
for (i in seq_along(lead_snps)) {
  lead <- lead_snps[i]
  # Extract region around lead in both datasets
  ce <- eqtl[snp == lead]
  cg <- gwas[snp == lead]
  if (nrow(ce)==0 && nrow(cg)==0) next

  # Determine genomic window
  chr    <- if (nrow(cg) > 0) cg$chr[1] else ce$chr[1]
  center <- if (nrow(cg) > 0) cg$pos[1] else ce$pos[1]
  start  <- center - window_kb * 1000
  end    <- center + window_kb * 1000

  # Subset and merge SNPs in window
  re <- unique(eqtl[chr == chr & pos >= start & pos <= end], by = "snp")
  rg <- unique(gwas[chr == chr & pos >= start & pos <= end], by = "snp")
  m  <- merge(re, rg, by = "snp")
  # Keep complete cases for essential stats
  m <- m[complete.cases(m[, .(beta_eqtl, varbeta_eqtl, beta_gwas, varbeta_gwas)]), ]
  if (nrow(m) < 3) next

  gene <- unique(m$gene_name)[1]

  # Run coloc with dataset1 = eQTL and dataset2 = GWAS
  res <- coloc.abf(
    dataset1 = list(
      snp     = m$snp,
      beta    = m$beta_eqtl,
      varbeta = m$varbeta_eqtl,
      type    = "quant",
      N       = N_eqtl,
      sdY     = eqtl_sdY
    ),
    dataset2 = list(
      snp       = m$snp,
      beta      = m$beta_gwas,
      varbeta   = m$varbeta_gwas,
      type      = "cc",
      Ncases    = m$cases_gwas,
      Ncontrols = m$controls_gwas
    )
  )

  # Collect locus-level posterior probabilities
  pp <- as.numeric(res$summary[paste0("PP.H", 0:4, ".abf")])
  locus_list[[i]] <- data.table(
    lead_snp  = lead,
    chr       = chr,
    center    = center,
    gene_name = gene,
    PP.H0     = pp[1], PP.H1 = pp[2], PP.H2 = pp[3],
    PP.H3     = pp[4], PP.H4 = pp[5]
  )

  # Collect SNP-level posterior results
  tbl <- as.data.table(res$results)
  tbl[, `:=`(
    gene_name = gene,
    lead_snp  = lead,
    chr       = chr,
    center    = center
  )]
  snp_list[[i]] <- tbl

  message(sprintf("[coloc] %3d/%3d loci done", i, length(lead_snps)))
  flush.console()
}

# === Write summary outputs ===
fwrite(
  rbindlist(locus_list)[PP.H4 >= 0.8],
  file.path(output_dir, "coloc_hypotheses_summary.txt"),
  sep = "\t"
)


fwrite(rbindlist(snp_list),   file.path(output_dir, "coloc_snp_level.txt"),         sep = "\t")
message("Wrote hypothesis and SNP-level summaries.")

# === Candidate-gene extraction with 95% credible sets ===
h4_thresh    <- 0.8
locus_summary <- rbindlist(locus_list)
good_loci     <- locus_summary[PP.H4 >= h4_thresh, lead_snp]

cand_list <- vector("list", length(good_loci))
for (j in seq_along(good_loci)) {
  lead <- good_loci[j]
  sub  <- rbindlist(snp_list)[lead_snp == lead]
  sub  <- sub[order(-SNP.PP.H4)]
  sub[, cumPP := cumsum(SNP.PP.H4)]
  cred <- sub[cumPP <= 0.95]
  if (nrow(cred) == 0) cred <- sub[1]
  top  <- cred[which.max(SNP.PP.H4)]
  cand_list[[j]] <- data.table(
    lead_snp        = lead,
    PP.H4           = locus_summary[lead_snp==lead, PP.H4],
    credible_snps   = paste(cred$snp, collapse=";"),
    top_snp         = top$snp,
    candidate_gene  = top$gene_name,
    top_snp_PP      = top$SNP.PP.H4
  )
}

candidates <- rbindlist(cand_list)
fwrite(candidates, file.path(output_dir, "coloc_candidate_genes.txt"), sep = "\t")
message("Wrote candidate gene list.")

message("=== Finished coloc_gtex.R ===")
sink(type="message"); sink(type="output"); close(con)
# EOF
