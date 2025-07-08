#!/usr/bin/env Rscript
# -------------------------------------------------------------------
# Count & rank genes by how often they appear in eQTL, Coloc & TWAS
# -------------------------------------------------------------------
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript prioritize_genes.R <eqtl_file> <coloc_file> <fusion_dir> <output_file>")
}
eqtl_file   <- args[1]
coloc_file  <- args[2]
fusion_dir  <- args[3]
out_file    <- args[4]

# 1) eQTL genes: use gene_name column
eqtl <- fread(eqtl_file, sep="\t", header=TRUE)
eqtl_genes <- eqtl$gene_name

# 2) Coloc candidate genes (may be semicolon–separated)
coloc <- fread(coloc_file, sep="\t", header=TRUE)
coloc_genes <- unlist(strsplit(coloc$candidate_genes, ";"))

# 3) Fusion TWAS genes: extract prefix of ID column before underscore
fusion_files <- list.files(path=fusion_dir,
                           pattern="*_results.txt$",
                           recursive=TRUE, full.names=TRUE)
fusion_list <- lapply(fusion_files, function(f){
  dt <- fread(f, sep="\t", header=TRUE)
  # split ID on underscore and keep prefix
  sapply(dt$ID, function(x) sub("_.*$", "", x))
})
fusion_genes <- unlist(fusion_list)

# 4) Build count table
all_genes <- unique(c(eqtl_genes, coloc_genes, fusion_genes))
dt <- data.table(gene_id = all_genes)
dt[, eqtl_count   := sum(gene_id == eqtl_genes),   by=gene_id]
dt[, coloc_count  := sum(gene_id == coloc_genes),  by=gene_id]
dt[, fusion_count := sum(gene_id == fusion_genes), by=gene_id]
dt[, total_count  := eqtl_count + coloc_count + fusion_count]

# 5) Rank and order
setorder(dt, -total_count)
dt[, rank := .I]

# 6) Write out
fwrite(dt, out_file, sep="\t", quote=FALSE)
