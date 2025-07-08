#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(coloc)
})

# === Parse args ===
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 6) 
  stop("Usage: coloc_gtex.R eqtl gwas outdir eqtl_p gwas_p window_kb")
eqtl_path        <- args[1]
gwas_path        <- args[2]
output_dir       <- args[3]
eqtl_pval_thresh <- as.numeric(args[4])
gwas_pval_thresh <- as.numeric(args[5])
window_kb        <- as.numeric(args[6])

if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# === Logging ===
logfile <- file.path(output_dir, "coloc.log")
con     <- file(logfile, open="wt")
sink(con, type="output"); sink(con, type="message")
message("=== Starting coloc_gtex.R ==="); flush.console()

# === Read & prep ===
eqtl <- fread(eqtl_path); setnames(eqtl,"variant_id","snp")
gwas <- fread(gwas_path); setnames(gwas,"MarkerName","snp")
for(dt in list(eqtl,gwas)){
  dt[, c("chr","pos","ref","alt") := tstrsplit(snp, ":", fixed=TRUE)]
  dt[, pos := as.integer(pos)]
}
# lead SNPs
eqtl_leads <- eqtl[pval_nominal < eqtl_pval_thresh, unique(snp)]
gwas_leads <- gwas[`P-value`    < gwas_pval_thresh, unique(snp)]
lead_snps  <- unique(c(eqtl_leads, gwas_leads))

# precompute stats
eqtl[, `:=`(
  beta_eqtl    = as.numeric(slope),
  varbeta_eqtl = (slope_se)^2,
  N_eqtl       = as.numeric(num_var)
)]
gwas[, `:=`(
  beta_gwas    = as.numeric(Effect),
  varbeta_gwas = (StdErr)^2,
  N_gwas       = as.numeric(TotalSampleSize),
  cases_gwas   = as.numeric(TotalEvents)
)]
gwas[, controls_gwas := N_gwas - cases_gwas]

# prepare collectors
locus_list <- vector("list", length(lead_snps))
snp_list   <- vector("list", length(lead_snps))

# === Loop over loci ===
for(i in seq_along(lead_snps)) {
  lead   <- lead_snps[i]
  ce     <- eqtl[snp==lead]
  cg     <- gwas[snp==lead]
  if(nrow(ce)==0 && nrow(cg)==0) next

  # window boundaries
  chr    <- if(nrow(cg)>0) cg$chr[1] else ce$chr[1]
  center <- if(nrow(cg)>0) cg$pos[1] else ce$pos[1]
  start  <- center - window_kb*1000
  end    <- center + window_kb*1000

  # subset & dedupe
  re <- unique(eqtl[chr==chr & pos>=start & pos<=end], by="snp")
  rg <- unique(gwas[chr==chr & pos>=start & pos<=end], by="snp")
  m  <- merge(re, rg, by="snp")
  m  <- m[complete.cases(m[, .(beta_eqtl,varbeta_eqtl,beta_gwas,varbeta_gwas)]), ]
  if(nrow(m) < 3) next

  # pick the gene for annotation
  gene <- unique(m$gene_name)[1]

  # run coloc
  res <- coloc.abf(
    dataset1 = list(
      snp     = m$snp, beta    = m$beta_eqtl,
      varbeta = m$varbeta_eqtl, N = m$N_eqtl,
      type    = "quant", sdY = sd(m$beta_eqtl)
    ),
    dataset2 = list(
      snp     = m$snp, beta    = m$beta_gwas,
      varbeta = m$varbeta_gwas, N = m$N_gwas,
      s       = sum(m$cases_gwas)/(sum(m$cases_gwas)+sum(m$controls_gwas)),
      type    = "cc"
    )
  )

  # --- locus‐level summary ---
  pp_all <- res$summary[grep("^PP\\.H", names(res$summary))]
  pp     <- as.numeric(pp_all[
    c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")])
  locus_list[[i]] <- data.table(
    lead_snp  = lead,
    chr       = chr,
    center    = center,
    gene_name = gene,
    PP.H0     = pp[1],
    PP.H1     = pp[2],
    PP.H2     = pp[3],
    PP.H3     = pp[4],
    PP.H4     = pp[5]
  )

  # --- SNP‐level posterior table ---
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

# write locus‐level summary
locus_summary <- rbindlist(locus_list, use.names=TRUE, fill=TRUE)
fwrite(locus_summary,
       file.path(output_dir,"coloc_hypotheses_summary.txt"),
       sep="\t", quote=FALSE)
message("Wrote locus summary.")

# write SNP‐level results
snp_tbl <- rbindlist(snp_list, use.names=TRUE, fill=TRUE)
fwrite(snp_tbl,
       file.path(output_dir,"coloc_snp_level.txt"),
       sep="\t", quote=FALSE)
message("Wrote SNP‐level results.")

# === candidate‐gene extraction ===
h4_thresh <- 0.8
good_loci <- locus_summary[PP.H4 >= h4_thresh, lead_snp]
cand_list <- vector("list", length(good_loci))

for(j in seq_along(good_loci)) {
  lead <- good_loci[j]
  sub  <- snp_tbl[lead_snp==lead]
  o    <- order(sub$SNP.PP.H4, decreasing=TRUE)
  cump <- cumsum(sub$SNP.PP.H4[o])
  w    <- which(cump >= 0.95)[1]
  cred <- sub[o[1:w]]
  top  <- cred[SNP.PP.H4 == max(SNP.PP.H4)][, unique(gene_name)]
  cand_list[[j]] <- data.table(
    lead_snp        = lead,
    PP.H4           = locus_summary[lead_snp==lead, PP.H4],
    credible_snps   = paste(cred$snp, collapse=";"),
    candidate_genes = paste(top, collapse=";")
  )
}

candidates <- rbindlist(cand_list, use.names=TRUE, fill=TRUE)
fwrite(candidates,
       file.path(output_dir,"coloc_candidate_genes.txt"),
       sep="\t", quote=FALSE)
message("Wrote candidate gene list.")

message("=== Finished coloc_gtex.R ===")
flush.console()

# restore sinks
sink(type="message"); sink(type="output"); close(con)
# EOF
