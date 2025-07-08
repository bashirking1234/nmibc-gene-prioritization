#!/usr/bin/env Rscript

# ------------------------------
# Robust FUSION TWAS Script on hg38 (full 18‐column output)
# ------------------------------
.libPaths("/gpfs/home2/hbashir1/Rlibs")
suppressMessages({
  library(plink2R)
  library(optparse)
})

# Command-line options
opt_list <- list(
  make_option("--sumstats",    type="character", help="GWAS sumstats [required]"),
  make_option("--weights",     type="character", help="Weight list file [required]"),
  make_option("--weights_dir", type="character", help="Directory for weight files"),
  make_option("--ref_ld_chr",  type="character", help="PLINK LD reference prefix [required]"),
  make_option("--chr",         type="character", help="Chromosome to analyze [required]"),
  make_option("--out",         type="character", help="Output file [required]"),
  make_option("--hsq_p",       type="double",    default=0.01, help="Min hsq p-value")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# Allele‐QC helper
allele.qc <- function(a1,a2,ref1,ref2) {
  a1<-toupper(a1); a2<-toupper(a2)
  ref1<-toupper(ref1); ref2<-toupper(ref2)
  flip_ref <- function(x) { f<-x; f[x=="A"]<-"T"; f[x=="T"]<-"A"; f[x=="C"]<-"G"; f[x=="G"]<-"C"; f }
  flip1 <- flip_ref(ref1); flip2 <- flip_ref(ref2)
  keep <- a1%in%c("A","C","G","T") & a2%in%c("A","C","G","T") &
          !( (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|
             (a1=="C"&a2=="G")|(a1=="G"&a2=="C") )
  flip <- (a1==ref2 & a2==ref1) | (a1==flip2 & a2==flip1)
  list(keep=keep, flip=flip)
}

# 1) LOAD GWAS sumstats & filter to chromosome
raw <- read.table(opt$sumstats, header=TRUE, sep="\t", as.is=TRUE)
raw$MarkerName <- trimws(raw$MarkerName)
raw$chr <- sub(":.*","", raw$MarkerName)
raw_chr <- subset(raw, chr==paste0("chr",opt$chr))
sumstat <- data.frame(
  locus = trimws(raw_chr$MarkerName),
  A1    = toupper(raw_chr$Allele1),
  A2    = toupper(raw_chr$Allele2),
  Z     = raw_chr$Effect / raw_chr$StdErr,
  stringsAsFactors=FALSE
)

# 2) LOAD LD reference (hg38)
ld_pref <- paste0(opt$ref_ld_chr, opt$chr)
ld <- read_plink(ld_pref, impute="avg")
colnames(ld$bim) <- c("chr","locus_id","cM","pos","ref","alt")
ld$bim$locus <- with(ld$bim, paste0(chr, ":", pos, ":", ref, ":", alt))
ld_mat <- as.matrix(ld$bed)

# 3) LOAD weight list
wgtlist <- read.table(opt$weights, header=TRUE, stringsAsFactors=FALSE)
wgtlist <- subset(wgtlist, CHR==opt$chr)

# Prepare output with 18 columns
out <- data.frame(
  FILE         = wgtlist$WGT,
  ID           = wgtlist$ID,
  CHR          = wgtlist$CHR,
  P0           = wgtlist$P0,
  P1           = wgtlist$P1,
  HSQ          = NA_real_,
  BEST.GWAS.ID = NA_character_,
  BEST.GWAS.Z  = NA_real_,
  EQTL.ID      = NA_character_,
  EQTL.R2      = NA_real_,
  EQTL.Z       = NA_real_,
  EQTL.GWAS.Z  = NA_real_,
  NSNP         = NA_integer_,
  MODEL        = NA_character_,
  MODELCV.R2   = NA_real_,
  MODELCV.PV   = NA_real_,
  TWAS.Z       = NA_real_,
  TWAS.P       = NA_real_,
  stringsAsFactors=FALSE
)

# 4) PER‐GENE LOOP
for(i in seq_len(nrow(wgtlist))) {
  gid   <- wgtlist$ID[i]
  wfile <- file.path(opt$weights_dir, wgtlist$WGT[i])
  if (!file.exists(wfile)) next

  env <- new.env()
  load(wfile, envir=env)
  if (!"snps" %in% ls(env)) next
  snps <- env$snps
  # use hg38 id
  snps$id <- snps$id_hg38

  w_locus <- snps$id
  wmat    <- env$wgt.matrix
  cvp     <- env$cv.performance

  # dynamic window ±500kb around gene bounds
  orig_P0 <- wgtlist$P0[i]
  orig_P1 <- wgtlist$P1[i]
  P0 <- max(1, orig_P0 - 500000)
  P1 <- orig_P1 + 500000

  # subset LD panel
  sel   <- which(ld$bim$pos >= P0 & ld$bim$pos <= P1)
  if (length(sel)==0) next
  bim_sub <- ld$bim[sel,]
  bed_sub <- ld_mat[,sel,drop=FALSE]

  # use full GWAS (sumstat)
  sumstat_sub <- sumstat

  # match weights → LD by chr:pos
  w_coord <- sub("(:[^:]+){2}$","", w_locus)
  ld_coord <- paste0(bim_sub$chr, ":", bim_sub$pos)
  m_pos <- match(w_coord, ld_coord)
  n_pos <- sum(!is.na(m_pos))
  if (n_pos < 1) {
    out[i, "NSNP"] <- 0
    next
  }
  valid <- which(!is.na(m_pos))
  w2    <- wmat[valid,,drop=FALSE]
  curb  <- bim_sub[m_pos[valid],]
  geno2 <- bed_sub[,m_pos[valid],drop=FALSE]
  w_coord <- ld_coord[m_pos[valid]]

  # match LD → GWAS by chr:pos
  g_coord <- sub("(:[^:]+){2}$","", sumstat_sub$locus)
  m2      <- match(w_coord, g_coord)
  if (sum(!is.na(m2)) < 1) {
    out[i, "NSNP"] <- length(valid)
    next
  }
  z_all  <- sumstat_sub$Z[m2]
  a1_all <- sumstat_sub$A1[m2]
  a2_all <- sumstat_sub$A2[m2]
  ref_all <- curb$ref
  alt_all <- curb$alt

  # allele‐QC
  qc   <- allele.qc(a1_all,a2_all,ref_all,alt_all)
  flip <- qc$flip; flip[is.na(flip)] <- FALSE
  keep <- qc$keep; keep[is.na(keep)] <- FALSE

  z_all[flip]  <- -z_all[flip]
  a1_all[flip] <- ref_all[flip]
  a2_all[flip] <- alt_all[flip]

  w2    <- w2[keep,,drop=FALSE]
  geno2 <- geno2[,keep,drop=FALSE]
  z_all <- z_all[keep]

  if (length(z_all) < 1) {
    out[i, "NSNP"] <- sum(keep)
    next
  }

  # TWAS stat
  geno2_s <- scale(geno2)
  LDmat   <- crossprod(geno2_s) / (nrow(geno2_s)-1)

  # choose best model
  fm <- if ("force_model" %in% names(opt)) opt$force_model else NA_character_
  if (!is.character(fm) || length(fm)==0 || is.na(fm) || fm=="") {
    pidx <- grep("pval",rownames(cvp)); if (!length(pidx)) pidx <- grep("^p",rownames(cvp))
    best <- which.min( apply(cvp[pidx,,drop=FALSE],2,min,na.rm=TRUE) )
  } else {
    best <- which(colnames(w2)==fm)
  }

  wv    <- w2[,best]
  denom <- sqrt(as.numeric(t(wv)%*%LDmat%*%wv))
  tZ    <- as.numeric((wv %*% z_all)/denom)
  tP    <- 2*pnorm(-abs(tZ))

  # heritability
  model_hsq    <- env$hsq[1]
  model_hsq_pv <- env$hsq.pv
  hsq <- if (!is.na(model_hsq_pv) && model_hsq_pv <= opt$hsq_p) model_hsq else NA_real_

  # best GWAS SNP
  loci_win     <- sumstat_sub$locus[m2]
  best_g_idx   <- which.max(abs(z_all))
  best_gwas_id <- loci_win[best_g_idx]
  best_gwas_z  <- z_all[best_g_idx]

  # top1 eQTL
  top1_w      <- w2[,"top1"]
  eqtl_idx    <- which.max(abs(top1_w))
  eqtl_id     <- w_coord[eqtl_idx]
  eqtl_r2     <- cvp["rsq","top1"]
  eqtl_z      <- top1_w[eqtl_idx]
  eqtl_gwas_z <- z_all[eqtl_idx]

  # model CV stats
  model       <- colnames(w2)[best]
  model_cv_r2 <- cvp["rsq", best]
  model_cv_pv <- cvp["pval",best]

  # fill output
  out[i,] <- list(
    wgtlist$WGT[i], wgtlist$ID[i], wgtlist$CHR[i],
    wgtlist$P0[i], wgtlist$P1[i],
    hsq,
    best_gwas_id, best_gwas_z,
    eqtl_id, eqtl_r2, eqtl_z, eqtl_gwas_z,
    nrow(LDmat),
    model, model_cv_r2, model_cv_pv,
    tZ, tP
  )
}

# WRITE RESULTS
write.table(out, opt$out, sep="\t", quote=FALSE, row.names=FALSE)
# End of script