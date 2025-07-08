#!/usr/bin/env Rscript

# Robust FUSION TWAS Script with Dynamic Chromosome Detection
.libPaths("/gpfs/home2/hbashir1/Rlibs")
suppressMessages({
  library(plink2R)
  library(optparse)
})

# Command-line options
option_list <- list(
  make_option("--sumstats",    type="character", help="GWAS sumstats (MarkerName,Allele1,Allele2,Effect,StdErr,TotalSampleSize)"),
  make_option("--weights",     type="character", help="Weight list file (WGT,ID,CHR,P0,P1)"),
  make_option("--weights_dir", type="character", default=NA, help="Directory for weight files"),
  make_option("--ref_ld_chr",  type="character", help="PLINK LD reference prefix"),
  make_option("--chr",         type="character", help="Chromosome number (or NA for auto)"),
  make_option("--out",         type="character", help="Output file"),
  make_option("--max_impute",  type="double",    default=0.5, help="Max missing fraction for Z-imputation"),
  make_option("--force_model", type="character", default=NA, help="Force model: blup,lasso,enet,top1"),
  make_option("--hsq_p",       type="double",    default=0.01, help="Minimum heritability p-value")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Enhanced allele QC
allele.qc <- function(a1, a2, ref1, ref2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  ref1 <- toupper(ref1); ref2 <- toupper(ref2)
  complement <- function(x) {
    ifelse(x=="A","T",
    ifelse(x=="T","A",
    ifelse(x=="C","G",
    ifelse(x=="G","C", NA))))
  }
  flip1 <- complement(ref1)
  flip2 <- complement(ref2)
  valid.alleles <- c("A","T","G","C")
  keep <- a1 %in% valid.alleles & a2 %in% valid.alleles &
          !((a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C"))
  flip <- (a1==ref2 & a2==ref1) | (a1==flip2 & a2==flip1)
  flip[is.na(flip)] <- FALSE
  keep[is.na(keep)] <- FALSE
  list(keep=keep, flip=flip)
}

# Sumstats loader (prefix chr)
sumstat_loader <- function(path) {
  raw <- read.table(path, header=TRUE, sep="\t", as.is=TRUE)
  req <- c("MarkerName","Allele1","Allele2","Effect","StdErr","TotalSampleSize")
  if (!all(req %in% colnames(raw))) {
    stop("Missing cols: ", paste(setdiff(req, colnames(raw)), collapse=", "))
  }
  parts <- strsplit(raw$MarkerName, ":")
  df <- data.frame(
    SNP  = raw$MarkerName,
    chr  = paste0("chr", vapply(parts, `[`, 1, FUN.VALUE="")),
    pos  = as.integer(vapply(parts, `[`, 2, FUN.VALUE="")),
    A1   = toupper(raw$Allele1),
    A2   = toupper(raw$Allele2),
    Z    = raw$Effect / raw$StdErr,
    N    = raw$TotalSampleSize,
    stringsAsFactors = FALSE
  )
  df$locus <- with(df, paste0(chr,":",pos,":",A1,":",A2))
  return(df)
}

# Weight list processor
get_valid_weights <- function(weight_file, weights_dir, target_chr) {
  wgtlist <- tryCatch({
    df <- read.table(weight_file, header=TRUE, stringsAsFactors=FALSE)
    if (!all(c("WGT","ID","CHR","P0","P1") %in% colnames(df))) {
      df2 <- read.table(weight_file, header=FALSE, stringsAsFactors=FALSE)
      if (ncol(df2)==5) setNames(df2, c("WGT","ID","CHR","P0","P1")) else
        stop("Weight list must have 5 columns")
    } else df
  }, error=function(e) stop("Weight list read failed: ", e$message))

  valid <- lapply(seq_len(nrow(wgtlist)), function(i) {
    wfile <- file.path(weights_dir, basename(wgtlist$WGT[i]))
    if (!file.exists(wfile)) {
      warning("Missing: ", wfile); return(NULL)
    }
    env <- new.env(); load(wfile, envir=env)
    if (!"snps" %in% ls(env)) { warning("No snps in: ", wfile); return(NULL) }
    chrs <- unique(env$snps[,1])
    if (paste0("chr", target_chr) %in% as.character(chrs)) wgtlist[i,,drop=FALSE] else NULL
  })
  res <- do.call(rbind, Filter(Negate(is.null), valid))
  if (nrow(res)==0) stop("No valid weights for chr ", target_chr)
  return(res)
}

# Model processor
model_processor <- function(wgt.mat, cv.perf, sumstat.df, bim.sub, force.model, max.imp, hsq.p) {
  model.name <- if (!is.na(force.model)) force.model else rownames(cv.perf)[which.min(cv.perf$p)]
  wgt <- wgt.mat[, model.name]

  m.idx <- match(bim.sub$locus, sumstat.df$locus)
  keep1 <- !is.na(m.idx)
  if (sum(keep1) < 2) stop("<2 SNPs overlap")
  df.idx <- m.idx[keep1]
  wgt <- wgt[keep1]
  z   <- sumstat.df$Z[df.idx]

  aqc <- allele.qc(a1=sumstat.df$A1[df.idx], a2=sumstat.df$A2[df.idx],
                   ref1=bim.sub$ref[keep1], ref2=bim.sub$alt[keep1])
  z[aqc$flip] <- -z[aqc$flip]
  keep2 <- aqc$keep
  wgt <- wgt[keep2]; z <- z[keep2]

  bed.sub <- genos$bed[, bim.sub$index[keep2], drop=FALSE]
  R <- cor(bed.sub, use="pairwise.complete.obs")

  miss.frac <- sum(is.na(wgt)) / length(wgt)
  if (miss.frac > max.imp) stop("Too many missing weights: ", round(miss.frac,3))

  denom <- sqrt(as.numeric(t(wgt) %*% R %*% wgt))
  twasZ <- as.numeric((wgt %*% z) / denom)
  pval  <- 2 * pnorm(-abs(twasZ))

  hsq <- cv.perf$hsq[model.name]
  if (cv.perf$p[model.name] > hsq.p) hsq <- NA

  list(Z=twasZ, P=pval, model=model.name, hsq=hsq,
       nsnp=length(z), nwgt=sum(!is.na(wgt)))
}

# Main pipeline
sumstat <- sumstat_loader(opt$sumstats)
if (is.na(opt$chr)) cur_chrs <- unique(sub("chr", "", sumstat$chr)) else cur_chrs <- opt$chr

# Load reference once
genos <- tryCatch({ read_plink(opt$ref_ld_chr, impute="avg") },
                 error=function(e) stop("LD load failed: ", e$message))
colnames(genos$bim) <- c("chr","id","cM","pos","ref","alt")
genos$bim$chr <- paste0("chr", sub(".*?([0-9XYM]+)$", "\1", genos$bim$chr))
genos$bim$locus <- with(genos$bim, paste0(chr,":",pos,":",ref,":",alt))

final.list <- list()
for (chr in cur_chrs) {
  message("--- Processing chr ", chr, " ---")
  wlist <- get_valid_weights(opt$weights, opt$weights_dir, chr)
  out <- data.frame(FILE=wlist$WGT, ID=wlist$ID, CHR=wlist$CHR,
                    P0=wlist$P0, P1=wlist$P1, HSQ=NA, NSNP=NA,
                    NWGT=NA, TWAS.Z=NA, TWAS.P=NA, MODEL=NA,
                    stringsAsFactors=FALSE)

  for (i in seq_len(nrow(wlist))) {
    win <- c(wlist$P0[i], wlist$P1[i])
    idx <- which(genos$bim$pos >= win[1] & genos$bim$pos <= win[2])
    if (length(idx)==0) { warning("No SNPs in window for gene ", wlist$ID[i]); next }
    bim_sub <- genos$bim[idx, ]
    bim_sub$index <- idx

    ss.df <- subset(sumstat, chr == paste0("chr", chr) & pos >= win[1] & pos <= win[2])
    if (nrow(ss.df)==0) { warning("No sumstats in window for gene ", wlist$ID[i]); next }
    ss.df$locus <- with(ss.df, paste0(chr,":",pos,":",A1,":",A2))

    env <- new.env(); load(file.path(opt$weights_dir, basename(wlist$WGT[i])), envir=env)
    res <- tryCatch({
      model_processor(env$wgt.matrix, env$cv.performance,
                      ss.df, bim_sub, opt$force_model,
                      opt$max_impute, opt$hsq_p)
    }, error=function(e) { warning("Gene ", wlist$ID[i], " failed: ", e$message); NULL })
    if (!is.null(res)) {
      out[i, c("TWAS.Z","TWAS.P","MODEL","HSQ","NSNP","NWGT")] <-
        c(res$Z, res$P, res$model, res$hsq, res$nsnp, res$nwgt)
      message(sprintf("%s: Z=%.2f P=%.2e", wlist$ID[i], res$Z, res$P))
    }
  }
  final.list[[as.character(chr)]] <- out
}
all.res <- do.call(rbind, final.list)
valid <- subset(all.res, !is.na(TWAS.Z))
message(sprintf("Done: %d/%d genes succeeded (%.1f%%)", nrow(valid), nrow(all.res), 100*nrow(valid)/nrow(all.res)))
write.table(valid, opt$out, sep="\t", quote=FALSE, row.names=FALSE)
