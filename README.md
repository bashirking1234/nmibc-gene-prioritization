# nmibc-geneprio

**nmibc-geneprio** is a modular and reproducible Nextflow pipeline to identify and prioritize candidate genes associated with prognosis in non-muscle invasive bladder cancer (NMIBC).  
It integrates eQTL-lookup, colocalization, and TWAS and gene ranking.


> ⚠️ This pipeline is built using the [nf-core pipeline template](https://nf-co.re) but is **not an official nf-core pipeline**.

---

## 🔍 Pipeline Overview

### Goals:
- Identify SNPs associated with NMIBC prognosis
- Link SNPs to genes through functional genomics
- Prioritize candidate genes for further validation

### Main modules:
1. **eQTL lookup** – Match SNPs to expression QTLs
2. **Colocalization** – Test for shared causal signals
3. **TWAS** – Predict gene-trait associations
4. **Prioritization** – Rank genes based on multi-omic evidence

---

