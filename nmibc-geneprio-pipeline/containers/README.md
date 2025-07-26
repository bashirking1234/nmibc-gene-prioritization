# Singularity Container for NMIBC-GenePrio Pipeline

This container provides a reproducible, portable R environment for running the **nmibc-geneprio** pipeline. It includes preinstalled system libraries, essential R packages, and specialized tools to support bioinformatics workflows in high-performance computing (HPC) environments.

---

## Base Image

- **[`rocker/r-base:4.3.1`](https://hub.docker.com/r/rocker/r-base)** — A minimal R environment based on Debian.

---

## Features

### R Packages for GWAS, eQTL, and Visualization

- `Rcpp`, `RcppEigen`
- `coloc`, `susieR`
- `ggplot2`, `viridis` (for plotting)
- `optparse`, `R.utils`, `data.table`, `dplyr` (for data manipulation and scripting)

### Additional Tools

- **Precompiled [`plink2R`](https://github.com/gabraham/plink2R)** for working with PLINK genotype files directly in R.
- **System dependencies** needed to compile R packages:
  - `libxml2-dev`, `libssl-dev`, `gfortran`, `libcurl4-openssl-dev`, etc.

---

## Environment Variables

| Variable         | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| `R_LIBS_USER`     | Points to the default R library path in the container (`/usr/local/lib/R/site-library`) |
| `PATH`            | Includes `/usr/local/bin` for consistent `Rscript` execution               |
| `LC_ALL`          | Set to `C` to avoid locale-related issues                                  |

---

## 🛠️ Building the Container

To build this container locally with Singularity:

```bash
sudo singularity build my_container.sif Singularity


License
This container is based on rocker/r-base and includes open-source R packages. Please refer to each package’s license for more details.

