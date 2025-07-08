# Fusion TWAS Singularity Container

This container provides a reproducible, portable R environment for performing Transcriptome-Wide Association Studies (TWAS) using FUSION. It includes preinstalled system libraries, R packages, and key scripts to support bioinformatics workflows in high-performance computing (HPC) environments.

---

##  Build Instructions

To build the container:

```bash
singularity build fusion.sif singularity.def
```

To save a fresh fixed version (e.g., after modifying the definition):

```bash
singularity build fusion_fixed.sif singularity.def
```

---

## Contents of the Container

### Copied Files

| Source Path                                         | Destination in Image    |
| --------------------------------------------------- | ----------------------- |
| `/gpfs/home2/hbashir1/tools/plink2R/plink2R`        | `/opt/plink2R`          |
| `fusion_twas/FUSION.assoc_test.R`                   | `/opt/fusion/`          |
| `fusion_twas/FUSION.compute_weights.R`              | `/opt/fusion/`          |
| `fusion_twas/FUSION.post_process.R`                 | `/opt/fusion/`          |
| `fusion_twas/gcta_nr_robust`                        | `/opt/fusion/`          |
| `fusion_twas/utils/` *(directory)*                  | `/opt/fusion/utils/`    |

### Symlinks Created

Inside the container’s `%post` section:

```bash
ln -s /opt/fusion/FUSION.assoc_test.R      /usr/local/bin/FUSION.assoc_test.R
ln -s /opt/fusion/FUSION.compute_weights.R /usr/local/bin/FUSION.compute_weights.R
ln -s /opt/fusion/FUSION.post_process.R    /usr/local/bin/FUSION.post_process.R
ln -s /opt/fusion/gcta_nr_robust           /usr/local/bin/gcta64
```

---

## Environment Variables

Set automatically at container runtime:

```bash
export R_LIBS_USER=/usr/local/lib/R/site-library
export PATH=/opt/fusion:/usr/local/bin:$PATH
export LC_ALL=C
```

---

## System & R Package Installation

###  System Dependencies

Installed via `apt`:

- `wget`, `unzip`, `curl`, `git`
- `libcurl4-openssl-dev`, `libxml2-dev`, `libssl-dev`
- `libblas-dev`, `liblapack-dev`, `libpng-dev`
- `libgfortran5`, `libxt6`, `fonts-dejavu`
- `libharfbuzz-dev`, `libfribidi-dev`, `libfreetype6-dev`
- `libtiff5-dev`, `libjpeg-dev`, `libcairo2-dev`, `libfontconfig1-dev`
- `build-essential`, `gfortran`

### R Packages Installed via CRAN

In your build script or `install_r_packages.R`:

```r
install.packages(c(
  'remotes', 'R.utils', 'coloc', 'susieR', 'ggplot2',
  'viridis', 'optparse', 'data.table', 'Rcpp',
  'RcppEigen', 'dplyr'
), repos = 'https://cloud.r-project.org')
```

---

## `plink2R` Installation Attempt

The build script tries to install `plink2R` from a local directory:

```bash
R -e "install.packages('/opt/plink2R', repos = NULL, type = 'source')" \
  || echo 'FOUT: plink2R lokaal faalde'
```

If this package is missing or improperly structured, the build continues without halting due to the `|| echo` fallback.

### Alternative Fix

For a more reliable approach, replace the above line in your definition with:

```bash
R -e "remotes::install_github('gabraham/plink2R')"
```

Then rebuild the container:

```bash
singularity build fusion_fixed.sif singularity.def
```

---

## Usage Examples

- **Run a FUSION R script inside the container:**

  ```bash
  singularity exec fusion.sif Rscript /opt/fusion/FUSION.assoc_test.R
  ```

- **Open an interactive R shell:**

  ```bash
  singularity shell fusion.sif
  ```

- **Test if `plink2R` was successfully installed:**

  ```bash
  singularity exec fusion.sif Rscript -e "library(plink2R)"
  ```

---

## Project Structure

```
.
├── fusion.sif
├── fusion_fixed.sif          # (optional future build)
├── singularity.def
├── build.log
├── install_r_packages.R
├── fusion_twas/
│   ├── FUSION.assoc_test.R
│   ├── FUSION.compute_weights.R
│   ├── FUSION.post_process.R
│   ├── gcta_nr_robust
│   └── utils/
└── tools/
    └── plink2R/
```

---

## License Notes

- R and its packages are governed by the **GNU General Public License (GPL)**.
- Docker base image from the **Rocker Project**.
- FUSION scripts are subject to their original licensing terms.
