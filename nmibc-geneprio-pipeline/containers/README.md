# Singularity Container

This container provides a reproducible, portable R environment for performing runnig the nmibc-geneprio pipeline. It includes preinstalled system libraries, R packages, and key scripts to support bioinformatics workflows in high-performance computing (HPC) environments.

Base Image
rocker/r-base:4.3.1 — a minimal R environment based on Debian.

Features
Preinstalled R packages for GWAS and eQTL workflows:

Rcpp, RcppEigen

coloc, susieR

ggplot2, viridis for plotting

optparse, R.utils, data.table, dplyr for data manipulation and scripting

Precompiled plink2R for handling PLINK genotype files in R

All system dependencies required for R package compilation (e.g., libxml2-dev, libssl-dev, gfortran, etc.)

Environment Variables
R_LIBS_USER: Points to the default R library path in the container (/usr/local/lib/R/site-library)

PATH: Includes /usr/local/bin (standard Rscript location)

LC_ALL=C: Prevents locale issues

Building the Container
To build this container with Singularity:

bash
sudo singularity build my_container.sif Singularity
