[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)


### Setting up required JSON file
```
nf-core schema build -d <pipeline_directory> --no-prompts
```

### Running from command line
```bash
nextflow run main.nf -c nextflow.config -work-dir /your/work/dir -profile docker
```

### Running as a jobscript
Make sure to have the cores specified in ``-pe smp`` the same as in ``--cpus``.
In case nextflow is in your anaconda/miniconda environment, it is required to use ``eval "$(conda shell.bash hook)"`` prior to ``conda activate`` and ``nextflow run``.
```bash
#!/bin/bash
#$ -N jobname
#$ -o jobname.out
#$ -e jobname.err
#$ -V
#$ -pe smp 5

eval "$(conda shell.bash hook)" 
conda activate nextflow
cd $HOME
nextflow run main.nf -c nextflow.config -work-dir /your/work/dir -profile docker

```

# main.nf:Het hoofdscript waar je de volgorde van de drie modules definieert (eQTL → coloc → TWAS).
# nextflow.config: Voor algemene instellingen zoals output directory, aantal cores per proces, werkdirectory, enz.
# modules/ : Hier zet je je eigen process-modules, zoals eqtl_lookup.nf, coloc.nf, twas.nf.
# conf/: Handig voor profielen (bijv. lokaal vs. Snellius) en voor resources per proces.
# README.md 





