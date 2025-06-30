# NMIBC Gene Prioritization Pipeline

Een Nextflow-gestuurde workflow voor de integratie van GWAS-, eQTL- en TWAS-analyses om genen te prioriteren in niet–spierinvasief blaaskanker (NMIBC). Ontworpen voor uitvoering op een HPC-cluster (bijv. Snellius) met Singularity-containers.

---

## Inhoudsopgave

1. [Projectstructuur](#projectstructuur)
2. [Functies](#functies)
3. [Vereisten](#vereisten)
4. [Installatie](#installatie)
5. [Configuratie](#configuratie)
6. [Workflowcomponenten](#workflowcomponenten)

   * [scripts/](#scripts)
   * [modules/](#modules)
   * [main.nf](#mainnf)
   * [nextflow.config](#nextflowconfig)
7. [Uitvoering](#uitvoering)
8. [Outputbestanden](#outputbestanden)
9. [Logbestanden](#logbestanden)
10. [Licentie](#licentie)

---

## Projectstructuur

```
├── conf/                # Extra Nextflow- of cluster-instellingen
├── containers/          # Singularity-images
├── Data/                # Inputdata (custom)
├── GTEx/                # Gecontroleerde eQTL-data
├── scripts/             # R-scripts voor analyses
├── modules/             # Nextflow-modules per type analyse
├── r_Scripts/           # Hulpscripts in R
├── TCGA/                # TCGA-expressiemodellen en pos-bestanden
├── TWAS/                # Overige TWAS-data
├── main.nf              # Hoofd-workflow
├── nextflow.config      # Parameter- en containerconfiguratie
├── README.md            # Dit bestand
├── results/             # Gekopieerde output van de workflow
└── work/                # Nextflow-werkdirectory (intermediair)
```

---

## Functies

* **GWAS-filtering** op p-waarde (aanpasbare drempel)
* **eQTL-lookup**: harmonisatie en matching tussen GWAS- en eQTL-varianten
* **COLOC-analyse** voor colocalisatie tussen eQTL en GWAS
* **FUSION TWAS** voor transcriptome-wide associatiestudies
* **Reproduceerbaarheid** via Singularity-containers

---

## Vereisten

* [Nextflow](https://www.nextflow.io) (v20+)
* Java (OpenJDK 11 of hoger)
* Singularity (geïnstalleerd op HPC)
* R (3.6+) met packages: `coloc`, `data.table`, `dplyr`, `R.utils`

---

## Installatie

1. **Clone de repository**

   ```bash
   git clone https://github.com/<gebruikersnaam>/nmibc-geneprio.git
   cd nmibc-geneprio
   ```
2. **Zorg dat Nextflow en Singularity beschikbaar zijn**
3. **Pas `nextflow.config` aan** (data- en containerpaden)

---

## Configuratie

Alle inputpaden, parameters en container-instellingen staan in `nextflow.config`:

```groovy
params {
  gwas_file        = '<pad>/MA_STDERR_gwas_prognose.tbl'
  eqtl_lookup_file = '<pad>/GTEx_Analysis_v10_eQTL/Bladder.v10.eGenes.txt.gz'
  eqtl_coloc_file  = '<pad>/filtered_merged_data_eqtl_gtex_results'
  weights_dir      = '<pad>/TCGA-BLCA.TUMOR'
  ld_dir           = '<pad>/fusion_LD'
  pos_file         = '<pad>/TCGA-BLCA.TUMOR.pos'
  pvalue_threshold = 0.1
}

singularity.enabled = true
singularity.autoMounts = true

process {
  container = '<pad>/containers/fusion.sif'
  containerOptions = '--bind <pad>/containers/fusion_twas:/fusion_twas'
}

profiles {
  singularity {
    process.withName: 'COLOC_ANALYSIS' {
      container = '<pad>/containers/fusion.sif'
      containerOptions = '--bind <pad>/containers/fusion_twas:/fusion_twas'
    }
  }
}
```

> Pas `<pad>` aan naar je lokale paden op de HPC.

---

## Workflowcomponenten

### scripts/

Alle R-scripts die in de modules worden aangeroepen:

* **coloc.R**: voert `coloc.abf` uit en genereert `coloc_result.txt` + credible set (`coloc_cs_snps.txt`).
* **eqtl\_lookup.R**: harmoniseert en matcht GWAS- en eQTL-varianten, schrijft `eqtl_lookup_result.txt`.
* **fusion\_twas.R**: runt FUSION TWAS per gen, maakt `*.fusion.out` en een samengevoegd bestand.

### modules/

Nextflow-processen per analyse:

* **coloc.nf** → `COLOC_ANALYSIS`
* **eqtl\_lookup.nf** → `EQTL_LOOKUP`
* **fusion\_twas.nf** → `FUSION_TWAS`

Elk proces leest input (`script`, inputbestanden), voert Rscript uit en kopieert output naar `results/`.

### main.nf

Orkestreert de volledige workflow:

1. **GWAS-filtering** (`FILTER_GWAS`): filtert op p-waarde, telt SNPs.
2. **EQTL\_LOOKUP**: harmoniseert gefilterde GWAS met eQTL.
3. **COLOC\_ANALYSIS**: colocalisatie-analyse met GTEx eQTL.
4. **FUSION\_TWAS**: transcriptome-wide associatiestudie met FUSION.

### nextflow\.config

Zie [Configuration](#configuratie).

---

## Uitvoering

```bash
nextflow run main.nf -profile singularity
```

---

## Outputbestanden

Alle workflows schrijven naar `results/`:

* `results/gwas/filtered_gwas.tbl`, `snp_count.log`
* `results/eqtl_lookup/eqtl_lookup_result.txt`
* `results/coloc/coloc_result.txt`, `coloc_cs_snps.txt`
* `results/fusion/` (per-gen `.fusion.out` en `merged_fusion_results.txt`)

---

## Logbestanden

* Interne logs in `work/`
* Console-output in `logs/` (Slurm stdout)

## Licentie

Dit project valt onder de MIT License. Zie [LICENSE](LICENSE) voor details.
