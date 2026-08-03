# benchmarking-supplementary-scripts

This repository contains the scripts used for the analyses presented in **"Metabarcoding of high-diversity microbial communities exposes fundamental trade-offs in clustering and denoising of 16S amplicon data under controlled benchmarking conditions"** by Stamsaas et al. 2026

---

## Repository Structure

| File  | Description  |
|-----------|-----------|
| `dada2_pipeline.Rmd`                 | R Markdown notebook containing the code required to run the DADA2 pipeline. | 
| `prepare_sim_data.Rmd`               | R Markdown notebook for preparing input data for sequence simulation. | 
| `readMachine_amplicon_illumina.sh`   | Shell script for simulating Illumina amplicon sequencing reads. |
| `real_environment_data_analysis.Rmd` | R Markdown notebook for analysing real environmental sequencing data. |
| `rsearch_pipelines.Rmd`              | R Markdown notebook containing the code required to run the Rsearch pipelines, including cluster_size and UNOISE. |
| `swarm.sh`                           | Shell script for running Swarm clustering. |
| `swarm_postprocessing.Rmd`           | R Markdown notebook for post-processing and organising the Swarm clustering results. |

---

## Dependencies

| Name      | Type      | Version   | Link      |
|-----------|-----------|-----------|-----------|
| R          | Software  | 4.5.1     | https://docs.posit.co/ide/user/ |
| Rsearch    | R package | 1.1.0     | https://CRAN.R-project.org/package=Rsearch |
| tidyverse  | R package | 2.0.0     | https://CRAN.R-project.org/package=tidyverse |
| microseq   | R package | 2.1.7     | https://CRAN.R-project.org/package=microseq |
| dada2      | R package | 1.36      | https://bioconductor.posit.co/packages/3.23/bioc/html/dada2.html |
| data.table | R package | 1.18.4    | https://CRAN.R-project.org/package=data.table |
| Swarm      | Software  | 3.1.6     | https://github.com/torognes/swarm |
| VSEARCH    | Software  | 2.30.0    | https://github.com/torognes/vsearch |

---

## Data Availability

The input data required to reproduce the analyses can be obtained as described in the article.

---

## Citation

When using these scripts, please cite:

Stamsaas, C., et al. (2026). *Metabarcoding of high-diversity microbial communities exposes fundamental trade-offs in clustering and denoising of 16S amplicon data under controlled benchmarking conditions*. [Journal]. [DOI]

---

## Contact

Cassandra Stamsaas<br> 
Norwegian University of Life Sciences<br>
Email: cassandra.stamsaas@nmbu.no
<br>
<br>
Hilde Vinje<br> 
Norwegian University of Life Sciences<br>
Email: hilde.vinje@nmbu.no
