---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# Motif2Site

The goal of Motif2Site is to detect transcription factor binding sites using 
motifs IUPAC sequence or bed coordinates and ChIP-seq experiments in bed or bam
format. It also Combines/compares binding sites across experiments, tissues, or 
conditions.

## Installation

You can install the released version of Motif2Site from 
[github](https://github.com/ManchesterBioinference/Motif2Site) with:

``` 
devtools::install_github("ManchesterBioinference/Motif2Site")
```

## Example

This is a basic example which shows Motif2Site detects binding sites across 
tissue and combines/compares them:

```
library(Motif2Site)
library(BSgenome.Mmusculus.UCSC.mm10)

# HOX candidate motifs in Chr19 mouse
HOXMotifs = system.file("extdata", "HOXHomer_chr19.bed",
                        package="Motif2Site")

# HOXA2 BA2 peak calling from bed file
# ChIP-seq datasets in bam paired end format
IPs <- c(system.file("extdata", "HOXA2_BA2_rep1_IP_chr19.bam",
                     package="Motif2Site"),
         system.file("extdata", "HOXA2_BA2_rep2_IP_chr19.bam",
            package="Motif2Site")
         )
INPUTs <- c(system.file("extdata", "HOXA2_BA2_rep1_INPUT_chr19.bam",
                        package="Motif2Site"),
            system.file("extdata", "HOXA2_BA2_rep2_INPUT_chr19.bam",
            package="Motif2Site")
            )
HOXA2BA2_statistics <- DetectBindingSitesBed(BedFile= HOXMotifs,
                                             IPfiles = IPs,
                                             BackgroundFiles = INPUTs,
                                             genome="Mmusculus",
                                             genomeBuild="mm10",
                                             format="BAMPE",
                                             expName = "HOXA2_BA2"
                                             )

# HOXA3 PBA2 peak calling from bed file motif
# ChIP-seq datasets in bed single end format
IPs <- c(system.file("extdata", "HOXA3_PBA_rep1_IP_chr19.bed.bz2",
                     package="Motif2Site"),
         system.file("extdata", "HOXA3_PBA_rep2_IP_chr19.bed.bz2",
                     package="Motif2Site")
         )
INPUTs <- c(system.file("extdata", "HOXA3_PBA_rep1_INPUT_chr19.bed.bz2",
                        package="Motif2Site"),
            system.file("extdata", "HOXA3_PBA_rep1_INPUT_chr19.bed.bz2",
                        package="Motif2Site")
            )
HOXA3PBA_statistics <- DetectBindingSitesBed(BedFile= HOXMotifs,
                                             IPfiles = IPs,
                                             BackgroundFiles = INPUTs,
                                             genome="Mmusculus",
                                             genomeBuild="mm10",
                                             format="BEDSE",
                                             expName = "HOXA3_PBA"
                                             )

# HOXA2 PBA peak calling from bed file motif
# ChIP-seq datasets in bam paired end format
IPs <- system.file("extdata","HOXA2_PBA_IP_chr19.bam", package="Motif2Site")
INPUTs <- system.file("extdata","HOXA2_PBA_INPUT_chr19.bam",
                      package="Motif2Site")
HOXA2PBA_statistics <- DetectBindingSitesBed(BedFile= HOXMotifs,
                                             IPfiles = IPs,
                                             BackgroundFiles = INPUTs,
                                             genome="Mmusculus",
                                             genomeBuild="mm10",
                                             format="BAMPE",
                                             expName = "HOXA2_PBA",
                                              fdrValue = 0.1
                                              )

# Combine all HOX binding sites into one table
corMAT <- recenterBindingSitesAcrossExperiments(
  expLocations=c("HOXA2_BA2","HOXA2_PBA", "HOXA3_PBA"),
  experimentNames=c("HOXA2BA2","HOXA2PBA", "HOXA3PBA"),
  expName="combinedHOX",
  fdrValue=0.05
  )
corMAT
matFile <- paste0(getwd(), "/combinedHOX/CombinedMatrix")
HOXTable <- read.table(matFile, header = TRUE, check.names = FALSE)
head(HOXTable)


```