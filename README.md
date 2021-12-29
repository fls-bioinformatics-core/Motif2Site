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
library(BSgenome.Ecoli.NCBI.20080805)

# FUR candidate motifs in NC_000913 E. coli
FurMotifs = system.file("extdata", "FurMotifs.bed", package="Motif2Site")

# ChIP-seq FUR fe datasets binding sites from user provided bed file 
# ChIP-seq datasets in bed single end format

IPFe <- c(system.file("extdata", "FUR_fe1.bed", package="Motif2Site"),
          system.file("extdata", "FUR_fe2.bed", package="Motif2Site"))
Inputs <- c(system.file("extdata", "Input1.bed", package="Motif2Site"),
            system.file("extdata", "Input2.bed", package="Motif2Site"))
FURfeBedInputStats <- 
  DetectBindingSitesBed(BedFile=FurMotifs,
                        IPfiles=IPFe, 
                        BackgroundFiles=Inputs, 
                        genome="Ecoli",
                        genomeBuild="20080805",
                        DB="NCBI",
                        expName="FUR_Fe_BedInput",
                        format="BEDSE"
                        )

# ChIP-seq FUR dpd datasets binding sites from user provided bed file 
# ChIP-seq datasets in bed single end format

IPDpd <- c(system.file("extdata", "FUR_dpd1.bed", package="Motif2Site"),
           system.file("extdata", "FUR_dpd2.bed", package="Motif2Site"))
FURdpdBedInputStats <- 
  DetectBindingSitesBed(BedFile=FurMotifs,
                        IPfiles=IPDpd, 
                        BackgroundFiles=Inputs, 
                        genome="Ecoli",
                        genomeBuild="20080805",
                        DB="NCBI",
                        expName="FUR_Dpd_BedInput",
                        format="BEDSE"
                        )

# Combine FUR binding sites from bed input into one table 
corMAT <- recenterBindingSitesAcrossExperiments(
  expLocations=c("FUR_Fe_BedInput","FUR_Dpd_BedInput"),
  experimentNames=c("FUR_Fe","FUR_Dpd"),
  expName="combinedFUR"
  )
corMAT

FurTable <- 
  read.table(file.path("combinedFUR","CombinedMatrix"), 
             header = TRUE,
             check.names = FALSE
             )
FurBindingTotal <- 
  GRanges(seqnames=Rle(FurTable[,1]), 
          ranges = IRanges(FurTable[,2], FurTable[,3])
          )
FurFe <- FurBindingTotal[which((FurTable$FUR_Fe_binding =="Binding")==TRUE)]
FurDpd <- FurBindingTotal[which((FurTable$FUR_Dpd_binding =="Binding")==TRUE)]
findOverlaps(FurFe,FurDpd) 


# Differential binding sites across FUR conditions fe vs dpd
diffFUR <- pairwisDifferential(tableOfCountsDir="combinedFUR",
                              exp1="FUR_Fe",
                              exp2="FUR_Dpd",
                              FDRcutoff = 0.05,
                              logFCcuttoff = 1
                              )
FeUp <- diffFUR[[1]]
DpdUp <- diffFUR[[2]]
TotalComparison <- diffFUR[[3]]
head(TotalComparison)


```