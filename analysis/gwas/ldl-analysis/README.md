# Perform cross-ancestry LDL cholesterol analysis using GLGC datasets 

## Setup

Make sure you have `bcftools` installed and available on the path e.g. on bc4

```
module add apps/bcftools-1.9-74/1.9-74
module add apps/plink/1.90
```

## Clump cross-ancestry GWAS

```
sh glgc_cross_ancestry_ldl_ld_clumping.sh
```

## Organise data and extract relevant variants

```
Rscript organise_data.r
```

Note paths are hardcoded - need to update.

## Perform analysis
Using LDLc clumped SNPs from the meta-analysis
```
quarto render ldl_analyse.qmd
```
Using LDLc clumped SNPs specific to LDLc drug targets 
```
quarto render ldl_analyses.drugtarget.qmd
```

## Follow-up analysis examining heterogenous SNPs across ancestries
GLGC-Examination-of-heteogenous-SNPs-across-populations.html 
```
