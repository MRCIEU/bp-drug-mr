# Perform cross-ancestry LDL cholesterol analysis

## Setup

Make sure you have `bcftools` installed and available on the path e.g. on bc4

```
module add apps/bcftools-1.9-74/1.9-74
module add apps/plink/1.90
```

## Organise data and extract relevant variants

```
Rscript organise_data.r
```

Note paths are hardcoded - need to update.

## Perform analysis

```
quarto render analyse.qmd
```
