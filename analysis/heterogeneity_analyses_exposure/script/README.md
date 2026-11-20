# Cross-population heterogeneity of exposure SNP effects

`heterogeneity_exposures.Rmd` tests whether the effects of exposure-associated SNPs differ between populations. It runs the analysis for every exposure, then again for the subset of SNPs in drug-target genes. Rendering it produces `heterogeneity_exposures.html`, which has one section per trait.

## Traits analysed

The report picks up every file in `../data/` whose name ends in `_pop_extract_common_snps_clumped_corrected.rds`. The trait name is the part of the file name before that suffix. The data currently has six traits:

| Trait | Input file prefix | SNPs per population | Drug classes (target genes) in the drug-target subset |
| --- | --- | --- | --- |
| Systolic blood pressure | `sbp` | 154 | CCB (*CACNA1C*, *CACNA1D*, *CACNB2*, *CACNB3*), BB (*ADRB1*), ACEI (*ACE*) |
| Diastolic blood pressure | `dbp` | 141 | CCB (*CACNA1C*, *CACNA1D*, *CACNA2D2*, *CACNB2*, *CACNB3*), BB (*ADRB1*), ACEI (*ACE*) |
| Pulse pressure | `pp` | 188 | CCB (*CACNA1C*, *CACNA1D*, *CACNB2*, *CACNB3*, *CACNG1*), BB (*ADRB1*) |
| LDL cholesterol | `ldlc` | 297 | Statins (*HMGCR*, *PCSK9*), Ezetimibe (*NPC1L1*) |
| HDL cholesterol | `hdl` | 141 | Statins (*HMGCR*, *PCSK9*), Ezetimibe (*NPC1L1*) |
| Triglycerides | `tg` | 121 | Statins (*HMGCR*, *PCSK9*), Ezetimibe (*NPC1L1*) |

CCB = calcium channel blockers, BB = beta blockers, ACEI = ACE inhibitors. The drug class and target gene come from each file's `Drug` and `Target gene` columns.

Each SNP has an estimate from six populations: AFR, AMR, EAS, EUR, SAS and UGR (Ugandan Genome Resource). See the main repository README for data sources. The lipid files also contain `CROSS_ANCESTRY` (meta-analysed) rows. These are left out so that each population counts once. To change which populations are used, edit the `populations` vector in the setup chunk.

## Analyses

The report runs the following steps for each trait:

1. **Cochran's Q heterogeneity test.** For each SNP, the population estimates are combined with an inverse-variance-weighted (IVW) fixed-effect meta-analysis. Cochran's Q then tests whether the population effects differ more than chance would allow (df = number of populations − 1).
2. **Multiple-testing correction.** The Bonferroni threshold is 0.05 divided by the number of SNPs tested. Benjamini–Hochberg FDR-adjusted p-values are reported alongside.
3. **Per-population Q contributions.** Each population's contribution to Q, `w_j(β_fixed − β_j)²`, is tested on 1 df. This shows which population drives the heterogeneity at a SNP. The results are shown as tables, a contribution plot, and QQ plots with genomic inflation (λ) for each population.
4. **Forest plots of Bonferroni-significant SNPs.** These show each population's estimate with its 95% CI.
5. **Expected versus observed replication.** EUR is the discovery population, using SNPs with EUR p < 5×10⁻⁸. Based on the EUR effect and the other population's standard error, the report works out how many SNPs should have the same direction of effect and how many should reach p < 0.05 in that population. It compares this with what is actually seen, using a binomial test. It also lists SNPs that were expected to replicate (probability > 0.95) but did not.
6. **Drug-target subset.** Steps 1–5 are repeated for SNPs that have a `Drug` annotation, with these differences:
   - The Bonferroni threshold is based on the number of drug-target SNPs.
   - Results are summarised by drug class.
   - The replication analysis uses EUR p < 0.05 to select discovery SNPs.
   - The forest plots of drug-target Bonferroni hits show:
     - one panel per SNP, labelled `chr:pos_<target gene>`, with the most significant SNP first
     - an extra "Combined" row with the IVW fixed-effect estimate (triangle)
     - a square for the population with the largest Q contribution
     - at most three SNPs per figure (set by `snps_per_page`).

The QQ plots and λ values are descriptive only. Most traits have between about 100 and 300 SNPs, so these figures are not very precise.

## Running the report

To analyse only some traits, set the YAML parameter at the top of the Rmd. Leave it as `traits: null` (the default) to run all traits.

```yaml
params:
  traits: ["sbp", "dbp"]
```

Render from the repository root:

```r
rmarkdown::render(
  "analysis/heterogeneity_analyses_exposure/script/heterogeneity_exposures.Rmd",
  output_file = "heterogeneity_exposures.html"
)
```

Each input file must have the columns `snpid`, `pop`, `beta` and `se`. The replication analysis also needs `pval`. The drug-target analysis needs `Drug`, and the gene labels on its forest plots come from `Target gene`. Required packages are `rmarkdown`, `knitr`, `dplyr`, `tidyr` and `ggplot2`.
