# Two-sample Mendelian randomization by population

`mr.R` runs two-sample Mendelian randomization (MR) of blood pressure and LDL cholesterol on cardiovascular outcomes, separately in each population, using the [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/) package. It then pools the population-specific estimates into a **cross-ancestry** estimate by meta-analysis. The analyses are run once with all instruments and again for each drug class (e.g. beta-blockers, calcium channel blockers, ACE inhibitors, statins, ezetimibe). At the end it renders `mr_report.Rmd` into an HTML report.

## Files

| File | Purpose |
| --- | --- |
| `script/mr.R` | Main analysis script: formats, harmonises, runs MR and sensitivity analyses (all instruments and by drug class), cross-ancestry meta-analysis, writes CSVs, PDFs and the HTML report |
| `script/mr_report.Rmd` | Report template. `mr.R` renders it from `results/mr_all_results.rds` |
| `data/pop_extract_<exposure>_<outcome>_clumped.rds` | Outcome SNP effects for each exposure's instruments, by population |
| `../heterogeneity_analyses_exposure/data/<exposure>_pop_extract_common_snps_clumped_corrected.rds` | Exposure SNP effects, by population |
| `results/` | All outputs (created when the script runs) |

## Analyses

The script finds every outcome file in `data/` and pairs it with the matching exposure file. The pairing uses the exposure name in the outcome file name. This gives 12 exposure–outcome pairs:

| Exposure | Outcomes |
| --- | --- |
| Systolic blood pressure (`sbp`) | Coronary heart disease (`chd`), hypertension (`hp`), stroke (`stroke`) |
| Diastolic blood pressure (`dbp`) | `chd`, `hp`, `stroke` |
| Pulse pressure (`pp`) | `chd`, `hp`, `stroke` |
| LDL cholesterol (`ldlc`) | `chd`, `hp`, `stroke` |

HDL cholesterol and triglycerides have exposure files but no outcome files, so they are not analysed.

### Populations

| Exposure population | Outcome population used |
| --- | --- |
| AFR, AMR, EAS, EUR, SAS | Same population |
| UGR (Ugandan Genome Resource) | AFR. There is no UGR outcome GWAS |

`CROSS_ANCESTRY` rows in the input files are never used. "Cross-ancestry" in the results always means the meta-analysis of the six populations above.

### Steps

1. **Format.** For each pair and each population, the exposure and outcome rows for that population are reformatted with `format_data()`:
   - `snp_col = "rsid"`, `ea`/`oa` alleles, `eaf`, `beta`, `se`, `pval`;
   - the outcome is restricted to the exposure SNPs with `snps = exposure_dat$SNP`.
2. **Harmonise.** `harmonise_data()` runs with its default settings (`action = 2`):
   - alleles are aligned to the exposure effect allele;
   - palindromic SNPs are aligned using allele frequencies;
   - palindromic SNPs with intermediate allele frequencies are flagged `mr_keep = FALSE` and excluded from MR.

   The exposure instruments are already LD clumped, so there is no clumping step.
3. **MR.** `mr()` runs inverse-variance weighted (IVW), MR-Egger, weighted median and weighted mode. The Wald ratio is used if only one SNP is available. The outcomes are binary, so `generate_odds_ratios()` adds odds ratios per unit increase in the exposure.
4. **Sensitivity analyses.**
   - `mr_heterogeneity()`: Cochran's Q across SNPs.
   - `mr_pleiotropy_test()`: MR-Egger intercept.
   - `mr_singlesnp()`: single-SNP estimates.
   - `mr_leaveoneout()`: leave-one-out estimates.
5. **Cross-ancestry meta-analysis.** For each pair, the six population IVW estimates are pooled by inverse-variance weighting. Each population's Wald ratio is used instead if it had only one SNP.
   - Fixed effect: β = Σwⱼβⱼ / Σwⱼ, with wⱼ = 1/SEⱼ².
   - Random effects: DerSimonian–Laird.
   - Between-population heterogeneity: Cochran's Q (df = k − 1) and I².

6. **Drug-class MR.** Steps 3–5 are repeated using only the instruments annotated with each drug class in the exposure file's `Drug` column (trailing spaces removed). The target genes come from the `Target gene` column.

   | Exposure | Drug classes (target genes) |
   | --- | --- |
   | SBP | CCB (*CACNA1C*, *CACNA1D*, *CACNB2*, *CACNB3*), BB (*ADRB1*), ACEI (*ACE*) |
   | DBP | CCB (*CACNA1C*, *CACNA1D*, *CACNA2D2*, *CACNB2*, *CACNB3*), BB (*ADRB1*), ACEI (*ACE*) |
   | PP | CCB (*CACNA1C*, *CACNA1D*, *CACNB2*, *CACNB3*, *CACNG1*), BB (*ADRB1*) |
   | LDL-C | Statins (*HMGCR*, *PCSK9*), Ezetimibe (*NPC1L1*) |

   CCB = calcium channel blockers, BB = beta-blockers, ACEI = ACE inhibitors. Several classes have only one or two SNPs, giving a Wald ratio or IVW estimate only. MR-Egger, weighted median and weighted mode need at least three SNPs. Not every drug SNP is in every outcome file. For example, SBP → CHD has no ACEI or BB SNPs, so those analyses are absent.

## Running the script

Set your R module and run from any directory. Paths are worked out from the script's location.

```bash
Rscript analysis/mendelian_randomization/script/mr.R                              # all pairs
Rscript analysis/mendelian_randomization/script/mr.R exposures=sbp outcomes=chd   # one pair
Rscript analysis/mendelian_randomization/script/mr.R exposures=sbp,dbp populations=AFR,UGR
```

Command-line options are `key=value`, with comma-separated values: `exposures`, `outcomes`, `populations`. Progress messages with timestamps go to the job log. A full run takes about 6 minutes. Most of the time goes on the weighted mode method.

Other settings are in the CONFIGURATION block at the top of `mr.R`:

| Setting | Default | Meaning |
| --- | --- | --- |
| `populations` | AFR, AMR, EAS, EUR, SAS, UGR | Populations to analyse |
| `outcome_population_map` | `c(UGR = "AFR")` | Outcome population to use when it differs from the exposure population |
| `mr_methods` | Wald ratio, IVW, MR-Egger, weighted median, weighted mode | Methods passed to `mr()` |
| `sig_threshold` | `0.05` | p-value threshold highlighted in the IVW forest plots |
| `run_drug_class` | `TRUE` | Run the drug-class MR |
| `drug_labels` | BB, CCB, ACEI, Statins, Ezetimibe | Full names for the drug-class codes |
| `make_plots` | `TRUE` | Write PDF plots |
| `make_html_report` | `TRUE` | Render `mr_report.Rmd` to `results/mr_report.html` |

Required packages are `TwoSampleMR`, `dplyr`, `purrr` and `ggplot2`. The HTML report also needs `rmarkdown`, `knitr` and pandoc. If pandoc is not available on the HPC, set `make_html_report <- FALSE`.

> **Running pairs as separate jobs:** if each job runs a subset of the pairs, the combined files below are overwritten by each job. The per-pair files are not affected. Combined files: `mr_all_results.rds`, `cross_ancestry_all.csv`, `harmonisation_summary.csv` (in both `results/` and `results/drug_class/`), the PDFs and `mr_report.html`. Run all pairs in one job when you need these.

## Outputs (`results/`)

Files named `*_<exposure>_<outcome>.*` hold one exposure–outcome pair, with one row per population (and per method, where relevant). The drug-class results are in `results/drug_class/`, with the same files named `*_<exposure>_<outcome>_<drug class>.*` (e.g. `mr_sbp_hp_CCB.csv`).

| File | Contents |
| --- | --- |
| `harmonised_<exp>_<out>.rds` | Harmonised data, with `population`, `outcome_population`, `drug_class` and `target_gene` columns |
| `mr_<exp>_<out>.csv` | `mr()` estimates (b, se, pval, OR and 95% CI) by population and method |
| `heterogeneity_<exp>_<out>.csv` | Cochran's Q across SNPs (IVW and MR-Egger) |
| `pleiotropy_<exp>_<out>.csv` | MR-Egger intercept test |
| `singlesnp_<exp>_<out>.csv` | Single-SNP MR estimates |
| `leaveoneout_<exp>_<out>.csv` | Leave-one-out MR estimates |
| `cross_ancestry_<exp>_<out>.csv` | Cross-ancestry fixed- and random-effects estimates, Q, τ², I² |
| `cross_ancestry_all.csv` | Cross-ancestry estimates for all pairs |
| `harmonisation_summary.csv` | SNPs harmonised, kept for MR, palindromic and dropped, by pair and population |
| `mr_all_results.rds` | List of all results (`harmonised`, `mr`, `heterogeneity`, `pleiotropy`, `singlesnp`, `leaveoneout`, `cross_ancestry`), with the drug-class results in `drug_class` |
| `plots/ivw_cross_ancestry_forest.pdf` | IVW estimate for each population plus cross-ancestry estimates, all pairs |
| `plots/mr_methods_forest.pdf` | All MR methods by population, one page per pair |
| `plots/mr_scatter.pdf` | SNP-outcome against SNP-exposure scatter plots with MR fits, one page per pair and population |
| `drug_class/plots/ivw_cross_ancestry_forest_drug_class.pdf` | As above by drug class, one page per exposure (drug classes as rows) |
| `drug_class/plots/mr_methods_forest.pdf`, `drug_class/plots/mr_scatter.pdf` | As above, for each drug class |
| `mr_report.html` | HTML report (see below) |

In the IVW forest plots:

- **red** marks estimates with p < `sig_threshold` (0.05 by default) and grey marks the rest;
- **circles** are population estimates and **diamonds** are cross-ancestry estimates;
- very wide 95% CIs (from single-SNP estimates in small samples) are truncated at 3× beyond the range of point estimates in that panel, and the plot caption says so.

## HTML report (`mr_report.html`)

The report reads `mr_all_results.rds` and contains:

1. **Overview.** The methods in brief and how to read the plots.
2. **Harmonisation.** For each pair and population: SNPs harmonised, kept for MR, palindromic, and ambiguous palindromic SNPs dropped.
3. **Cross-ancestry IVW estimates: all instruments.** A table of cross-ancestry fixed- and random-effects odds ratios with Q, its p-value and I². A forest plot shows each population's IVW estimate and the two cross-ancestry estimates, with exposures as rows and outcomes as columns.
4. **Results by exposure–outcome pair.** One tab per pair, containing:
   - estimates from every MR method for each population;
   - Cochran's Q across SNPs;
   - the MR-Egger intercept test;
   - a forest plot of the methods for each population.
5. **Drug-class MR.**
   - a table of the drug-class instruments (target genes and SNPs per population);
   - the cross-ancestry estimates by drug class;
   - one IVW forest plot per exposure, with drug classes as rows;
   - one tab per exposure–outcome–drug class with the same tables and plots as section 4.
6. **Scatter plots.** One plot per pair and population, for all instruments.

To re-render the report from existing results without re-running the analysis:

```r
rmarkdown::render("analysis/mendelian_randomization/script/mr_report.Rmd",
                  params = list(results = "analysis/mendelian_randomization/results/mr_all_results.rds"),
                  output_dir = "analysis/mendelian_randomization/results")
```

## Interpretation notes

- **UGR and AFR are not independent.** Both use the AFR outcome sample. Including both in the cross-ancestry meta-analysis means its standard errors are somewhat too small.
- **Fixed versus random effects.** Heterogeneity between populations is high for most pairs (I² 68–97% across the 12 pairs), so the random-effects cross-ancestry estimate is usually the more appropriate one to report.
- **Weighted mode estimates.** In populations with small samples (e.g. SAS, UGR), the weighted mode estimates can have very large standard errors and uninformative confidence intervals. Do not read these as evidence of no effect.
- **Effect scale.** Odds ratios are per unit of the exposure as provided in the exposure files.
- **Drug-class estimates.** These rest on few SNPs, often in one gene region, and many come from one SNP. Their CIs are wide in the smaller populations and they have little power for pleiotropy checks. The p < 0.05 highlighting is nominal and not corrected for the number of drug-class tests.
