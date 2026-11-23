
# Two-sample Mendelian randomization by population

`mr.R` runs two-sample Mendelian randomization (MR) of blood pressure and LDL cholesterol on cardiovascular outcomes, separately in each population, using the [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/) package. It then pools the population-specific estimates in a meta-analysis. At the end it renders `mr_report.Rmd` into an HTML report.

## Files

| File | Purpose |
| --- | --- |
| `script/mr.R` | Main analysis script: formats, harmonises, runs MR and sensitivity analyses, meta-analyses, writes CSVs, PDFs and the HTML report |
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

`CROSS_ANCESTRY` rows are never used.

### Steps

1. **Format.** For each pair and each population, the exposure and outcome rows for that population are reformatted with `format_data()`:
   - `snp_col = "rsid"`, `ea`/`oa` alleles, `eaf`, `beta`, `se`, `pval`;
   - the outcome is restricted to the exposure SNPs with `snps = exposure_dat$SNP`.
2. **Harmonise.** `harmonise_data()` runs with its default settings (`action = 2`):
   - alleles are aligned to the exposure effect allele;
   - palindromic SNPs are aligned using allele frequencies;
   - palindromic SNPs with intermediate allele frequencies are flagged `mr_keep = FALSE` and excluded from MR.

   The exposure instruments are already LD clumped, so there is no clumping step.
3. **MR.** `mr()` runs inverse-variance weighted (IVW), MR-Egger, weighted median, simple mode and weighted mode. The Wald ratio is used if only one SNP is available. The outcomes are binary, so `generate_odds_ratios()` adds odds ratios per unit increase in the exposure.
4. **Sensitivity analyses.**
   - `mr_heterogeneity()`: Cochran's Q across SNPs.
   - `mr_pleiotropy_test()`: MR-Egger intercept.
   - `mr_singlesnp()`: single-SNP estimates.
   - `mr_leaveoneout()`: leave-one-out estimates.
5. **Meta-analysis across populations.** For each pair, the six population IVW estimates are pooled by inverse-variance weighting. Each population's Wald ratio is used instead if it had only one SNP.
   - Fixed effect: β = Σwⱼβⱼ / Σwⱼ, with wⱼ = 1/SEⱼ².
   - Random effects: DerSimonian–Laird.
   - Between-population heterogeneity: Cochran's Q (df = k − 1) and I².

## Running the script

Set your R module and run from any directory. Paths are worked out from the script's location.

```bash
Rscript analysis/mendelian_randomization/script/mr.R                              # all pairs
Rscript analysis/mendelian_randomization/script/mr.R exposures=sbp outcomes=chd   # one pair
Rscript analysis/mendelian_randomization/script/mr.R exposures=sbp,dbp populations=AFR,UGR
```

Command-line options are `key=value`, with comma-separated values: `exposures`, `outcomes`, `populations`. Progress messages with timestamps go to the job log. A full run takes about 5 minutes, most of it for the two mode-based methods.

Other settings are in the CONFIGURATION block at the top of `mr.R`:

| Setting | Default | Meaning |
| --- | --- | --- |
| `populations` | AFR, AMR, EAS, EUR, SAS, UGR | Populations to analyse |
| `outcome_population_map` | `c(UGR = "AFR")` | Outcome population to use when it differs from the exposure population |
| `mr_methods` | 6 methods above | Methods passed to `mr()` |
| `make_plots` | `TRUE` | Write PDF plots |
| `make_html_report` | `TRUE` | Render `mr_report.Rmd` to `results/mr_report.html` |

Required packages are `TwoSampleMR`, `dplyr`, `purrr` and `ggplot2`. The HTML report also needs `rmarkdown`, `knitr` and pandoc. If pandoc is not available on the HPC, set `make_html_report <- FALSE`.

> **Running pairs as separate jobs:** if each job runs a subset of the pairs, the combined files below are overwritten by each job. The per-pair files are not affected. Combined files: `mr_all_results.rds`, `meta_ivw_all.csv`, `harmonisation_summary.csv`, the PDFs and `mr_report.html`. Run all pairs in one job when you need these.

## Outputs (`results/`)

Files named `*_<exposure>_<outcome>.*` hold one exposure–outcome pair, with one row per population (and per method, where relevant).

| File | Contents |
| --- | --- |
| `harmonised_<exp>_<out>.rds` | Harmonised data, with `population` and `outcome_population` columns |
| `mr_<exp>_<out>.csv` | `mr()` estimates (b, se, pval, OR and 95% CI) by population and method |
| `heterogeneity_<exp>_<out>.csv` | Cochran's Q across SNPs (IVW and MR-Egger) |
| `pleiotropy_<exp>_<out>.csv` | MR-Egger intercept test |
| `singlesnp_<exp>_<out>.csv` | Single-SNP MR estimates |
| `leaveoneout_<exp>_<out>.csv` | Leave-one-out MR estimates |
| `meta_ivw_<exp>_<out>.csv` | Pooled fixed- and random-effects estimates, Q, τ², I² |
| `meta_ivw_all.csv` | Pooled estimates for all pairs |
| `harmonisation_summary.csv` | SNPs harmonised, kept for MR, palindromic and dropped, by pair and population |
| `mr_all_results.rds` | List of all of the above (`harmonised`, `mr`, `heterogeneity`, `pleiotropy`, `singlesnp`, `leaveoneout`, `meta_ivw`) |
| `plots/meta_ivw_forest.pdf` | IVW estimate for each population plus pooled estimates, all pairs |
| `plots/mr_methods_forest.pdf` | All MR methods by population, one page per pair |
| `plots/mr_scatter.pdf` | SNP-outcome against SNP-exposure scatter plots with MR fits, one page per pair and population |
| `mr_report.html` | HTML report (see below) |

## HTML report (`mr_report.html`)

The report reads `mr_all_results.rds` and contains:

1. **Overview.** The methods in brief.
2. **Harmonisation.** For each pair and population: SNPs harmonised, kept for MR, palindromic, and ambiguous palindromic SNPs dropped.
3. **Meta-analysis of IVW estimates across populations.** A table of pooled fixed- and random-effects odds ratios with Q, its p-value and I². A forest plot shows each population's IVW estimate (grey) with the two pooled estimates (red diamonds), with exposures as rows and outcomes as columns.
4. **Results by exposure–outcome pair.** One tab per pair, containing:
   - estimates from every MR method for each population;
   - Cochran's Q across SNPs;
   - the MR-Egger intercept test;
   - a forest plot of the methods for each population.
5. **Scatter plots.** One plot per pair and population.

To re-render the report from existing results without re-running the analysis:

```r
rmarkdown::render("analysis/mendelian_randomization/script/mr_report.Rmd",
                  params = list(results = "analysis/mendelian_randomization/results/mr_all_results.rds"),
                  output_dir = "analysis/mendelian_randomization/results")
```

## Interpretation notes

- **UGR and AFR are not independent.** Both use the AFR outcome sample. Including both in the meta-analysis means the pooled standard errors are somewhat too small.
- **Fixed versus random effects.** Heterogeneity between populations is high for most pairs (I² 68–97% across the 12 pairs), so the random-effects pooled estimate is usually the more appropriate one to report.
- **Mode-based estimates.** In populations with small samples (e.g. SAS, UGR), the simple and weighted mode estimates can have very large standard errors and uninformative confidence intervals. Do not read these as evidence of no effect.
- **Effect scale.** Odds ratios are per unit of the exposure as provided in the exposure files.
