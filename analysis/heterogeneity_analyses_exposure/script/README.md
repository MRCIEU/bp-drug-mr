# Cross-population heterogeneity of exposure SNP effects

`heterogeneity_exposures.Rmd` tests whether the effects of exposure-associated SNPs differ between populations, and how well European (EUR) SNP effects transfer to other populations. It runs the analysis for every exposure, then again for the subset of SNPs in drug-target genes. Rendering it produces `heterogeneity_exposures.html`, which has one section per trait, and saves cross-ancestry scatter figures to `../figures/`.

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

### About the input data

- **Alleles are harmonised to EUR.** Every non-EUR estimate uses the same effect allele as EUR for the same SNP. The cross-ancestry figures and slopes rely on this.
- **Non-EUR estimates are rescaled.** In the `_corrected` files, each non-EUR population's betas and standard errors have been divided by its regression slope against EUR (`beta / slope`, `se / slope`; see the rescaling script). This keeps each SNP's z-statistic, direction and p-value unchanged, but it puts all populations on the EUR scale. Slopes close to 1 in the figures are therefore partly forced by the rescaling, and are not evidence on their own that effect sizes agree.

## Analyses

The report runs the following steps for each trait:

1. **Cochran's Q heterogeneity test.** For each SNP, the population estimates are combined with an inverse-variance-weighted (IVW) fixed-effect meta-analysis. Cochran's Q then tests whether the population effects differ more than chance would allow (df = number of populations − 1).
2. **Multiple-testing correction.** The Bonferroni threshold is 0.05 divided by the number of SNPs tested. Benjamini–Hochberg FDR-adjusted p-values are reported alongside.
3. **Per-population contributions.** These show which population drives the heterogeneity at a SNP. Two methods are available, chosen with `contribution_method` in the setup chunk:
   - `"loo"` (default): a leave-one-out standardised residual. Each population's estimate is compared with the IVW estimate of all the other populations, `(β_j − β_−j) / √(se_j² + se_−j²)`. Under no heterogeneity this follows a standard normal distribution for every population.
   - `"qj"`: the original contribution `w_j(β_fixed − β_j)²` tested on 1 df. This is too small for populations that dominate the pooled estimate. In simulations with no heterogeneity, it gave λ ≈ 0.5 for EUR and ≈ 0.8 for AMR and EAS, compared with ≈ 1.0 for the leave-one-out method.

   The results are shown as tables, a contribution plot, and QQ plots with genomic inflation (λ) for each population. λ is the median of the chi-square statistics divided by its expected value under the null (0.455). For each population the report also gives:
   - **Standard error of λ.** By default this is a bootstrap: SNPs are resampled 1,000 times and the SE is the SD of the resampled λ values (`lambda_se_method <- "bootstrap"`). The alternative, `"asymptotic"`, uses the large-sample SE of a median, about 2.33 λ / √n. In simulations with 150 SNPs, both matched the true spread of λ and gave close to 95% coverage.
   - **95% CI and a test of λ = 1.** These are calculated on the log scale (SE of log λ = SE / λ), so the interval stays positive.
   - **Plots.** The QQ plot legend shows each λ with its 95% CI, and a grey band marks the pointwise 95% range expected under the null. A second plot shows each λ with its CI against a dashed line at λ = 1.
4. **Transferability versus EUR.** Each SNP is oriented to the EUR effect-increasing allele. For each non-EUR population, the report gives:
   - the IVW slope through the origin of non-EUR on EUR betas (weights 1/se²), with its 95% CI
   - the proportion of SNPs with the same direction of effect as EUR (sign concordance)
   - the Pearson correlation.
5. **Forest plots of Bonferroni-significant SNPs.** These show each population's estimate with its 95% CI.
6. **Expected versus observed replication.** EUR is the discovery population, using SNPs with EUR p < 5×10⁻⁸. For each SNP, the report works out the chance that it has the same direction of effect in the other population, and the chance that it reaches two-sided p < 0.05 there (the power, given the EUR effect and the other population's standard error). It adds these up to get the expected counts and compares them with the observed counts using a binomial test. It also lists SNPs that were expected to replicate (probability > 0.95) but did not.
7. **Drug-target subset.** Steps 1–3, 5 and 6 are repeated for SNPs that have a `Drug` annotation, with these differences:
   - The Bonferroni threshold is based on the number of drug-target SNPs.
   - Results are summarised by drug class.
   - The replication analysis uses EUR p < 0.05 to select discovery SNPs.
   - The forest plots of drug-target Bonferroni hits show:
     - one panel per SNP, labelled `chr:pos_<target gene>`, with the most significant SNP first
     - an extra "Combined" row with the IVW fixed-effect estimate (triangle)
     - a square for the population with the largest contribution to heterogeneity
     - at most three SNPs per figure (set by `snps_per_page`).
8. **Cross-ancestry figures.** After all traits are analysed, the report draws one grid for blood pressure (SBP, DBP, PP) and one for lipids (LDL-C, HDL-C, triglycerides). Rows are traits and columns are non-EUR populations. In each panel:
   - each point is one SNP, oriented to the EUR effect-increasing allele, so SNPs with the same direction as EUR sit above zero
   - the dashed line is y = x and the blue line is the IVW slope
   - grey bars are 95% CIs for the non-EUR estimate, cut off at the panel edge (turn off with `show_ci <- FALSE`)
   - the text gives the IVW slope (95% CI), sign concordance and number of SNPs
   - panels in the same row share axis limits.

   Which traits go in each figure is set by `figure_groups`, and their labels by `trait_labels`, both in the setup chunk.

The QQ plots and λ values are descriptive. With 100–300 SNPs per trait, the 95% CI for λ is roughly ±40% of its value. The drug-target subsets have about 10–140 SNPs, so their λ estimates are much less precise, and those with fewer than about 50 SNPs should not be interpreted.

## Output

### `heterogeneity_exposures.html`

The report has four parts:

1. **Analysis overview**: the input file used for each trait.
2. **Analysis functions**: the code for each step above. It is hidden by default; click *Code* to show it.
3. **Results by exposure**: one section per trait. Each section contains:
   - **Variant counts by population**
   - **Overall heterogeneity summary**: SNPs tested, Bonferroni threshold, and the number of SNPs with Q p < 0.05, Bonferroni-significant and FDR < 0.05
   - **Transferability versus EUR**: `n_snps`, `ivw_slope` with `ivw_slope_se`, `ivw_slope_lci` and `ivw_slope_uci`, `sign_concordance` and `pearson_r`
   - **SNPs passing the Bonferroni threshold**: the IVW estimate (`beta_fixed`, `se_fixed`), `Q`, `Qdf`, `Qpval` and `fdr`
   - **Per-population contributions**: p-values for every SNP with Q p < 0.05, a contribution plot, a QQ plot with λ and its 95% CI in the legend, a plot of λ with 95% CIs, and a table giving `n_snps`, `lambda`, `lambda_se`, `lambda_lci`, `lambda_uci` and `p_lambda_eq_1` for each population
   - **Forest plots** of the Bonferroni-significant SNPs
   - **Population-level contribution summary**: the proportion of SNPs where each population's contribution has p < 0.05
   - **Expected versus observed replication: all exposure SNPs**: expected and observed counts for sign and significance, with binomial `pdiff`, and the SNPs that were expected but failed to replicate (`expected_power`, `observed_sig`, `expected_sign`, `observed_sign`)
   - **Drug-target subset**: the same outputs for drug-target SNPs, summarised by drug class, followed by its own replication table.
4. **Cross-ancestry effect comparison figures**: the blood pressure and lipid grids, and a table of IVW slopes, sign concordance and correlations for every trait and population.

The report ends with the settings used and the R session information.

### `../figures/`

`cross_ancestry_blood_pressure.{pdf,png}` and `cross_ancestry_lipids.{pdf,png}` are the same grids as in the report, saved as files. To save them somewhere else, change `figure_dir`.

### Interpreting the results

- **Replication**: expected counts assume the true effect in each population equals the EUR estimate. EUR SNPs are chosen because they are significant, so their effects may be overestimated (winner's curse). Expected counts are therefore likely to be upper bounds.
- **Heterogeneity in AMR and EAS**: these populations have the smallest standard errors among the non-EUR populations, so even small differences from the other populations are flagged.
- **Heterogeneity in EUR**: if λ is high for EUR, EUR differs from the rescaled non-EUR estimates as a group. This can happen when the rescaling slope (unweighted OLS) is different from the IVW slope.

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

Paths are relative to the `script/` folder, so inputs are read from `../data/` and figures are written to `../figures/` wherever you render from.

### Settings

These are set in the setup chunk:

| Setting | Default | Meaning |
| --- | --- | --- |
| `populations` | AFR, AMR, EAS, EUR, SAS, UGR | Populations included |
| `alpha` | 0.05 | Significance level for Q, replication and Bonferroni |
| `replication_discovery_p` | 5e-8 | EUR p-value threshold for discovery SNPs (all SNPs) |
| `drug_replication_discovery_p` | 0.05 | EUR p-value threshold for discovery SNPs (drug-target subset) |
| `contribution_method` | `"loo"` | Per-population contribution method (`"loo"` or `"qj"`) |
| `figure_groups` | blood pressure; lipids | Traits shown in each cross-ancestry figure |
| `figure_dir` | `../figures` | Where the figure files are saved |
| `show_ci` | `TRUE` | Whether to draw 95% CIs in the cross-ancestry figures |
| `lambda_se_method` | `"bootstrap"` | How the SE of λ is calculated (`"bootstrap"` or `"asymptotic"`) |
| `lambda_boot_n` | 1000 | Number of bootstrap resamples for the SE of λ |
| `lambda_boot_seed` | 20261127 | Random seed for the bootstrap, so results can be reproduced |

### Requirements

Each input file must have the columns `snpid`, `pop`, `beta` and `se`. The replication analysis also needs `pval`. The drug-target analysis needs `Drug`, and the gene labels on its forest plots come from `Target gene`.

Required packages are `rmarkdown`, `knitr`, `dplyr`, `tidyr`, `ggplot2` and `scales`.

## Changes from the previous version

- **Replication power fixed.** In the expected-significance calculation, the second term computed the wrong tail of the normal distribution, so expected counts were close to zero. The observed replication p-value is now two-sided, to match.
- **Leave-one-out contributions.** These are now the default way to measure each population's contribution to heterogeneity.
- **New outputs.** The transferability table, cross-ancestry figures and saved figure files are new.
- **Uncertainty in λ.** Each λ now comes with an SE, a 95% CI and a test of λ = 1. These appear in the QQ plot legend, in a new λ plot and in a new table, and the QQ plot has a 95% null band.
- **Renamed columns.** In the replication failure table, `expected_pval` and `observed_pval` are now `expected_power` and `observed_sig`.
