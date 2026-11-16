# Cross-population heterogeneity reports

`heterogeneity_exposures.Rmd` discovers all RDS inputs in the sibling `data/` directory whose filenames end in `_pop_extract_common_snps_clumped_corrected.rds` and produces one HTML report covering each exposure. Currently available traits are taken from the filename prefix, so newly added matching files are included automatically.

To analyse only selected exposures, edit the YAML parameter at the top of the Rmd, for example:

```yaml
params:
  traits: ["sbp", "dbp"]
```

To analyse all matching exposures, set `traits: null` (the default). Population labels can be changed in the `populations` vector in the setup chunk. The report checks for the columns `snpid`, `pop`, `beta`, and `se`; drug-specific summaries are produced only when a `Drug` column is available.

Render from the repository root with:

```r
rmarkdown::render(
  "analysis/heterogeneity_analyses_exposure/script/heterogeneity_exposures.Rmd",
  output_file = "heterogeneity_exposures.html"
)
```

The report requires `rmarkdown`, `knitr`, `dplyr`, `tidyr`, and `ggplot2`. Rendering requires a local R installation and the input RDS files; the HTML report is created locally by the render command.
