---
title: "LDL analysis - Adjustment of medication status"
output: html_notebook
---

# Background 
- There is currently no approach to fully adjust for medication status in lipid traits like LDL cholesterol. 
- Predominantly, studies have excluded participants based on medication use associated with the phenotype of interest. However, issues such as collider bias and measurement error can arise as a result of unadjustment or adjustment of the observed values. 

# Simulation 

Load relevant packages 
```{r}
library(dplyr)
library(ggplot2)
```

Generate data 
```{r}
nid <- 100000
g <- rnorm(nid)
g_other <- rep(0.1, nid)
mean_ldl <- 3.56
sd_ldl <- 0.87
med_effect <- 0.8
b_omed <- 0.5

generate_data <- function(nid, mean_ldl, sd_ldl,b_gy, b_gother) {
  mean_ldl <- mean_ldl
  sd_ldl <- sd_ldl
  ldl <- rnorm(nid, mean_ldl, sd_ldl)
  b_geno_y <- rnorm(nid, b_gy)
  b_gother <- rnorm(nid, b_gother)
    dat <- data.frame(
    nid = rep(1:nid, each = 1),
    mean_ldl = mean_ldl,
    sd_ldl = sd_ldl,
    b_gy = b_gy,
    b_gother = b_gother,
    ldl = ldl
  )
  
  return(dat)
}
```

Medication use - Calculate the probability of medication use based on prevalence 
```{r}
medication_prevalence <- 0.13

med_prob <- function(g_other, b_omed) {
  plogis(qlogis(medication_prevalence) + g_other * b_omed)
}

med_use_prob <- med_prob(g_other, b_omed)

medication_use <- rbinom(nid, 1, med_use_prob)

mean_medication_use <- mean(medication_use)
print(mean_medication_use)
```

Estimate genetic effects of LDL cholesterol
```{r}
est_effs <- function(dat, med_effect, medication_use) {
  y <- dat$ldl + g * dat$b_gy + rnorm(nrow(dat), sd = sqrt(dat$sd_ldl^2 - dat$b_gy^2))
  y_obs <- y
  y_obs[as.logical(medication_use)] <- y_obs[as.logical(medication_use)] * med_effect
  y_adj_true <- y_obs
  y_adj_true[as.logical(medication_use) & med_effect != 0] <- y_obs[as.logical(medication_use) & med_effect != 0] / med_effect
  y_adj_approx <- y_obs
  y_adj_approx[as.logical(medication_use)] <- y_obs[as.logical(medication_use)] + 20
  
  lm_y <- lm(y ~ g_other)
  lm_y_obs <- lm(y_obs ~ g_other)
  lm_y_adj_true <- lm(y_adj_true ~ g_other)
  lm_y_adj_approx <- lm(y_adj_approx ~ g_other)
  
  results <- data.frame(
    measure = c("y", "y_obs", "y_adj_true", "y_adj_approx"),
    coefficient = c(coef(lm_y)[2], coef(lm_y_obs)[2], coef(lm_y_adj_true)[2], coef(lm_y_adj_approx)[2])
  )
  
  return(results)
}

```

Create parameter grid 
```{r}
param <- expand.grid(
  nsim = 1:100,
  nid = 10000,
  b_gy = c(0.5, 1),
  b_gother = c(0.5, 1)
)
```

Apply data generation and estimation to each parameter combination
```{r}
all_results <- list()

all_results <- lapply(1:nrow(param), function(i) {
  d <- generate_data(
    nid = param$nid[i],
    mean_ldl = 3.55595,
    sd_ldl = 0.870798,
    b_gy = param$b_gy[i],
    b_gother = param$b_gother[i]
  )
  
  # Generate medication_use vector within each iteration
  med_use_prob <- med_prob(g_other, b_omed)
  medication_use <- rbinom(nid, 1, med_use_prob)
  
  r <- est_effs(d, med_effect, medication_use)
  
  # Set 'nsim' attribute for each data frame in the 'all_results' list
  attr(r, "nsim") <- param$nsim[i]
  
  return(r)
})

```

Organise results
```{r}
# Add nsim, b_gy, and b_gother as attributes to each data frame in the 'all_results' list
for (i in seq_along(all_results)) {
  attr(all_results[[i]], "nsim") <- param$nsim[i]
  attr(all_results[[i]], "Estimate") <- param$Estimate[i]
  attr(all_results[[i]], "Std. Error") <- param$Std.Error[i]
  attr(all_results[[i]], "measure") <- param$measure[i]
}

# Convert the list of results to a single data frame
all_results_df <- bind_rows(all_results)

# Calculate standard errors based on the number of simulations (nsim)
all_results_df <- all_results_df %>%
  group_by(measure) %>%
  mutate(se = sd(Estimate) / sqrt(attr(all_results_df, "nsim")))

# Calculate the mean estimated coefficient and standard error for each measure
summary_results <- all_results_df %>%
  group_by(measure) %>%
  summarise(mean_estimate = mean(Estimate),
            mean_se = mean(se),
            nsim = unique(attr(all_results_df, "nsim")))

# Order of measures for the forest plot
measure_order <- c("y", "y_obs", "y_adj_approx", "y_adj_true")

# Create a forest plot
ggplot(summary_results, aes(x = mean_estimate, y = measure)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbarh(aes(xmin = mean_estimate - 1.96 * mean_se, xmax = mean_estimate + 1.96 * mean_se), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", size = 1) +
  labs(title = "Forest Plot of Estimated Coefficients",
       subtitle = paste("Based on", unique(summary_results$nsim), "simulations"),
       x = "Estimated Coefficient",
       y = "Measure") +
  scale_y_discrete(limits = rev(measure_order)) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
```
