# MRC Network project

This repo has the code for the MRC network multi-ancestry drug MR project.

## Organisation

- `data_generation` - use this folder for the per-study / per-ancestry phenotype harmonisation and GWAS analysis scripts
- `analysis` - use this folder for the upstream and downstream analysis components of the project (e.g. using the GWAS results for downstream analysis)

## Analysis 
### Genome-wide association analyses 
Phenotypes of interest: 
- LDL cholesterol
- HDL cholesterol
- Triglycerides
- Systolic blood pressure
- Diastolic blood pressure
- Pulse pressure
- Coronary heart disease
- Hypertension
- Stroke

| Data source  | Ancestry | GWAS software |
| --- | --- | --- |
| UK Biobank | AFR | GENESIS |
| | EAS| GENESIS  |
| | EUR | BOLT-LMM |
| | SAS | GENESIS |
| Ugandan Genome Resource | Continental AFR | GENESIS |

### Meta-analyses using METAL 
Analysis: 
- Inverse variance weighted approach weighted on sample size

### Comparison of effects across ancestry 
Analysis: 

