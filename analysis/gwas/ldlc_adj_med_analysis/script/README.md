# Adjusting observed LDL cholesterol values by medication status 

## Background 
In individuals taking LDL-lowering medication, their observed LDL cholesterol (LDLc) values will appear to be lower when truly these values would be high. Therefore, to restore their LDLc values to approximately their true values we perform relative and absolute medication adjustment 

## Simulation 
```

```

## Organise data 
Relative medication adjustment 
```
Rscript LDL_adj_by_med_status_rel.r
```
Absolute medication adjustment 
```
Rscript LDL_adj_by_med_status_abs.r
```

## Conduct GWAS 
Relative medication adjustment
```
Rscipt rel_LDLc_GWAS.r
```
Absolute medication adjustment 
```
Rscript abs_LDLc_GWAS.r
```
## LD clump GWAS to identify number of independent SNPs 
```

```
