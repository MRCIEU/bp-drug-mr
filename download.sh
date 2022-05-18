#!/bin/bash

dl="/mnt/storage/scratch/gh13047/data/bp-drug-mr"
mkdir -p ${dl}/raw
mkdir -p ${dl}/extract

# MESA eQTL data
# from https://www.dropbox.com/sh/f6un5evevyvvyl9/AAA3sfa1DgqY67tx4q36P341a?dl=0
wget https://www.dropbox.com/sh/f6un5evevyvvyl9/AACdHg6dZOHTkgMEkC83MfZKa/AFA_cis_eqtl_summary_statistics.txt.gz?dl=1 -O ${dl}/raw/AFA_cis_eqtl_summary_statistics.txt.gz
wget https://www.dropbox.com/sh/f6un5evevyvvyl9/AAA6CDpnzSu-wjZiRNHyE3q-a/AFHI_cis_eqtl_summary_statistics.txt.gz?dl=1 -O ${dl}/raw/AFHI_cis_eqtl_summary_statistics.txt.gz
wget https://www.dropbox.com/sh/f6un5evevyvvyl9/AACA18cEUeKmaXvr6rwkrVXEa/ALL_cis_eqtl_summary_statistics.txt.gz?dl=1 -O ${dl}/raw/ALL_cis_eqtl_summary_statistics.txt.gz
wget https://www.dropbox.com/sh/f6un5evevyvvyl9/AADPa1-ZVcudv9U57CY5mxK-a/CAU_cis_eqtl_summary_statistics.txt.gz?dl=1 -O ${dl}/raw/CAU_cis_eqtl_summary_statistics.txt.gz
wget https://www.dropbox.com/sh/f6un5evevyvvyl9/AADHdeIGfdWmwNwWfLkyaD9wa/HIS_cis_eqtl_summary_statistics.txt.gz?dl=1 -O ${dl}/raw/HIS_cis_eqtl_summary_statistics.txt.gz

# ARIC pQTL data
wget https://jhupwas.s3.amazonaws.com/summary_data/EA.zip -O ${dl}/raw/EA.zip
wget https://jhupwas.s3.amazonaws.com/summary_data/AA.zip -O ${dl}/raw/AA.zip
unzip ${dl}/raw/EA.zip
mv seqid.txt EA
mv EA ${dl}/raw

unzip ${dl}/raw/AA.zip
mv seqid.txt AA
mv AA ${dl}/raw

## Extract

zgrep -wf genes.txt ${dl}/raw/AFA_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/AFA_cis_eqtl_summary_statistics.txt
zgrep -wf genes.txt ${dl}/raw/AFHI_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/AFHI_cis_eqtl_summary_statistics.txt
zgrep -wf genes.txt ${dl}/raw/ALL_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/ALL_cis_eqtl_summary_statistics.txt
zgrep -wf genes.txt ${dl}/raw/CAU_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/CAU_cis_eqtl_summary_statistics.txt
zgrep -wf genes.txt ${dl}/raw/HIS_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/HIS_cis_eqtl_summary_statistics.txt

zgrep -wf gene_alias.txt ${dl}/raw/AFA_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/AFA_cis_eqtl_summary_statistics_alias.txt
zgrep -wf gene_alias.txt ${dl}/raw/AFHI_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/AFHI_cis_eqtl_summary_statistics_alias.txt
zgrep -wf gene_alias.txt ${dl}/raw/ALL_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/ALL_cis_eqtl_summary_statistics_alias.txt
zgrep -wf gene_alias.txt ${dl}/raw/CAU_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/CAU_cis_eqtl_summary_statistics_alias.txt
zgrep -wf gene_alias.txt ${dl}/raw/HIS_cis_eqtl_summary_statistics.txt.gz > ${dl}/extract/HIS_cis_eqtl_summary_statistics_alias.txt



