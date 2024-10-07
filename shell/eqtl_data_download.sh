#!/bin/bash

# Download publicly available eQTL data


### eqtlgen

mkdir eqtl_data_eqtlgen

wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz -P eqtl_data_eqtlgen

wget
https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz -P eqtl_data_eqtlgen

