#!/bin/bash

# Process the GWAS data for Mendelian randomization analysis


###GLIDE

echo "replication_EUR_nteeth_excl_HCHSSOL
replication_EUR_DMFS_excl_HCHSSOL
replication_EUR_DFSS_excl_HCHSSOL" > GLIDE.txt


## generate read_outcome_data scripts for progression
while read OUTCOME; do
    cat ./mr_druggable_genome_pd/R/read_outcome_data_GLIDE.R > ./mr_druggable_genome_pd/R/read_outcome_data_${OUTCOME}.R
done < GLIDE.txt



###UKB

echo "EUR_dentures_prep
EUR_toothache_prep
EUR_ulcers_prep
EUR_painfulgums_prep
EUR_bleedinggums_prep
EUR_looseteeth_prep" > UKB.txt


## generate read_outcome_data scripts for progression
while read OUTCOME; do
    cat ./mr_druggable_genome_pd/R/read_outcome_data_UKB.R > ./mr_druggable_genome_pd/R/read_outcome_data_${OUTCOME}.R
done < UKB.txt


###Finngen

echo "finngen_R9_K11_periodontitis_prep
finngen_R9_K11_caries_prep
finngen_R9_K11_nTeeth_prep
replication_finngen_R9_K11_APHTA_RECUR_INCLAVO_prep
replication_EUR_perio_excl_HCHSSOL" > Finngen.txt


## generate read_outcome_data scripts for progression
while read OUTCOME; do
    cat ./mr_druggable_genome_pd/R/read_outcome_data_finngen.R > ./mr_druggable_genome_pd/R/read_outcome_data_${OUTCOME}.R
done < Finngen.txt



# pQTL data #占内存太大运行速度太慢，手动运行一下
#nohup Rscript ./mr_druggable_genome_pd/R/08data_prep_pqtl.R &> nohup_data_data_prep_pqtl.log &
