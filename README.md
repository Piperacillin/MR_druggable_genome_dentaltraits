# MR_druggable_genome_dentaltraits
Adapted from the original publication from a Mendelian randomization analysis for the druggable genome in Parkinson's disease.(https://www.nature.com/articles/s41467-021-26280-1).

20231216 Update: The manuscript is currently being prepared for submission to a peer-reviewed journal. Upon formal publication, the scripts will be made available to the readers.

## Pipeline Overview


1. Install the required tools

```bash
bash ./mr_druggable_genome_pd/shell/installing_tools.sh
```


2. Download/upload your QTL and disease GWAS data. eQTL data from the PsychENCODE and eQTLGen consortia can be accessed as shown below. If you are looking for a comprehensive list of pQTL data available, check out [this very helpful blog post](http://www.metabolomix.com/a-table-of-all-published-gwas-with-proteomics/).

```bash
bash ./mr_druggable_genome_pd/shell/eqtl_data_download.sh
```

3. Select exposure data and outcome data.

Specify your exposure data. Example below.
```bash
echo "eqtlgen" > exposure_data.txt
```


For discovery analyses, specify your outcomes. Example below.
```bash
echo "EUR_dentures_prep
EUR_toothache_prep
EUR_ulcers_prep
EUR_painfulgums_prep
EUR_bleedinggums_prep
EUR_looseteeth_prep
finngen_R9_K11_caries_prep
finngen_R9_K11_nTeeth_prep
finngen_R9_K11_periodontitis_prep" > outcomes.txt

echo "" > discovery_outcomes.txt
```


4. Prepare the data for the Mendelian randomization analysis. Note: if you are using your own GWAS data you will need to edit these files and create appropriate `read_exposure_data` and `read_outcome_data` files
```bash
bash ./mr_druggable_genome_pd/shell/06data_prep.sh
```

5. Generate scripts that can be run in parallel
```bash

while read EXPOSURE_DATA; do
    while read OUTCOME; do

            export EXPOSURE_DATA=${EXPOSURE_DATA}
            export OUTCOME=${OUTCOME}
            mkdir ${EXPOSURE_DATA}_${OUTCOME}

            while read DISCOVERY_OUTCOME; do
                    export DISCOVERY_OUTCOME=${DISCOVERY_OUTCOME}
                    bash ./mr_druggable_genome_pd/shell/07generate_parallel_scripts.sh
            done < discovery_outcomes.txt

    done < outcomes.txt
done < exposure_data.txt

```

6. Run Mendelian randomization for the druggable genome using all the exposure data and outcome data you have specified.
```bash
echo "exposure_data,exposure,outcome" > exposures_to_remove.txt

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ./mr_druggable_genome_pd/shell/run_liberal_scripts_all_nohup.sh &> ./mr_druggable_genome_pd/shell/nohup_run_liberal_scripts_all.log &
    done < outcomes.txt
done < exposure_data.txt

```

7. Remove genes that throw errors. Note: you may need to rerun this step a few times until all genes that cause an error are removed.
```bash
nohup bash ./mr_druggable_genome_pd/shell/run_liberal_scripts_failed_nohup.sh &> ./mr_druggable_genome_pd/shell/nohup_run_liberal_scripts_failed.log &
```

8. For each exposure-data-outcome combination, put all the results into one results file.
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        cd ./${EXPOSURE_DATA}_${OUTCOME}/results
        nohup Rscript ../../mr_druggable_genome_pd/R/combine_results_liberal_r2_0.2.R &> ../../${EXPOSURE_DATA}_${OUTCOME}/nohup_combine_results_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.log &
        cd ../..
    done < outcomes.txt
done < exposure_data.txt
```

9. Rerun the analysis for all significant genes using the clumping threshold r2 = 0.001.
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ./${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.sh &> ./${EXPOSURE_DATA}_${OUTCOME}/nohup_script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.log &
    done < outcomes.txt
done < exposure_data.txt
```



10.For replication analyses, the outcome should begin with "replication_", and you need to specify your discovery outcome. This will make sure your pvalue threshold for significance is the non-adjusted pvalue, and that you only run the replication step for discovered genes. Examples below.
```bash
echo "replication_finngen_R9_ulcers_prep"  > outcomes.txt

echo "EUR_ulcers_prep" > discovery_outcomes.txt
```

```bash
echo "replication_EUR_perio"  > outcomes.txt

echo "finngen_R9_K11_periodontitis_prep" > discovery_outcomes.txt
```


```bash
echo "replication_EUR_DMFS
replication_EUR_DFSS"  > outcomes.txt

echo "finngen_R9_K11_caries_prep" > discovery_outcomes.txt
```

For pQTL replication
```bash
echo "replication_pqtl" > exposure_data.txt
echo "EUR_dentures_prep
EUR_toothache_prep
EUR_ulcers_prep
EUR_painfulgums_prep
EUR_bleedinggums_prep
EUR_looseteeth_prep
finngen_R9_K11_caries_prep
finngen_R9_K11_nTeeth_prep
finngen_R9_K11_periodontitis_prep" > outcomes.txt
echo "" > discovery_outcomes.txt
```

11. Prepare the data for the Mendelian randomization analysis. Note: if you are using your own GWAS data you will need to edit these files and create appropriate `read_exposure_data` and `read_outcome_data` files
```bash
bash ./mr_druggable_genome_pd/shell/06data_prep.sh
```

12. Generate scripts that can be run in parallel
```bash

while read EXPOSURE_DATA; do
    while read OUTCOME; do

            export EXPOSURE_DATA=${EXPOSURE_DATA}
            export OUTCOME=${OUTCOME}
            mkdir ${EXPOSURE_DATA}_${OUTCOME}

            while read DISCOVERY_OUTCOME; do
                    export DISCOVERY_OUTCOME=${DISCOVERY_OUTCOME}
                    bash ./mr_druggable_genome_pd/shell/07generate_parallel_scripts.sh
            done < discovery_outcomes.txt

    done < outcomes.txt
done < exposure_data.txt

```

13. Run Mendelian randomization for the druggable genome using all the exposure data and outcome data you have specified.
```bash
echo "exposure_data,exposure,outcome" > exposures_to_remove.txt

while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ./mr_druggable_genome_pd/shell/run_liberal_scripts_all_nohup.sh &> ./mr_druggable_genome_pd/shell/nohup_run_liberal_scripts_all.log &
    done < outcomes.txt
done < exposure_data.txt

```

14. Remove genes that throw errors. Note: you may need to rerun this step a few times until all genes that cause an error are removed.
```bash
nohup bash ./mr_druggable_genome_pd/shell/run_liberal_scripts_failed_nohup.sh &> ./mr_druggable_genome_pd/shell/nohup_run_liberal_scripts_failed.log &

15. For each exposure-data-outcome combination, put all the results into one results file.
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        cd ./${EXPOSURE_DATA}_${OUTCOME}/results
        nohup Rscript ../../mr_druggable_genome_pd/R/combine_results_liberal_r2_0.2.R &> ../../${EXPOSURE_DATA}_${OUTCOME}/nohup_combine_results_liberal_r2_0.2_${EXPOSURE_DATA}_${OUTCOME}.log &
        cd ../..
    done < outcomes.txt
done < exposure_data.txt
```

16. Rerun the analysis for all significant genes using the clumping threshold r2 = 0.001.
```bash
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup bash ./${EXPOSURE_DATA}_${OUTCOME}/script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.sh &> ./${EXPOSURE_DATA}_${OUTCOME}/nohup_script_conservative_r2_0.001_${EXPOSURE_DATA}_${OUTCOME}.log &
    done < outcomes.txt
done < exposure_data.txt
```

17. Some data formatting steps.
```bash

mkdir full_results

bash ./mr_druggable_genome_pd/shell/final_results_report.sh

Rscript ./mr_druggable_genome_pd/R/combine_dat_steiger_pQTL.R

Rscript ./mr_druggable_genome_pd/R/format_supplement_pQTL.R



bash ./mr_druggable_genome_pd/shell/final_results_report.sh

bash ./mr_druggable_genome_pd/shell/final_results_report.sh

bash ./mr_druggable_genome_pd/shell/final_results_report.sh



Rscript ./mr_druggable_genome_pd/R/combine_dat_steiger.R


Rscript ./mr_druggable_genome_pd/R/format_supplement.R



```

18. Check if the direction of effect is consistent between any discovery cohorts and replication cohorts.


```bash
while read OUTCOME; do
    export OUTCOME=${OUTCOME}

    while read DISCOVERY_OUTCOME; do
        export DISCOVERY_OUTCOME=${DISCOVERY_OUTCOME}
        nohup Rscript ./mr_druggable_genome_pd/R/check_direction_of_effect.R &> full_results/metric_check_direction_of_effect_${OUTCOME}_${DISCOVERY_OUTCOME}.log
    done < discovery_outcomes.txt

    
done < outcomes.txt
```


19. Generate forest plots. This script is not generic for any exposure/outcome data.

```bash
mkdir figures

Rscript ./mr_druggable_genome_pd/R/make_forest_plot.R
```

20. Colocalization analysis.

```bash
mkdir coloc

Rscript ./mr_druggable_genome_pd/R/data_prep_coloc.R

Rscript ./mr_druggable_genome_pd/R/run_coloc.R


```



