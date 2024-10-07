

### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                         filename = str_c("outcome_data/",OUTCOME,".txt"),
                         sep = "\t",
                         snp_col = 'SNP',
                         beta_col = "BETA",
                         se_col = "SE",
                         effect_allele_col = "ALLELE1",
                         other_allele_col = "ALLELE0",
                         pval_col = "P_BOLT_LMM_INF",
                         samplesize_col = "total_n"
                         )
out$outcome <- OUTCOME