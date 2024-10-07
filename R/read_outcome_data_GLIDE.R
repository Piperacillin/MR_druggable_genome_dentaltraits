

### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                         filename = str_c("outcome_data/",OUTCOME,".txt"),
                         sep = "\t",
                         snp_col = 'SNP',
                         beta_col = "Effect",
                         se_col = "StdErr",
                         effect_allele_col = "Allele1",
                         other_allele_col = "Allele2",
                         pval_col = "P-value",
                         eaf_col = "Freq1",
                         samplesize_col = "N"
)
out$outcome <- OUTCOME