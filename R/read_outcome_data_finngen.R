

### load outcome data
out <- read_outcome_data(snps = exp$SNP,
                         filename = str_c("outcome_data/",OUTCOME,".txt"),
                         sep = "\t",
                         snp_col = 'SNP',
                         beta_col = "beta.outcome",
                         se_col = "se.outcome",
                         effect_allele_col = "effect_allele.outcome",
                         other_allele_col = "other_allele.outcome",
                         eaf_col = "eaf.outcome",
                         pval_col = "pval.outcome",
                         samplesize_col = "samplesize.outcome",
                         phenotype_col = "outcomedata"
)

out$outcome <- OUTCOME

