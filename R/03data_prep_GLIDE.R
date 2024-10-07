
library("dplyr")
library("readr")
library("stringr")
library("data.table")



folder_path <- "./outcome_data"
glide_files <- list.files(folder_path, pattern = "HCHSSOL.txt$", full.names = TRUE)
#file <- glide_files[4]
for (file in glide_files) {
# read in gwas data
pd_risk_discovery <- fread(file)
head(pd_risk_discovery)
colnames(pd_risk_discovery)[1] <- "chr_pos"
pd_risk_discovery$chr_pos <- apply(pd_risk_discovery, 1, function(row) paste0("chr", row["chr_pos"]))

pd_risk_discovery$outcome <- sub(".*/.*_.*_(.*?)_.*", "\\1", file)

# read in eqtl data
eqtlgen <- read.table("eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

eqtlgen <- distinct(eqtlgen[,c("SNP", "SNPChr", "SNPPos")])

names(eqtlgen) <- c("rsid_eqtlgen", "chr", "position")

eqtlgen$chr_pos <- str_c("chr",eqtlgen$chr, ":",eqtlgen$position)

# pqtl
#暂时无 pqtl <- read.table("pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt", sep = "\t",colClasses = "character", header = T)
# 
# pqtl$chr_pos <- str_c("chr", pqtl$chr, ":", pqtl$pos)
# names(pqtl)[1] <- "rsid_replication_pqtl"

# join to add rsids to the gwas data
pd_risk_discovery_1 <- left_join(pd_risk_discovery, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = "chr_pos")
#暂时无pd_risk_discovery_2 <- left_join(pd_risk_discovery_1, pqtl[,c("chr_pos", "rsid_replication_pqtl")], by = "chr_pos")

saved_file_name <- paste0(tools::file_path_sans_ext(file), "_prep.txt")

fwrite(pd_risk_discovery_1, file = saved_file_name, sep = "\t")
print(saved_file_name)
}


setwd('/home/piperacillin/projects/DENTALTRAITS')
glide <- fread('./outcome_data/newGLIDE/EUR_perio_excl_HCHSSOL.txt_matched.txt')
head(glide)

library(TwoSampleMR)
exp0 <- read_exposure_data(
  filename = "eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  phenotype_col = "GeneSymbol",
  samplesize_col = "NrSamples",
  min_pval = 1e-400
)

out_dat2 <- read_outcome_data(snps = exp0$SNP,
  filename = './outcome_data/newGLIDE/EUR_perio_excl_HCHSSOL.txt_matched.txt', ##你的数据
  sep= "\t",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col ="Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  samplesize_col = "N"
)


glide1<-glide[glide$`P-value` < 0.05, ]

out_dat <- format_data(glide1,
                        type = 'exposure',
                        snp_col = "SNP",
                        beta_col = "Effect",
                        se_col = "StdErr",
                        effect_allele_col ="Allele1",
                        other_allele_col = "Allele2",
                        pval_col = "P-value",
                        samplesize_col = "N")

source('./rawinstruments/geteaffrom1000G.R')

newdat <- get_eaf_from_1000G(out_dat, "./rawinstruments", type = "exposure")
setwd('/home/piperacillin/projects/DENTALTRAITS')
head(newdat)

head(glide)
merged_df <- merge(glide, newdat[c("SNP", "eaf.exposure")], by.x = "SNP",by.y = "SNP", all.x = TRUE)
head(merged_df)
names(merged_df)[names(merged_df) == "eaf.exposure"] <- "Freq1"

newdat2 <- get_eaf_from_1000G(out_dat2, "./rawinstruments", type = "exposure") #整合暴露数据
head(newdat2)
setwd('/home/piperacillin/projects/DENTALTRAITS')
# 合并两个数据框
merged_df2 <- merge(merged_df, newdat2[, c("SNP", "eaf.outcome")], by = "SNP", all.x = TRUE)

# 使用eaf.outcome的值填充Freq1中的NA值
merged_df2$Freq1 <- ifelse(is.na(merged_df2$Freq1), merged_df2$eaf.outcome, merged_df2$Freq1)

# 移除eaf.outcome列
merged_df2 <- subset(merged_df2, select = -eaf.outcome)

head(merged_df2)
sum(!is.na(newdat$eaf.exposure))
sum(!is.na(merged_df2$Freq1))


fwrite(merged_df2, file = './outcome_data/newGLIDE/EUR_perio_excl_HCHSSOL.txt_matched_EAF.txt', sep = "\t")
