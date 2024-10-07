rm(list=ls())
library("dplyr")
library("readr")
library("stringr")
library("data.table")
library(TwoSampleMR)
getwd()

finn_info <- fread('./Finn_R9_data.csv',data.table = F)

# read in eqtl data – will get rs IDs for the SNPs from here
eqtlgen <- read.table("eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

eqtlgen <- distinct(eqtlgen[,c("SNP", "SNPChr", "SNPPos")])

names(eqtlgen) <- c("rsid_eqtlgen", "chr", "position")

eqtlgen$chr_pos <- str_c("chr",eqtlgen$chr, ":",eqtlgen$position)


folder_path <- "./outcome_data"
finn_files <- list.files(folder_path, pattern = "^finngen_R9_.*\\.gz$", full.names = TRUE)


for (file in finn_files) {
# read in replication risk data
replication_data_finn <- fread(file,data.table = F)


finn <- replication_data_finn
finndata <- format_data(finn,type = 'outcome',snp_col = "rsids", beta_col = "beta", se_col = "sebeta", 
                        eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref", 
                        pval_col = "pval", gene_col = "nearest_genes", chr_col = "#chrom", 
                        pos_col = "pos")
trait_row <- finn_info[grepl(paste0(tools::file_path_sans_ext(basename(file)), ".gz"),
                             finn_info$path_https),]

finndata$ncase.outcome <- trait_row$num_cases
finndata$ncontrol.outcome <- trait_row$num_controls
finndata$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls

outcome <- trait_row$name
finndata$outcome <- outcome

replication_data_finn <- finndata
rm(finn)
rm(finndata)
replication_data_finn$chr_pos <- str_c("chr",replication_data_finn$chr.outcome, ":",replication_data_finn$pos.outcome)
replication_data_finn$outcomedata <- sub(".*_(.*)", "\\1", sub("\\.txt", "", file))
head(replication_data_finn)

# join to add rsids to the gwas data
replication_data_finn <- left_join(replication_data_finn, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = "chr_pos")

saved_file_name <- paste0(tools::file_path_sans_ext(file), "_prep.txt")

fwrite(replication_data_finn, file = saved_file_name, sep = "\t")
print(saved_file_name)
}


print("mission complete")
