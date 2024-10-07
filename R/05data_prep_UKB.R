
library("dplyr")
library("readr")
library("stringr")
library("data.table")

setwd('G:/rawdata/MR/DENTALTRAITS')
# # read in allele data
# alleles <- read.table(gzfile("outcome_data/reference.txt.gz"), sep = ",", stringsAsFactors = F, header = T)
# names(alleles) <- c("chrpos", "rsid", "chr", "start", "other_allele", "effect_allele", "maf", "func", "near_gene")

# read in QTL data

# eqtlgen
eqtlgen <- read.table("eqtl_data_eqtlgen/eqtlgen_exposure_dat_snps_5kb_window.txt", sep = "\t",colClasses = "character", header = T)

eqtlgen <- distinct(eqtlgen[,c("SNP", "SNPChr", "SNPPos")])

names(eqtlgen) <- c("rsid_eqtlgen", "chr", "position")

eqtlgen$chr_pos <- str_c(eqtlgen$chr, ":",eqtlgen$position)


# ## pqtl
#暂时无 pqtl <- read.table("pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt", sep = "\t",colClasses = "character", header = T)
# 
# pqtl$chr_pos <- str_c(pqtl$chr, ":", pqtl$pos)
# names(pqtl)[1] <- "rsid_replication_pqtl"


# read in  each file and populate a data frame
folder_path <- "./outcome_data"
file_list_res <- list.files(folder_path, 
                            pattern = "^EUR", 
                            ignore.case = TRUE, 
                            full.names = TRUE)[!grepl("HCHSSOL", 
                                                      list.files(folder_path, 
                                                                 pattern = "^EUR", 
                                                                 ignore.case = TRUE))]

for (file in file_list_res) {
  raw_data <- fread(file,data.table = F)
  
  raw_data$outcome <- gsub(".*/EUR_(.*)\\.txt", "\\1", file)
  raw_data$chr_pos <- str_c("chr",raw_data$CHR, ":",raw_data$BP)
  
  raw_data_with_alleles_1 <- left_join(raw_data, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = "chr_pos")
  
  #暂时无raw_data_with_alleles_1 <- left_join(raw_data_with_alleles_1, pqtl[,c("chr_pos", "rsid_replication_pqtl")], by = c("SNP" = "chr_pos"))
  
  saved_file_name <- paste0(tools::file_path_sans_ext(file), "_prep.txt")
  
  fwrite(raw_data_with_alleles_1, file = saved_file_name, sep = "\t")
  print(saved_file_name)
}

print("mission complete")
