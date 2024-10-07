rm(list=ls())
library(TwoSampleMR)
library(data.table)


folder_path <- "./"
gz_files <- list.files(folder_path, pattern = ".gz$", full.names = TRUE)
finn_info <- fread('./Finn_R9_data.csv',data.table = F)

#paste0(tools::file_path_sans_ext(basename(gz_files)), ".gz") basename截取最后一段

# 循环读取并保存成CSV文件
for (file in gz_files) {
  
  finn <- fread(file,data.table = F)
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
  
  # 构建保存的文件路径
  rdata_file <- paste0(tools::file_path_sans_ext(file), ".txt")
  
  # 保存为R.data
  fwrite(finndata, file = rdata_file, sep = "\t")
}

