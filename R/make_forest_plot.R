## this script will put the results and quality control metrics for all significant outcomes into one data frame
getwd()

# prep
library(dplyr)
library(readr)
library(forestplot)
library(stringr)
library(tidyverse)
rm(list=ls())
significant_res <- as.data.frame(read_tsv("full_results/significant_genes_results_all_outcomes.txt"))
significant_res <- unique(significant_res)

head(significant_res)

#### make forest plot Ulcers ####
discovery <- significant_res[significant_res$outcome == "EUR_ulcers_prep",]
discovery <- discovery[!grepl("_",discovery$exposure),]
discovery$exposure_tissue <- str_c(discovery$exposure, "_", discovery$tissue)
replication <- significant_res[significant_res$outcome == "replication_finngen_R9_K11_APHTA_RECUR_INCLAVO_prep",]
replication$exposure_tissue <- str_c(replication$exposure, "_", replication$tissue)
discovery <- subset(discovery, discovery$exposure_tissue %in% replication$exposure_tissue)


data <- as.data.frame(subset(discovery, discovery$clump_thresh == "0.2"))
data$outcome <- tools::toTitleCase(gsub("^.*_(.*?)_.*$", "\\1", data$outcome)) #提取中间部分并大写
data <- data[order(data$exposure),]

data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
#data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
#data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"
data$se <- (log(data$or_uci95) - log(data$or))/1.96
data$or_ci <- ifelse(is.na(data$se), "",
                     sprintf("%.2f (%.2f to %.2f)",
                             data$or, data$or_lci95, data$or_uci95))
data$fdr_qval <- sprintf("%.2f", data$fdr_qval)
#data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers_p)

colnames(data)[colnames(data) %in% c("exposure", "nsnp", "or_ci", "fdr_qval", "druggability_tier")] <- c("Gene", "No. SNPs","FDR-corrected P", "Druggability tier", "OR (95%CI)")

data$` ` <- paste(rep(" ", 20), collapse = " ")

library(grid)
library(forestploter)

tm <- forest_theme(base_size = 10,
                   refline_col = "grey20",
                   ci_col = c("darkred"),
                   arrow_type = "closed",
                   footnote_col = "royalblue",
                   ci_lwd = 1.5,
                   ci_Theight = 0.2)


pdf("figures/forest_ulcers.pdf", width = 10, height = 4, onefile=FALSE)
forest(data[,c(1,2, 3 , 23, 22, 9, 11)],
       est = data$or,
       lower = data$or_lci95, 
       upper = data$or_uci95,
       sizes = 0.3,
       ci_column = 4,
       ref_line = 1,
       arrow_lab = c("Protective effect", "Harmful effect"),
       xlim = c(0.95, 1.05),
       ticks_at = c(0.95, 1, 1.05),
       #footnote = "This is the demo data. Please feel free to change\nanything you want.",
       theme = tm)
dev.off()



#### make forest plot periodontitis ####
discovery <- significant_res[significant_res$outcome == "finngen_R9_K11_periodontitis_prep",]
discovery <- discovery[!grepl("_",discovery$exposure),]
discovery$exposure_tissue <- str_c(discovery$exposure, "_", discovery$tissue)
replication <- significant_res[significant_res$outcome == "replication_EUR_perio_excl_HCHSSOL",]
replication$exposure_tissue <- str_c(replication$exposure, "_", replication$tissue)
discovery <- subset(discovery, discovery$exposure_tissue %in% replication$exposure_tissue)


data <- as.data.frame(subset(discovery, discovery$clump_thresh == "0.2"))
data$outcome <- sapply(strsplit(data$outcome, "_"), function(x) x[4])
data$outcome <- tools::toTitleCase(data$outcome)

data <- data[order(data$exposure),]

data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
#data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
#data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"
data$se <- (log(data$or_uci95) - log(data$or))/1.96
data$or_ci <- ifelse(is.na(data$se), "",
                     sprintf("%.2f (%.2f to %.2f)",
                             data$or, data$or_lci95, data$or_uci95))
data$fdr_qval <- sprintf("%.2f", data$fdr_qval)
#data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers_p)

colnames(data)[colnames(data) %in% c("exposure", "nsnp", "or_ci", "fdr_qval", "druggability_tier")] <- c("Gene", "No. SNPs","FDR-corrected P", "Druggability tier", "OR (95%CI)")
data$` ` <- paste(rep(" ", 20), collapse = " ")

library(grid)
library(forestploter)

tm <- forest_theme(base_size = 10,
                   refline_col = "grey20",
                   ci_col = c("darkred"),
                   arrow_type = "closed",
                   footnote_col = "royalblue",
                   ci_lwd = 1.5,
                   ci_Theight = 0.2)


pdf("figures/forest_periodontitis.pdf", width = 10, height = 4, onefile=FALSE)
forest(data[,c(1,2, 3 , 23, 22, 9, 11)],
       est = data$or,
       lower = data$or_lci95, 
       upper = data$or_uci95,
       sizes = 0.3,
       ci_column = 4,
       ref_line = 1,
       arrow_lab = c("Protective effect", "Harmful effect"),
       xlim = c(0.5, 1.50),
       ticks_at = c(0.5, 0.75, 1, 1.25,1.5),
       #footnote = "This is the demo data. Please feel free to change\nanything you want.",
       theme = tm)
dev.off()


#### make forest plot other traits ####
# 选择outcome列非replication开头的行
data_forest <- subset(significant_res, !grepl("^replication", outcome))

# 排除EUR_ulcers_prep或finngen_R9_K11_periodontitis_prep的行
data_forest <- subset(data_forest, !(outcome %in% c("EUR_ulcers_prep", "finngen_R9_K11_periodontitis_prep")))

data_forest$outcome <- tools::toTitleCase(gsub("^.*_(.*?)_.*$", "\\1", data_forest$outcome)) #提取中间部分并大写

data_forest <- data_forest[!grepl("_",data_forest$exposure),]
data <- as.data.frame(subset(data_forest, data_forest$clump_thresh == "0.2"))
#data[data$outcome == "pd_age_at_onset","outcome"] <- "Age at onset"
data <- data[order(data$outcome, data$exposure),]

data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]

data$se <- (log(data$or_uci95) - log(data$or))/1.96
data$or_ci <- ifelse(is.na(data$se), "",
                     sprintf("%.2f (%.2f to %.2f)",
                             data$or, data$or_lci95, data$or_uci95))
data$fdr_qval <- sprintf("%.2f", data$fdr_qval)

colnames(data)[colnames(data) %in% c("exposure", "nsnp", "or_ci", "fdr_qval", "druggability_tier")] <- c("Gene", "No. SNPs","FDR-corrected P", "Druggability tier", "OR (95%CI)")
head(data)
data$` ` <- paste(rep(" ", 20), collapse = " ")

library(grid)
library(forestploter)

tm <- forest_theme(base_size = 10,
                   refline_col = "grey20",
                   ci_col = c("darkred"),
                   arrow_type = "closed",
                   footnote_col = "royalblue",
                   ci_lwd = 1.5,
                   ci_Theight = 0.2)


pdf("figures/forest_othertraits.pdf", width = 10, height = 10, onefile=FALSE)
forest(data[,c(1,2, 3, 22, 21, 9, 11)],
       est = data$or,
       lower = data$or_lci95, 
       upper = data$or_uci95,
       sizes = 0.3,
       ci_column = 4,
       ref_line = 1,
       arrow_lab = c("Protective effect", "Harmful effect"),
       xlim = c(0.85, 1.10),
       ticks_at = c(0.85, 0.9, 0.95, 1, 1.05,1.1),
       #footnote = "This is the demo data. Please feel free to change\nanything you want.",
       theme = tm)
dev.off()




#### PQTL DATA ####
data <- as.data.frame(read_tsv("full_results/significant_genes_results_all_outcomes_pQTL.txt"))
data <- as.data.frame(subset(data, data$clump_thresh == "0.2"))
data$outcome <- tools::toTitleCase(gsub("^.*_(.*?)_.*$", "\\1", data$outcome)) #提取中间部分并大写

data <- data[order(data$outcome, data$exposure),]
data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]

data$pqtl_study0 <- gsub(".*\\_", "", data$exposure)
data$pqtl_study[which(data$pqtl_study0 == "sun2018")] <- "Sun et al. 2018"
data$pqtl_study[which(data$pqtl_study0 == "suhre2017")] <- "Suhre et al. 2017"
data$pqtl_study[which(data$pqtl_study0 == "pietzner2021")] <- "Pietzner et al. 2021"
data$pqtl_study[which(data$pqtl_study0 == "gudjonsson2022")] <- "Gudjonsson et al. 2022"

data$gene <- gsub("\\_.*", "", data$exposure)

data$se <- (log(data$or_uci95) - log(data$or))/1.96
data$or_ci <- ifelse(is.na(data$se), "",
                     sprintf("%.2f (%.2f to %.2f)",
                             data$or, data$or_lci95, data$or_uci95))
data$fdr_qval <- sprintf("%.2f", data$fdr_qval)

colnames(data)[colnames(data) %in% c("nsnp", "or_ci", "fdr_qval", "druggability_tier", "pqtl_study", "gene")] <- c("No. SNPs","FDR-corrected P", "Druggability tier", "pQTL source", "Gene", "OR (95%CI)")
head(data)
data$` ` <- paste(rep(" ", 20), collapse = " ")


data <- data[order(data$outcome),]

library(grid)
library(forestploter)

tm <- forest_theme(base_size = 10,
                   refline_col = "grey20",
                   ci_col = c("royalblue"),
                   arrow_type = "closed",
                   footnote_col = "royalblue",
                   ci_lwd = 1.5,
                   ci_Theight = 0.2)


pdf("figures/forest_pQTL_replication.pdf", width = 10, height = 10, onefile=FALSE)
forest(data[,c(19, 18, 2, 3, 21, 20, 9)],
       est = data$or,
       lower = data$or_lci95, 
       upper = data$or_uci95,
       sizes = 0.3,
       ci_column = 5,
       ref_line = 1,
       arrow_lab = c("Protective effect", "Harmful effect"),
       xlim = c(0.95, 1.05),
       ticks_at = c(0.95, 1, 1.05),
       #footnote = "This is the demo data. Please feel free to change\nanything you want.",
       theme = tm)
dev.off()




dt <- read.csv(system.file("extdata", "example_data.csv", package = "forestploter"))
dt <- dt[1:7, ]
# Indent the subgroup if there is a number in the placebo column
dt$Subgroup <- ifelse(is.na(dt$Placebo), 
                      dt$Subgroup,
                      paste0("   ", dt$Subgroup))

# NA to blank or NA will be transformed to carachter.
dt$n1 <- ifelse(is.na(dt$Treatment), "", dt$Treatment)
dt$n2 <- ifelse(is.na(dt$Placebo), "", dt$Placebo)

# Add two blank columns for CI
dt$`CVD outcome` <- paste(rep(" ", 20), collapse = " ")
dt$`COPD outcome` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
dt$ci1 <- paste(sprintf("%.1f (%.1f, %.1f)", dt$est_gp1, dt$low_gp1, dt$hi_gp1),
                sprintf("%.1f (%.1f, %.1f)", dt$est_gp3, dt$low_gp3, dt$hi_gp3),
                sep = "\n")
dt$ci1[grepl("NA", dt$ci1)] <- "" # Any NA to blank

dt$ci2 <- paste(sprintf("%.1f (%.1f, %.1f)", dt$est_gp2, dt$low_gp2, dt$hi_gp2),
                sprintf("%.1f (%.1f, %.1f)", dt$est_gp4, dt$low_gp4, dt$hi_gp4),
                sep = "\n")
dt$ci2[grepl("NA", dt$ci2)] <- ""



data <- data.frame(
  exposure = rep(c("Exposure 1", "Exposure 2", "Exposure 3"), each = 3),
  Method = rep(c("IVW", "Egger", "MaxLik"), times = 3),
  beta = c(0.2, -0.1, 0.5, 0.3, -0.2, 0.6, 0.4, -0.3, 0.7),
  uci95 = c(0.25, -0.05, 0.55, 0.35, -0.15, 0.65, 0.45, -0.25, 0.75),
  lci95 = c(0.15, -0.15, 0.45, 0.25, -0.25, 0.55, 0.35, -0.35, 0.65),
  p = c(0.01, 0.05, 0.001, 0.02, 0.04, 0.002, 0.03, 0.03, 0.003)
)

# 使用 tidyr 包中的 pivot_wider 函数将数据框转换为宽格式
data_wide <- data %>%
  pivot_wider(names_from = Method, values_from = c(beta, uci95, lci95, p), names_glue = "{.value}_{Method}")

# 选择 Method 为 IVW 的行
data_ivw <- data %>%
  filter(Method == "IVW") %>%
  select(-Method)

# 将宽格式的数据与 Method 为 IVW 的数据合并
final_data <- data_ivw %>%
  left_join(data_wide, by = "exposure")
