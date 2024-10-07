# pqtl data prep for druggable genome replication


library(dplyr)
library(stringr)
library(tidyverse)

#setwd('/home/piperacillin/projects/DENTALTRAITS')
getwd()
# genes to replicate
druggable_genome <- read.csv("druggable_genome_new.txt", sep = "\t", header = T, stringsAsFactors = F)
sign <- as.data.frame(read_tsv("full_results/significant_genes_results_all_outcomes.txt"))
sign <- subset(sign, sign$outcome != "pd_risk_discovery") #后面要改
unique(sign$exposure)

#druggable_genome_sign <- druggable_genome
druggable_genome_sign <- subset(druggable_genome, druggable_genome$gene_display_label %in% sign$exposure)

length(unique(druggable_genome_sign$gene_display_label))


# read in pqtl data
# obtained from supplementary material of the relevant papers
# please see the original publications for these papers

# Suhre et al. 2017. Nat Comm. Supplementary Dataset 1. KORA
suhre_2017 <- read.csv("pqtl_data/suhre_2017.csv", header = T, stringsAsFactors = F)
suhre_2017_replicate <- subset(suhre_2017, suhre_2017$EntrezGeneSymbol %in% druggable_genome_sign$gene_display_label)
suhre_2017_replicate$gene_and_source <- str_c(suhre_2017_replicate$EntrezGeneSymbol, "_suhre2017")
suhre_2017_replicate$position <- gsub(",","",suhre_2017_replicate$Position)
suhre_2017_replicate$maf <- gsub("%", "", suhre_2017_replicate$maf)
suhre_2017_replicate$maf <- as.numeric(suhre_2017_replicate$maf)/100
suhre_2017_replicate_keep <- suhre_2017_replicate[,c("SNP", "Chr","position","Minor.Allele", "Major.Allele","maf","Beta.inv","S.E.","P.Value..inv.","EntrezGeneSymbol","gene_and_source","N")]
names(suhre_2017_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")

# Sun et al. 2018. Nature. Supplementary Table 4.
sun_2018 <- read.csv("pqtl_data/Sun_2018.csv", header = T, stringsAsFactors = F)
sun_2018_replicate <- subset(sun_2018, sun_2018$Mapped.gene %in% druggable_genome_sign$gene_display_label)
# sun_2018_replicate$Target[sun_2018_replicate$Target == "Cathepsin B"] <- "CTSB"
# sun_2018_replicate$Target[sun_2018_replicate$Target == "DHPR"] <- "QDPR"
sun_2018_replicate$gene_and_source <- str_c(sun_2018_replicate$Mapped.gene, "_sun2018")
sun_2018_replicate$position <- gsub(",","",sun_2018_replicate$position)
sun_2018_replicate$n <- 3301
sun_2018_replicate_keep <- sun_2018_replicate[,c("Sentinel.variant.", "Chr","position","Effect.Allele..EA.", "Other.Allele..OA.","EAF","meta_beta","meta_se","meta_p","Mapped.gene","gene_and_source","n")]
names(sun_2018_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")

#Pietzner et al. 2021. Fenland study
pietzner_2021 <- read.csv("pqtl_data/pietzner_2021.csv", header = T, stringsAsFactors = F)
pietzner_2021$HGNC.symbol.protein <- sub("^'", "", pietzner_2021$HGNC.symbol.protein)
pietzner_2021_replicate <- subset(pietzner_2021, pietzner_2021$HGNC.symbol.protein %in% druggable_genome_sign$gene_display_label)
pietzner_2021_replicate$gene_and_source <- str_c(pietzner_2021_replicate$HGNC.symbol.protein, "_pietzner2021")
pietzner_2021_replicate$n <- 12084
head(pietzner_2021_replicate)
head(sun_2018_replicate_keep) 
pietzner_2021_replicate_keep <- pietzner_2021_replicate[,c("rsID", "Chr","Position","EA", "NEA","EAF","Effect","SE","P.value","HGNC.symbol.protein","gene_and_source","n")]
names(pietzner_2021_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")

#Ferkingstad et al. 2021. NG deCODE study
ferkingstad_2021 <- read.csv("pqtl_data/Ferkingstad_2021.csv", header = T, stringsAsFactors = F)
ferkingstad_2021_replicate <- subset(ferkingstad_2021, ferkingstad_2021$gene...prot.. %in% druggable_genome_sign$gene_display_label)
ferkingstad_2021_replicate$gene_and_source <- str_c(ferkingstad_2021_replicate$gene...prot.., "_ferkingstad2021")
ferkingstad_2021_replicate$n <- 35559
ferkingstad_2021_replicate <- ferkingstad_2021_replicate %>%
  mutate(SE = sqrt(variance.expl..var.. / 35559))
ferkingstad_2021_replicate$chr..var.. <- sub("^chr", "", ferkingstad_2021_replicate$chr..var..)
ferkingstad_2021_replicate$MAF.... <- as.numeric(ferkingstad_2021_replicate$MAF....)/100 #这个要换算成EAF
ferkingstad_2021_replicate$P <- 10^ferkingstad_2021_replicate$X.Log10.P...adj..
ferkingstad_2021_replicate_keep <- ferkingstad_2021_replicate[,c("variant", "chr..var..","pos..var..","Amaj", "Amin","MAF....","beta..adj..","SE","P","gene...prot..","gene_and_source","n")]
names(ferkingstad_2021_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")

#Gudjonsson et al. 2022 NC islander cohort2
gudjonsson_2022 <- read.csv("pqtl_data/Gudjonsson_2022.csv", header = T, stringsAsFactors = F)

gudjonsson_2022_replicate <- subset(gudjonsson_2022, gudjonsson_2022$Protein..Entrez.symbol. %in% druggable_genome_sign$gene_display_label)
head(gudjonsson_2022_replicate)
gudjonsson_2022_replicate$gene_and_source <- str_c(gudjonsson_2022_replicate$Protein..Entrez.symbol., "_gudjonsson2022")
gudjonsson_2022_replicate$n <- 5368
library(data.table)
SNP <- fread('./hg19_avsnp150_cleaned.txt.gz')
gudjonsson_2022_replicate <- left_join(gudjonsson_2022_replicate, SNP, by = c("rsID" = "SNP"))
head(gudjonsson_2022_replicate)
gudjonsson_2022_replicate_keep <- gudjonsson_2022_replicate[,c("rsID", "Chr","Pos..GRCh37.","Allele1", "Allele2.y","EAF","beta","se","P","Protein..Entrez.symbol.","gene_and_source","n")]
names(gudjonsson_2022_replicate_keep) <- c("snp", "chr","pos","effect_allele", "other_allele","eaf","beta","se","pval","protein","gene_and_source","sample_size")

# combine
complete_pqtl_data_for_druggable_genome_replication <- suhre_2017_replicate_keep %>% rbind(sun_2018_replicate_keep) %>% rbind(pietzner_2021_replicate_keep) %>% rbind(ferkingstad_2021_replicate_keep) %>% rbind(gudjonsson_2022_replicate_keep) 
complete_pqtl_data_for_druggable_genome_replication <- distinct(complete_pqtl_data_for_druggable_genome_replication)
#complete_pqtl_data_for_druggable_genome_replication <- fread('./pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt')

# add gene info from druggable genome file
gene_chr_pos <- druggable_genome_sign[,c("gene_display_label","chr_name","gene_start","gene_end")]
names(gene_chr_pos) <- c("protein","gene_chr","gene_start","gene_end")

complete_pqtl_data_for_druggable_genome_replication1 <- left_join(complete_pqtl_data_for_druggable_genome_replication, gene_chr_pos, by = "protein")

# clarify if pQTL is acting in cis or trans
complete_pqtl_data_for_druggable_genome_replication1$cis_trans <- NA
complete_pqtl_data_for_druggable_genome_replication1$cis_trans[complete_pqtl_data_for_druggable_genome_replication1$chr == complete_pqtl_data_for_druggable_genome_replication1$gene_chr] <- "cis"
complete_pqtl_data_for_druggable_genome_replication1$cis_trans[complete_pqtl_data_for_druggable_genome_replication1$chr != complete_pqtl_data_for_druggable_genome_replication1$gene_chr] <- "trans"

# check how many proteins were covered
length(unique(complete_pqtl_data_for_druggable_genome_replication1$gene_and_source))

write.table(complete_pqtl_data_for_druggable_genome_replication, "pqtl_data/complete_pqtl_data_for_druggable_genome_replication.txt",sep ="\t",row.names = F)
write.table(complete_pqtl_data_for_druggable_genome_replication1, "pqtl_data/complete_pqtl_data_for_druggable_genome_replication1.txt",sep ="\t",row.names = F)
