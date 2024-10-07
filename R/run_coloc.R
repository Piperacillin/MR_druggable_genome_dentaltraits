
library("coloc")
library("stringr")
getwd()
rm(list=ls())
#GENE <- Sys.Getenv("GENE")
#EXPOSURE_DATA <- Sys.Getenv("EXPOSURE_DATA")
#OUTCOME <- Sys.Getenv("OUTCOME")


# set priors
p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5

# read in data
dat_harmonized4coloc_files <- list.files(path = "./coloc", pattern = "dat_harmonized4coloc", full.names = TRUE)
dat_harmonized4coloc_files <- subset(dat_harmonized4coloc_files, !grepl("_full.txt$", dat_harmonized4coloc_files))

coloc_results_gene <- list()

for (i in 1:length(dat_harmonized4coloc_files)){

    dat <- read.table(dat_harmonized4coloc_files[i],header=T,sep="\t",stringsAsFactors=F)
    
    dat <- dat[complete.cases(dat$eaf.outcome), ]
    
    GENE <- dat$exposure[1]
    OUTCOME <- dat$outcome[1]
    EXPOSURE_DATA <- dat$tissue[1]
    
    
    dat$var.exposure <- dat$se.exposure^2
    dat$var.outcome <- dat$se.outcome^2


    if (OUTCOME == "EUR_bleedinggums_prep") {

        PD_cc_ratio <- 60210/(60210+400821)

        # using eaf as maf

        coloc_results_prior <-
      coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                 varbeta = dat$var.exposure,  # standard error squared
                                 N = dat$samplesize.exposure,  # sample size # might work without it
                                 type = "quant",
                                 sdY = 1,
                                 #MAF = dat$eaf.exposure,
                                 snp = dat$SNP),
                dataset2 = list( beta = dat$beta.outcome,
                                 varbeta = dat$var.outcome,
                                 type = "cc", # case control
                                 N = dat$samplesize.outcome,
                                 s = PD_cc_ratio,
                                 MAF = dat$eaf.outcome,
                                 snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
                             )
    } else if (OUTCOME == "EUR_dentures_prep") {

        PD_cc_ratio <- 77714/(77714+383317)

        # using eaf as maf
        coloc_results_prior <-
        coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                 varbeta = dat$var.exposure,  # standard error squared
                                 N = dat$samplesize.exposure,  # sample size # might work without it
                                 type = "quant",
                                 sdY = 1,
                                 #MAF = dat$eaf.exposure,
                                 snp = dat$SNP),
                dataset2 = list( beta = dat$beta.outcome,
                                 varbeta = dat$var.outcome,
                                 type = "cc", # case control
                                 N = dat$samplesize.outcome,
                                 s = PD_cc_ratio,
                                 MAF = dat$eaf.outcome,
                                 snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
                             )

    } else if (OUTCOME == "EUR_looseteeth_prep") {
      
      PD_cc_ratio <- 18979/(18979+442052)
      
      # using eaf as maf
      coloc_results_prior <-
        coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                   varbeta = dat$var.exposure,  # standard error squared
                                   N = dat$samplesize.exposure,  # sample size # might work without it
                                   type = "quant",
                                   sdY = 1,
                                   #MAF = dat$eaf.exposure,
                                   snp = dat$SNP),
                  dataset2 = list( beta = dat$beta.outcome,
                                   varbeta = dat$var.outcome,
                                   type = "cc", # case control
                                   N = dat$samplesize.outcome,
                                   s = PD_cc_ratio,
                                   MAF = dat$eaf.outcome,
                                   snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
        )
      
    } else if (OUTCOME == "EUR_ulcers_prep") {
      
      PD_cc_ratio <- 47091/(47091+413940)
      
      # using eaf as maf
      coloc_results_prior <-
        coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                   varbeta = dat$var.exposure,  # standard error squared
                                   N = dat$samplesize.exposure,  # sample size # might work without it
                                   type = "quant",
                                   sdY = 1,
                                   #MAF = dat$eaf.exposure,
                                   snp = dat$SNP),
                  dataset2 = list( beta = dat$beta.outcome,
                                   varbeta = dat$var.outcome,
                                   type = "cc", # case control
                                   N = dat$samplesize.outcome,
                                   s = PD_cc_ratio,
                                   MAF = dat$eaf.outcome,
                                   snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
        )
      
    } else if (OUTCOME == "finngen_R9_K11_periodontitis_prep") {
      
      PD_cc_ratio <- 30377/(30377+346900)
      
      # using eaf as maf
      coloc_results_prior <-
        coloc.abf(dataset1 = list( beta = dat$beta.exposure, # beta
                                   varbeta = dat$var.exposure,  # standard error squared
                                   N = dat$samplesize.exposure,  # sample size # might work without it
                                   type = "quant",
                                   sdY = 1,
                                   #MAF = dat$eaf.exposure,
                                   snp = dat$SNP),
                  dataset2 = list( beta = dat$beta.outcome,
                                   varbeta = dat$var.outcome,
                                   type = "cc", # case control
                                   N = dat$samplesize.outcome,
                                   s = PD_cc_ratio,
                                   MAF = dat$eaf.outcome,
                                   snp = dat$SNP), p1 = p1, p2 = p2, p12 = p12
        )
      
    }

    
    # make a pretty table

    coloc_results_prior[["summary"]]["p1"] <- coloc_results_prior[["priors"]]["p1"]
    coloc_results_prior[["summary"]]["p2"] <- coloc_results_prior[["priors"]]["p2"]
    coloc_results_prior[["summary"]]["p12"] <- coloc_results_prior[["priors"]]["p12"]

    coloc_results_prior[["summary"]]["gene"] <- GENE
    coloc_results_prior[["summary"]]["outcome"] <- OUTCOME
    coloc_results_prior[["summary"]]["exposure_data"] <- EXPOSURE_DATA

    if(i == 1){

      coloc_results_gene_summary <- coloc_results_prior[["summary"]] %>% t() %>% tibble::as_tibble()

    } else {

      coloc_results_gene_summary <-
        dplyr::bind_rows(coloc_results_gene_summary,
                  coloc_results_prior[["summary"]] %>% t() %>% tibble::as_tibble())

    }
}


coloc_results_gene_summary <- subset(coloc_results_gene_summary, !duplicated(coloc_results_gene_summary)) #去除重复行

# keep only powered analyses
coloc_results_gene_summary$sum_pph3_pph4 <- as.numeric(coloc_results_gene_summary$PP.H3.abf) + as.numeric(coloc_results_gene_summary$PP.H4.abf)

coloc_results_gene_summary_powered <- coloc_results_gene_summary[which(coloc_results_gene_summary$sum_pph3_pph4 > 0.8),]

write.table(coloc_results_gene_summary_powered, str_c("coloc/res_coloc_all_outcomes_p1_",p1,"_p2_",p2,"_p12_",p12, ".txt"),sep="\t",quote=F,row.names=F)

print("mission_complete")
