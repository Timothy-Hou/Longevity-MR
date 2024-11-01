#####   多组学对长寿
setwd("C:/Users/hou/Desktop")
gc()

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)

##1.阳性对照
exposure_data <- extract_instruments(outcomes='ieu-a-798',r2 = 0.01) ##心梗作为暴露
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, 
  outcomes = c("ieu-a-1091","ebi-a-GCST006702","ebi-a-GCST003394","ieu-a-1094","ebi-a-GCST003395"))
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data)
subset_on_method(mr_results)
#2019 elife百万人：
outcome_data <- fread("2019-elife.tsv",header = T)
head(outcome_data)
outcome_data <- format_data(
  dat = outcome_data,
  type = "outcome",
  snps = exposure_data$SNP,
  snp_col = "rsid",
  beta_col = "beta1",
  se_col = "se",
  eaf_col = "freq1",
  effect_allele_col = "a1",
  other_allele_col = "a0",
  pval_col = "p",
  samplesize_col = "n",
  info_col = "info"
)
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data)
subset_on_method(mr_results)
gc()
#2019 nc l90：
outcome_data <- fread("l90.csv",header = T)
head(outcome_data)
outcome_data <- format_data(
  dat = outcome_data,
  type = "outcome",
  snps = exposure_data$SNP,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "FRQ",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  pval_col = "P",
  samplesize_col = "N",
)
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data)
subset_on_method(mr_results)
#2019 nc l99：
outcome_data <- fread("l99.csv",header = T)
head(outcome_data)
outcome_data <- format_data(
  dat = outcome_data,
  type = "outcome",
  snps = exposure_data$SNP,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "FRQ",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  pval_col = "P",
  samplesize_col = "N",
)
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data)
subset_on_method(mr_results)

################
##分别提取暴露并写出

tryCatch用法：
https://cloud.tencent.com/developer/article/1701469






