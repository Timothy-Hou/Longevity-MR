###decode附录文件转
gc()
setwd("C:/Users/hou/Desktop/")
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(ieugwasr)
library(readxl)
library(stringr)

pqtl <- fread("pqtl.csv",header = T)
pqtl$se <- abs(pqtl$beta/qnorm(pqtl$p/2))
pqtl$MAF <- pqtl$MAF/100

head(pqtl)
pqtlexp <- format_data(
  pqtl,
  type = "exposure",
  phenotype_col = "prot",
  gene_col = "protlong",
  snps = NULL,
  header = TRUE,
  snp_col = "variant",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "MAF",
  effect_allele_col = "Amin",
  other_allele_col = "Amaj",
  pval_col = "p",
  samplesize_col = 35599,
  info_col = "Info"
)
head(pqtlexp)
fwrite(pqtlexp,"pqtlexp.csv")  ##把gene_exposure所有/<>|?:*之类的都换成_
pqtlexp <- fread("pqtlexp.csv",header = T)

a <- split.data.frame(
  x = pqtlexp,
  f = as.factor(pqtlexp$gene_exposure))
exppath <- "C:/Users/hou/Desktop/pqtl/"
for (i in 1:length(a)) ###输出文件
{
  path <- paste0(exppath,a[[i]][["gene_exposure"]][1],".csv")
  fwrite(a[[i]],path)
}




