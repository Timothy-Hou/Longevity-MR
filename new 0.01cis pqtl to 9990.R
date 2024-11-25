#0.01 的蛋白组，对99和90做分析

#  NC99
rm(list=ls())
setwd("/Users/hou/Desktop")
gc()
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(ieugwasr)
library(beepr)


dir.create("resor")
dir.create("resple")
dir.create("h")
dir.create("het")
dir.create("notfound")

outcome_gwas <- fread("Results_99th_percentile_processed.txt",header = T) %>% as.data.frame()
setwd("/Users/hou/Desktop/gtex1e5_ld0.01/")    #pqtlexp    #newcispqtl0.01  #gtex1e5_ld0.01  #cis_decodepqtl
filename <- dir()
pqtl <- gsub(".csv","",filename)
orpath <- "/Users/hou/Desktop/resor/"##注意末尾需要有/
plepath <- "/Users/hou/Desktop/resple/"
hpath <- "/Users/hou/Desktop/h/"
hetpath <- "/Users/hou/Desktop/het/"

notfound <- "/Users/hou/Desktop/notfound/"  ##此处不再是weberror

errorval <- c("无法运算的包括")
misssnp <- data.frame()
temp <- 1   ##上次error到第几个了？
gc()
for(k in temp:length(filename))##此处是k，不能是i  length(filename)
{
  if(is.integer(k/20)){gc()}
  print(paste0("现在正在做第",k,"个蛋白，是",pqtl[k]))
  exposure_gwas <- data.frame()
  exposure_data <- data.frame()
  H_data <- data.frame() ##新建h和out，以免出错
  outcome_data <- data.frame()
  exposure_gwas <- fread(filename[k],header = T) %>% as.data.frame()
  exposure_data <- format_data(
    exposure_gwas,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    units_col = "units.exposure",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    samplesize_col = "samplesize.exposure",
    gene_col = "gene.exposure",
    id_col = "id",
    min_pval = 1e-200,
    z_col = "z",
    info_col = "info",
    chr_col = "chr.exposure",
    pos_col = "pos",
    log_pval = FALSE
  )
  tryCatch({   #尝试提取结局，若因网络中断，则停止提取，输出暴露
    outcome_data <- format_data(
      dat = outcome_gwas,type = "outcome",snps = exposure_data$SNP,header = T,
      snp_col = "name",
      beta_col = "Beta",
      se_col = "SE",
      eaf_col = "EAF",
      effect_allele_col = "EA",other_allele_col = "NEA",
      pval_col = "P-value",
      samplesize_col = "Effective_N",
      chr_col = "Chr",pos_col = "Position")},
    error = function(e){
      print(paste0("第",k,"个蛋白",pqtl[k],"无法提取有效暴露"))
      
      fwrite(exposure_data,paste0(notfound,pqtl[k]))
    })
  if(!is.null(outcome_data))   ###条件判断1:如果提取的不是空集（即找到了对应的snp）
  {
    if(nrow(outcome_data)!=0)   ##条件判断2:如果nrow确实找到了，且不是错误找到（里面真的有outcome）
    {
      row <- c(pqtl[k],nrow(exposure_data),nrow(outcome_data),(nrow(exposure_data)-nrow(outcome_data)))   #看看每一个蛋白分别少了几个snp？
      print(row)
      misssnp <- rbind(misssnp,row)
      H_data <- harmonise_data(
        exposure_dat = exposure_data, 
        outcome_dat = outcome_data,action = 1
      )
      if(any(H_data$mr_keep==T))  ##条件判断3:如果所有位点都因为eaf约等于0.5而被删掉，则进行mr，若执行会报错
      {
        H_data$id.exposure <- pqtl[k]##为了方便以后整合数据，把此处的名字写上
        H_data$exposure <- pqtl[k]##为了方便以后整合数据，把此处的名字写上
        H_data$outcome <- "nc-l99"
        H_data$id.outcome <- "nc-l99"
        mr_results<-mr(H_data)
        res <- generate_odds_ratios(mr_results)
        pleio <- mr_pleiotropy_test(H_data)
        het <- mr_heterogeneity(H_data)
        fwrite(res, paste(orpath,pqtl[k],"-OR.csv",sep = ""))
        fwrite(pleio, paste(plepath,pqtl[k],"-pleio.csv",sep = ""))
        fwrite(H_data,paste(hpath,pqtl[k],"-h_data.csv",sep = ""))
        fwrite(het,paste(hetpath,pqtl[k],"-het.csv",sep = ""))
        gc()
      }else{errorval <- append(errorval,pqtl[k])   ##条件判断3：这个循环表示：如果所有的h之后的都被删了，则输出这个无法分析的蛋白
      print(errorval)
      fwrite(exposure_data,paste0(notfound,pqtl[k]))
      }
    }else{errorval <- append(errorval,pqtl[k])   ##条件判断2：这个循环表示：如果是outcome是零行（而非空集），则输出该蛋白
    print(errorval)
    fwrite(exposure_data,paste0(notfound,pqtl[k]))
    }
  }else{errorval <- append(errorval,pqtl[k])   ##条件判断1：这个循环表示：如果是空集，则输出该蛋白
  print(errorval)
  fwrite(exposure_data,paste0(notfound,pqtl[k]))
  }
}

print(errorval)
fwrite(misssnp,"/Users/hou/Desktop/misssnp.csv")
write.table(errorval,"/Users/hou/Desktop/error.txt")
beep(4)


#####  NC90
rm(list=ls())
setwd("/Users/hou/Desktop")
gc()
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(ieugwasr)
library(beepr)


dir.create("resor")
dir.create("resple")
dir.create("h")
dir.create("het")
dir.create("notfound")

outcome_gwas <- fread("Results_90th_percentile_processed.txt",header = T) %>% as.data.frame()
setwd("/Users/hou/Desktop/gtex1e5_ld0.01/")    #pqtlexp    #newcispqtl0.01      cis_decodepqtl
filename <- dir()
pqtl <- gsub(".csv","",filename)
orpath <- "/Users/hou/Desktop/resor/"##注意末尾需要有/
plepath <- "/Users/hou/Desktop/resple/"
hpath <- "/Users/hou/Desktop/h/"
hetpath <- "/Users/hou/Desktop/het/"

notfound <- "/Users/hou/Desktop/notfound/"  ##此处不再是weberror

errorval <- c("无法运算的包括")
misssnp <- data.frame()
temp <- 1   ##上次error到第几个了？
gc()
for(k in temp:length(filename))##此处是k，不能是i  length(filename)
{
  if(is.integer(k/20)){gc()}
  print(paste0("现在正在做第",k,"个蛋白，是",pqtl[k]))
  exposure_gwas <- data.frame()
  exposure_data <- data.frame()
  H_data <- data.frame() ##新建h和out，以免出错
  outcome_data <- data.frame()
  exposure_gwas <- fread(filename[k],header = T) %>% as.data.frame()
  exposure_data <- format_data(
    exposure_gwas,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    units_col = "units.exposure",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    samplesize_col = "samplesize.exposure",
    gene_col = "gene.exposure",
    id_col = "id",
    min_pval = 1e-200,
    z_col = "z",
    info_col = "info",
    chr_col = "chr.exposure",
    pos_col = "pos",
    log_pval = FALSE
  )
  tryCatch({   #尝试提取结局，若因网络中断，则停止提取，输出暴露
    outcome_data <- format_data(
      dat = outcome_gwas,type = "outcome",snps = exposure_data$SNP,header = T,
      snp_col = "name",
      beta_col = "Beta",
      se_col = "SE",
      eaf_col = "EAF",
      effect_allele_col = "EA",other_allele_col = "NEA",
      pval_col = "P-value",
      samplesize_col = "Effective_N",
      chr_col = "Chr",pos_col = "Position")},
    error = function(e){
      print(paste0("第",k,"个蛋白",pqtl[k],"无法提取有效暴露"))
      
      fwrite(exposure_data,paste0(notfound,pqtl[k]))
    })
  if(!is.null(outcome_data))   ###条件判断1:如果提取的不是空集（即找到了对应的snp）
  {
    if(nrow(outcome_data)!=0)   ##条件判断2:如果nrow确实找到了，且不是错误找到（里面真的有outcome）
    {
      row <- c(pqtl[k],nrow(exposure_data),nrow(outcome_data),(nrow(exposure_data)-nrow(outcome_data)))   #看看每一个蛋白分别少了几个snp？
      print(row)
      misssnp <- rbind(misssnp,row)
      H_data <- harmonise_data(
        exposure_dat = exposure_data, 
        outcome_dat = outcome_data,action = 1
      )
      if(any(H_data$mr_keep==T))  ##条件判断3:如果所有位点都因为eaf约等于0.5而被删掉，则进行mr，若执行会报错
      {
        H_data$id.exposure <- pqtl[k]##为了方便以后整合数据，把此处的名字写上
        H_data$exposure <- pqtl[k]##为了方便以后整合数据，把此处的名字写上
        H_data$outcome <- "nc-l90"
        H_data$id.outcome <- "nc-l90"
        mr_results<-mr(H_data)
        res <- generate_odds_ratios(mr_results)
        pleio <- mr_pleiotropy_test(H_data)
        het <- mr_heterogeneity(H_data)
        fwrite(res, paste(orpath,pqtl[k],"-OR.csv",sep = ""))
        fwrite(pleio, paste(plepath,pqtl[k],"-pleio.csv",sep = ""))
        fwrite(H_data,paste(hpath,pqtl[k],"-h_data.csv",sep = ""))
        fwrite(het,paste(hetpath,pqtl[k],"-het.csv",sep = ""))
        gc()
      }else{errorval <- append(errorval,pqtl[k])   ##条件判断3：这个循环表示：如果所有的h之后的都被删了，则输出这个无法分析的蛋白
      print(errorval)
      fwrite(exposure_data,paste0(notfound,pqtl[k]))
      }
    }else{errorval <- append(errorval,pqtl[k])   ##条件判断2：这个循环表示：如果是outcome是零行（而非空集），则输出该蛋白
    print(errorval)
    fwrite(exposure_data,paste0(notfound,pqtl[k]))
    }
  }else{errorval <- append(errorval,pqtl[k])   ##条件判断1：这个循环表示：如果是空集，则输出该蛋白
  print(errorval)
  fwrite(exposure_data,paste0(notfound,pqtl[k]))
  }
}

print(errorval)
fwrite(misssnp,"/Users/hou/Desktop/misssnp.csv")
write.table(errorval,"/Users/hou/Desktop/error.txt")
beep(4)


##为什么总是消失的批量or
library(data.table)


setwd("/Users/hou/Desktop/resor")
filename <- dir()
zhenghe <- data.frame()
for(i in 1:length(filename))
{
  temp <- fread(filename[i],header = T)
  zhenghe <- rbind(zhenghe,temp)
}
fwrite(zhenghe,"/Users/hou/Desktop/or.csv")
gc()


setwd("/Users/hou/Desktop/resple")
filename <- dir()
zhenghe <- data.frame()
for(i in 1:length(filename))
{
  temp <- fread(filename[i],header = T)
  zhenghe <- rbind(zhenghe,temp)
}
fwrite(zhenghe,"/Users/hou/Desktop/ple.csv")



setwd("/Users/hou/Desktop/het")
filename <- dir()
zhenghe <- data.frame()
for(i in 1:length(filename))
{
  temp <- fread(filename[i],header = T)
  zhenghe <- rbind(zhenghe,temp)
}
fwrite(zhenghe,"/Users/hou/Desktop/het.csv")

