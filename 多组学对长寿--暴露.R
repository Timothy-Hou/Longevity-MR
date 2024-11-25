###长寿 整理ENSG暴露：5*10-8,0.01
rm(list=ls())
setwd("/Users/hou/Desktop")
gc()
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(ieugwasr)
library(beepr)

ensg <- fread("ensg.txt",header = T)
exppath <- "/Users/hou/Desktop/ensgexp/"
gc()

temp <- 1
errorvalue <- c(28,176,341,453,499,525,619,755,1191)
for(i in temp:nrow(ensg))  #nrow(ensg)
{
  exposure_data <- data.frame()
  print(paste0("第",i,"个基因，是",ensg[i]))
  tryCatch({
  exposure_data <- extract_instruments(
    outcomes = as.character(ensg[i]),r2 = 0.01)},
  error = function(e){
    print(paste0("第",i,"个基因",ensg[i],"遭遇网络故障～～"))
    beep(6)
    errorvalue[length(errorvalue)+1] <- i  ##函数内的变化无效
    print(errorvalue)
    })
  if(!is.null(exposure_data))
  {
    fwrite(exposure_data,paste0(exppath,as.character(ensg[i]),".csv"))
    print(paste0("第",i,"个基因",ensg[i],"找到",nrow(exposure_data),"个位点了"))
  }
  else(print(paste0("第",i,"个基因",ensg[i],"没能够找到位点！！")))
}
print(i)
temp <- i  ## 无法做的：28，176,341,453,499,525,619,755,1191,1235  第13898个基因eqtl-a-ENSG00000241493找到1个位点了"
gc()



temp <- 19907
#######失败的部分重做一遍

setwd("/Users/hou/Desktop/empty")
ensg <- gsub(".csv","",dir())
exppath <- "/Users/hou/Desktop/ensgexp/"
for(i in 1:length(ensg))  #nrow(ensg)
{
  exposure_data <- data.frame()
  print(paste0("第",i,"个基因，是",ensg[i]))
  tryCatch({
    exposure_data <- extract_instruments(
      outcomes = ensg[i],r2 = 0.01)},
    error = function(e){
      print(paste0("第",i,"个基因",ensg[i],"遭遇网络故障～～"))
      beep(6)
      print(errorvalue)
    })
  if(!is.null(exposure_data))
  {
    fwrite(exposure_data,paste0(exppath,as.character(ensg[i]),".csv"))
    print(paste0("第",i,"个基因",ensg[i],"找到",nrow(exposure_data),"个位点了"))
  }
  else(print(paste0("第",i,"个基因",ensg[i],"没能够找到位点！！")))
}
print(i)
gc()
########仅提取到1个基因



######用ieu的数据库做一下
rm(list=ls())
setwd("/Users/hou/Desktop/ensgexp/")
gc()
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(ieugwasr)
library(beepr)
filename <- dir()
ensg <- gsub(".csv","",filename)
orpath <- "/Users/hou/Desktop/resor/"##注意末尾需要有/
plepath <- "/Users/hou/Desktop/resple/"
hpath <- "/Users/hou/Desktop/h/"
hetpath <- "/Users/hou/Desktop/het/"

weberror <- "/Users/hou/Desktop/weberror/"

errorval <- c("无法运算的包括")
temp <- 1   ##上次error到第几个了？


for(k in temp:length(filename))##此处是k，不能是i  
{
  print(paste0("现在正在做第",k,"个基因，是",ensg[k]))
  exposure_gwas <- data.frame()
  exposure_data <- data.frame()
  H_data <- data.frame() ##新建h和out，以免出错
  outcome_data <- data.frame()
  exposure_gwas <- fread(filename[k],header = T)
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
    samplesize_col = "samplesize",
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
    outcome_data <- extract_outcome_data(
      snps = exposure_data$SNP, 
      outcomes = c("ebi-a-GCST003395", "ebi-a-GCST006702","ebi-a-GCST003394","ebi-a-GCST006702","ukb-b-3018")
  )},
  error = function(e){
    print(paste0("第",k,"个基因",ensg[k],"遭遇网络故障～～"))
    beep(6)
    fwrite(exposure_data,paste0(weberror,ensg[k]))
  })
  
  if(!is.null(outcome_data))   ###条件判断1:如果提取的不是空集（即找到了对应的snp）
  {
    if(nrow(outcome_data)!=0)   ##条件判断2:如果nrow确实找到了，且不是错误找到（里面真的有outcome）
      {
        H_data <- harmonise_data(
        exposure_dat = exposure_data, 
        outcome_dat = outcome_data
        )
  if(any(H_data$mr_keep==T))  ##条件判断3:如果所有位点都因为eaf约等于0.5而被删掉，则进行mr，若执行会报错
  {
    H_data$id.exposure <- ensg[k]##为了方便以后整合数据，把此处的名字写上
    H_data$exposure <- ensg[k]##为了方便以后整合数据，把此处的名字写上
    mr_results<-mr(H_data)
    res <- generate_odds_ratios(mr_results)
    pleio <- mr_pleiotropy_test(H_data)
    het <- mr_heterogeneity(H_data)
    fwrite(res, paste(orpath,ensg[k],"-OR.csv",sep = ""))
    fwrite(pleio, paste(plepath,ensg[k],"-pleio.csv",sep = ""))
    fwrite(H_data,paste(hpath,ensg[k],"-h_data.csv",sep = ""))
    fwrite(het,paste(hetpath,ensg[k],"-het.csv",sep = ""))
    gc()
  }else{errorval <- append(errorval,ensg[k])   ##条件判断3：这个循环表示：如果所有的h之后的都被删了，则输出这个无法分析的基因
        print(errorval)
        fwrite(exposure_data,paste0(weberror,ensg[k]))
        beep(6)}
  }else{errorval <- append(errorval,ensg[k])   ##条件判断2：这个循环表示：如果是outcome是零行（而非空集），则输出该基因
        print(errorval)
        fwrite(exposure_data,paste0(weberror,ensg[k]))
        beep(6)}
  }else{errorval <- append(errorval,ensg[k])   ##条件判断1：这个循环表示：如果是空集，则输出该基因
  print(errorval)
  fwrite(exposure_data,paste0(weberror,ensg[k]))
  beep(6)}
  
}

print(errorval)
beep(4)


##重算weberror的
rm(list=ls())
setwd("/Users/hou/Desktop/weberror/")
gc()
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(ieugwasr)
library(beepr)
filename <- dir()
ensg <- gsub(".csv","",filename)
orpath <- "/Users/hou/Desktop/resor/"##注意末尾需要有/
plepath <- "/Users/hou/Desktop/resple/"
hpath <- "/Users/hou/Desktop/h/"
hetpath <- "/Users/hou/Desktop/het/"

weberror <- "/Users/hou/Desktop/weberror2/"

errorval <- c("无法运算的包括")
temp <- 1   ##上次error到第几个了？


for(k in temp:length(filename))##此处是k，不能是i  
{
  print(paste0("现在正在做第",k,"个基因，是",ensg[k]))
  exposure_gwas <- data.frame()
  exposure_data <- data.frame()
  H_data <- data.frame() ##新建h和out，以免出错
  outcome_data <- data.frame()
  exposure_gwas <- fread(filename[k],header = T)
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
    samplesize_col = "samplesize",
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
    outcome_data <- extract_outcome_data(
      snps = exposure_data$SNP, 
      outcomes = c("ebi-a-GCST003395", "ebi-a-GCST006702","ebi-a-GCST003394","ebi-a-GCST006702","ukb-b-3018")
    )},
    error = function(e){
      print(paste0("第",k,"个基因",ensg[k],"遭遇网络故障～～"))
      beep(6)
      fwrite(exposure_data,paste0(weberror,ensg[k]))
    })
  
  if(!is.null(outcome_data))   ###条件判断1:如果提取的不是空集（即找到了对应的snp）
  {
    if(nrow(outcome_data)!=0)   ##条件判断2:如果nrow确实找到了，且不是错误找到（里面真的有outcome）
    {
      H_data <- harmonise_data(
        exposure_dat = exposure_data, 
        outcome_dat = outcome_data
      )
      if(any(H_data$mr_keep==T))  ##条件判断3:如果所有位点都因为eaf约等于0.5而被删掉，则进行mr，若执行会报错
      {
        H_data$id.exposure <- ensg[k]##为了方便以后整合数据，把此处的名字写上
        H_data$exposure <- ensg[k]##为了方便以后整合数据，把此处的名字写上
        mr_results<-mr(H_data)
        res <- generate_odds_ratios(mr_results)
        pleio <- mr_pleiotropy_test(H_data)
        het <- mr_heterogeneity(H_data)
        fwrite(res, paste(orpath,ensg[k],"-OR.csv",sep = ""))
        fwrite(pleio, paste(plepath,ensg[k],"-pleio.csv",sep = ""))
        fwrite(H_data,paste(hpath,ensg[k],"-h_data.csv",sep = ""))
        fwrite(het,paste(hetpath,ensg[k],"-het.csv",sep = ""))
        gc()
      }else{errorval <- append(errorval,ensg[k])   ##条件判断3：这个循环表示：如果所有的h之后的都被删了，则输出这个无法分析的基因
      print(errorval)
      fwrite(exposure_data,paste0(weberror,ensg[k]))
      beep(6)}
    }else{errorval <- append(errorval,ensg[k])   ##条件判断2：这个循环表示：如果是outcome是零行（而非空集），则输出该基因
    print(errorval)
    fwrite(exposure_data,paste0(weberror,ensg[k]))
    beep(6)}
  }else{errorval <- append(errorval,ensg[k])   ##条件判断1：这个循环表示：如果是空集，则输出该基因
  print(errorval)
  fwrite(exposure_data,paste0(weberror,ensg[k]))
  beep(6)}
  
}

print(errorval)
beep(4)



