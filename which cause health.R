#什么疾病影响寿命？

rm(list=ls())
setwd("C:/Users/hou/Desktop")
gc()
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(beepr)

path_exp <- "C:/Users/hou/Desktop/exp/"

exposure_id <- fread("exp.txt",header = T)

for (i in 1:nrow(exposure_id)) {   #length(exposure_id)
  tryCatch( 
    {
      gc()
      exposure_data <- data.frame()
      print(paste0("这是第",i,"个表型",as.character(exposure_id[i])))
      exposure_data <- extract_instruments(as.character(exposure_id[i]),r2 = 0.01)
      fwrite(exposure_data,paste0(path_exp,as.character(exposure_id[i]),".csv"))
    },error=function(e)
    {
      print(paste0("第",i,"个表型",as.character(exposure_id[i]),"无法做MR"))
    }
  )
  
}
beep(4)





rm(list=ls())
setwd("C:/Users/hou/Desktop")
gc()
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(beepr)

exposure_id <- dir("C:/Users/hou/Desktop/exp/")
path_res <- "C:/Users/hou/Desktop/long/"


for (i in 1:length(exposure_id)) {   #length(exposure_id)
  tryCatch( 
    {
      gc()
      exposure_data <- data.frame()
      exposure_data <- fread(paste0("C:/Users/hou/Desktop/exp/",exposure_id[i]) ,header = T)
      outcome_data <- data.frame()
      H_data <- data.frame()
      mr_results <- data.frame()
      print(paste0("这是第",i,"个表型",exposure_id[i]))
      outcome_data <- extract_outcome_data(
        snps = exposure_data$SNP, 
        outcomes = c("ebi-a-GCST006702"))
      H_data <- harmonise_data(
        exposure_dat = exposure_data, 
        outcome_dat = outcome_data
      )
      mr_results<-mr(H_data)
      mr_results
      print(mr_results)
      gc()
      subset_on_method(mr_results)
      fwrite(mr_results,paste0(path_res,exposure_id[i]))
    },error=function(e)
    {
      print(paste0("第",i,"个表型",exposure_id[i],"无法做MR"))
    }
  )
  
}
beep(4)


setwd("C:/Users/hou/Desktop/long/")
filename <- dir()
zhenghe <- data.frame()
for(i in 1:length(filename))
{
  temp <- fread(filename[i],header = T)
  zhenghe <- rbind(zhenghe,temp)
}
fwrite(zhenghe,"C:/Users/hou/Desktop/who_cause_health.csv")
gc()

