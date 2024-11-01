###长寿 整理ENSG暴露：5*10-8,0.01
gc()
setwd("C:/Users/hou/Desktop/")
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(ieugwasr)

ensg <- fread("ensg.txt",header = T)
exppath <- "C:/Users/hou/Desktop/ensgexp/"
gc()

temp <- 1

i <- temp
for(i in 1:nrow(ensg))  #nrow(ensg)
{
  exposure_data <- data.frame()
  exposure_data <- extract_instruments(
    outcomes = as.character(ensg[i]),r2 = 0.01)
  print(paste0("第",i,"个基因，是",ensg[i]))
    if(!is.null(exposure_data))
    {
      fwrite(exposure_data,paste0(exppath,as.character(ensg[i]),".csv"))
      print(paste0("第",i,"个基因",ensg[i],"找到位点了"))
    }
}
print(i)
temp <- i




