library("tidyverse")
library("data.table")
library(ggplot2)
library("Matrix")
library(Seurat)
library(viridis)
library(RColorBrewer)

# 加载需要的包
library(readr)
library(dplyr)

# 列出"data"文件夹中所有以".csv.gz"结尾的文件名
file_paths <- list.files(path = "./0_data/OAS_DB_COVID/OAS_DB_COVID_1",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)

# 使用 skip = 1 参数跳过每个文件的第一行
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))

# 将列表中的所有数据框合并为一个数据框
df_num <- bind_rows(df_list)

nrow(df_num)
colnames(df_num)

# 加载需要的包
library(readr)
library(dplyr)

# 列出"data"文件夹中所有以".csv.gz"结尾的文件名
file_paths <- list.files(path = "./0_data/OAS_DB_COVID/OAS_DB_COVID_ScoV1",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)

# 使用 skip = 1 参数跳过每个文件的第一行
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))

# 将列表中的所有数据框合并为一个数据框
df_ScoV1 <- bind_rows(df_list)

nrow(df_ScoV1)
colnames(df_ScoV1)

# 加载需要的包
library(readr)
library(dplyr)

# 列出"data"文件夹中所有以".csv.gz"结尾的文件名
file_paths <- list.files(path = "./0_data/OAS_DB_COVID/OAS_DB_COVID_SRR",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)

# 使用 skip = 1 参数跳过每个文件的第一行
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))

# 将列表中的所有数据框合并为一个数据框
df_SRR <- bind_rows(df_list)

nrow(df_SRR)
colnames(df_SRR)

nrow(df_num)
nrow(df_SRR)
nrow(df_ScoV1)

#验证列名格式是否一致
length(colnames(df_num))
length(colnames(df_SRR))
length(colnames(df_ScoV1))
sum(colnames(df_num) %in% colnames(df_SRR))
sum(colnames(df_num) %in% colnames(df_ScoV1))

BCR_COVID_data <- bind_rows(df_num,df_SRR,df_ScoV1)

nrow(BCR_COVID_data)

BCR_COVID_data$Virus_name <- "COVID"

write.csv(BCR_COVID_data,"2_BCR_COVID_data.csv")

# 加载需要的包
library(readr)
library(dplyr)

file_paths <- list.files(path = "./0_data/OAS_DB_other_virus/CMV",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))
df_CMV <- bind_rows(df_list)

              
file_paths <- list.files(path = "./0_data/OAS_DB_other_virus/HIV",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))
df_HIV <- bind_rows(df_list)
              
file_paths <- list.files(path = "./0_data/OAS_DB_other_virus/Obstructive_Sleep_Apnea",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))
df_Obstructive_Sleep_Apnea <- bind_rows(df_list)

file_paths <- list.files(path = "./0_data/OAS_DB_other_virus/Ton_Obs",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))
df_Ton_Obs <- bind_rows(df_list)

              
file_paths <- list.files(path = "./0_data/OAS_DB_other_virus/Tonsilitis",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))
df_Tonsilitis <- bind_rows(df_list)

file_paths <- list.files(path = "./0_data/OAS_DB_other_virus/Mulitple_sclerosis",
                         pattern = "\\.csv\\.gz$", full.names = TRUE)
df_list <- lapply(file_paths, function(fp) read_csv(fp, skip = 1))
df_Mulitple_sclerosis <- bind_rows(df_list)

#验证列名格式是否一致
length(colnames(df_CMV))
length(colnames(df_HIV))
length(colnames(df_Mulitple_sclerosis))
length(colnames(df_Obstructive_Sleep_Apnea))
length(colnames(df_Ton_Obs))
length(colnames(df_Tonsilitis))
sum(colnames(df_HIV) %in% colnames(df_CMV))
sum(colnames(df_Mulitple_sclerosis) %in% colnames(df_Obstructive_Sleep_Apnea))
sum(colnames(df_Mulitple_sclerosis) %in% colnames(df_HIV))
sum(colnames(df_Ton_Obs) %in% colnames(df_Tonsilitis))
sum(colnames(df_HIV) %in% colnames(df_Tonsilitis))

df_CMV_2 <- df_CMV %>% select(which(colnames(df_CMV) %in% colnames(df_HIV)))
df_Mulitple_sclerosis_2 <- df_Mulitple_sclerosis %>% select(which(colnames(df_Mulitple_sclerosis) %in% colnames(df_HIV)))

df_CMV_2$Virus_name <- "CMV"
df_HIV$Virus_name <- "HIV"
df_Mulitple_sclerosis_2$Virus_name <- "Mulitple_sclerosis"
df_Obstructive_Sleep_Apnea$Virus_name <- "Obstructive_Sleep_Apnea"
df_Ton_Obs$Virus_name <- "Ton_Obs"
df_Tonsilitis$Virus_name <- "Tonsilitis"

#验证列名格式是否一致
length(colnames(df_CMV_2))
length(colnames(df_HIV))
length(colnames(df_Mulitple_sclerosis_2))
length(colnames(df_Obstructive_Sleep_Apnea))
length(colnames(df_Ton_Obs))
length(colnames(df_Tonsilitis))
sum(colnames(df_HIV) %in% colnames(df_CMV_2))
sum(colnames(df_Mulitple_sclerosis_2) %in% colnames(df_Obstructive_Sleep_Apnea))
sum(colnames(df_Ton_Obs) %in% colnames(df_Tonsilitis))
sum(colnames(df_HIV) %in% colnames(df_Tonsilitis))

BCR_Other_virus_data <- bind_rows(df_CMV_2,df_HIV,df_Mulitple_sclerosis_2,
                           df_Obstructive_Sleep_Apnea,df_Ton_Obs,df_Tonsilitis)

nrow(BCR_Other_virus_data)

write.csv(BCR_Other_virus_data,"2_BCR_Other_virus_data.csv")

BCR_COVID_data$Virus_name <- "COVID"

sum(colnames(BCR_COVID_data) %in% colnames(BCR_Other_virus_data))

BCR_COVID_data <- BCR_COVID_data %>% select(which(colnames(BCR_COVID_data) %in% colnames(BCR_Other_virus_data)))

colnames(BCR_COVID_data)

All_virus_data <- bind_rows(BCR_COVID_data,BCR_Other_virus_data)
nrow(All_virus_data)

 select_info %in% colnames(df_HIV)

head(df_HIV$sequence_id_heavy,n=1)

select_info <- c("v_call_heavy","d_call_heavy","j_call_heavy","cdr3_heavy","cdr3_aa_heavy","Isotype_heavy",
   "v_call_light","d_call_light","j_call_light","cdr3_light","cdr3_aa_light","Isotype_light",'Virus_name')

select_info

All_virus_data <-All_virus_data %>% select(which(colnames(All_virus_data) %in% select_info))

colnames(All_virus_data)

write.csv(All_virus_data,"2_All_virus_data_select_info.csv")

nrow(All_virus_data)

Bcell <- readRDS("1_Bcell_BCR.rds")

All_virus_data <- read.csv("2_All_virus_data_select_info.csv")

table(All_virus_data$Virus_name)

head(All_virus_data)


