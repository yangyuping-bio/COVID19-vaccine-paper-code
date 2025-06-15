library("tidyverse")
library("data.table")
library(ggplot2)
library("Matrix")
library(Seurat)
library(viridis)
library(RColorBrewer)
library("ggpubr")
library(pacman)
library(fuzzyjoin)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(dplyr)
# 0.环境准备 --------------------------------------------------------------------
library("tidyverse")

VDJ_db<- read.table("data/SearchTable-2024-09-30.tsv", header = TRUE, sep = "\t", fill = TRUE)
VDJ_db <- VDJ_db%>% filter(Species == "HomoSapiens")
VDJ_db <- VDJ_db[, c(2:5, 11:12,15:17)]
#colnames(VDJ_db)[which(colnames(VDJ_db) == "Epitope.gene")] <- "Epitope.gene"
VDJ_db <- VDJ_db%>% filter(Epitope.species == "SARS-CoV-2")
nrow(VDJ_db)
table(VDJ_db$Gene)
table(VDJ_db$Score)
table(VDJ_db$Epitope.gene)
head(VDJ_db)

table(VDJ_db$Epitope.gene,VDJ_db$Score)

CD4T<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD4T_TCR.rds")
CD8T<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD8T_TCR.rds")

CD4T_data <- CD4T@meta.data %>%
                  filter(! is.na(tcr_cdr3))%>% 
                  select(tcr_barcode,ct_level_4,Donor,Condition,
                         tcr_chain,tcr_v_gene,tcr_d_gene,tcr_j_gene,
                         tcr_cdr3s_aa,tcr_cdr3_nt,tcr_cdr3)
CD8T_data <- CD8T@meta.data %>%
                  filter(! is.na(tcr_cdr3))%>% 
                  select(tcr_barcode,ct_level_4,Donor,Condition,
                         tcr_chain,tcr_v_gene,tcr_d_gene,tcr_j_gene,
                         tcr_cdr3s_aa,tcr_cdr3_nt,tcr_cdr3)

head(CD4T_data$tcr_cdr3s_aa)

# 提取第一条TRB链的函数
extract_first_trb <- function(cdr3_string) {
  # 使用正则表达式查找所有TRB链
  trb_matches <- unlist(regmatches(cdr3_string, gregexpr("TRB:[^;]+", cdr3_string)))
  # 如果存在TRB链，返回第一条；否则返回NA
  if (length(trb_matches) > 0) {
    return(trb_matches[1])
  } else {
    return(NA)
  }
}

extract_first_trb_or_tra <- function(cdr3_string) {
  # 查找所有TRB链
  trb_matches <- unlist(regmatches(cdr3_string, gregexpr("TRB:[^;]+", cdr3_string)))
  
  if (length(trb_matches) > 0) {
    return(trb_matches[1])
  } else {
    # 若没有TRB，查找第一条TRA链
    tra_matches <- unlist(regmatches(cdr3_string, gregexpr("TRA:[^;]+", cdr3_string)))
    if (length(tra_matches) > 0) {
      return(tra_matches[1])
    } else {
      return(NA)
    }
  }
}

# 应用函数提取所有数据中的第一条TRB链
CD4T_data$tcr_TRB_cdr3<- sapply(CD4T_data$tcr_cdr3s_aa, extract_first_trb)

# 查看结果
head(CD4T_data$tcr_TRB_cdr3)

# 检查是否存在 NA 值
any_na <- any(is.na(CD4T_data$tcr_TRB_cdr3))
cat("是否存在 NA 值：", any_na, "\n")


# 统计 NA 的总数
na_count <- sum(is.na(CD4T_data$tcr_TRB_cdr3))
cat("总共有", na_count, "个 NA 值\n")
total <- length(CD4T_data$tcr_TRB_cdr3)
cat("总行数：", total, "\n")
cat("NA数量：", na_count, "\n")
cat("NA比例：", round(na_count / total * 100, 2), "%\n")

# 应用函数提取所有数据中的第一条TRB链
CD4T_data$tcr_TRB_or_TRA_cdr3<- sapply(CD4T_data$tcr_cdr3s_aa, extract_first_trb_or_tra)

# 查看结果
head(CD4T_data$tcr_TRB_or_TRA_cdr3)

# 检查是否存在 NA 值
any_na <- any(is.na(CD4T_data$tcr_TRB_or_TRA_cdr3))
cat("是否存在 NA 值：", any_na, "\n")

# 若有 NA 值，可以查看具体哪些是 NA
if (any_na) {
  cat("以下行存在 NA 值：\n")
  print(which(is.na(CD4T_data$tcr_TRB_or_TRA_cdr3)))
}

library(tidyr)

# 假设你的列名是 tcr_TRB_or_TRA_cdr3
CD4T_data <- CD4T_data |>
  separate(tcr_TRB_or_TRA_cdr3, into = c("tcr_TRB_or_TRA_chain", "tcr_TRB_or_TRA_cdr3"), sep = ":", remove = FALSE)
head(CD4T_data,n=1)


#每一种序列的统计数
clonotype_num <- as.data.frame(table(CD4T_data$tcr_TRB_or_TRA_cdr3))
table(clonotype_num$Freq)
CD4T_data$Clonotype_num_b <- clonotype_num$Freq[match(CD4T_data$tcr_TRB_or_TRA_cdr3, clonotype_num$Var1)]


#每个 Condition 下每一种序列的统计数
clonotype_num_Condition <- CD4T_data[ ,c("tcr_TRB_or_TRA_cdr3","Condition")] %>% 
                  group_by(tcr_TRB_or_TRA_cdr3,Condition) %>%
                  mutate(Clonotype_num_b_Condition = n())
nrow(CD4T_data)
nrow(clonotype_num_Condition)
head(clonotype_num_Condition)

## 筛选
nrow(CD4T_data)
CD4T_data$Clontype_num_b_Condition <- clonotype_num_Condition$Clonotype_num_b_Condition
CD4T_data <-  CD4T_data%>%  filter (!(Condition =="A0" & Clontype_num_b_Condition >1))
nrow(CD4T_data)


group_index <- ceiling(seq(1, nrow(VDJ_db)) / 5000)
split_datasets <- split(VDJ_db, group_index) # 使用分组索引切割数据
process_data <- function(data) {
  CD4T_data %>% stringdist_join(data, by=c("tcr_TRB_or_TRA_cdr3"="CDR3"), method="lv",
                  mode="inner", max_dist=2, distance_col="Distance_cdr3")
}

results <- lapply(split_datasets, process_data)
res_distance_2 <- do.call(rbind, results)

res_distance_0 <- res_distance_2%>% filter(Distance_cdr3 ==0)
nrow(res_distance_0)
length(unique(res_distance_0$tcr_cdr3))
write.csv(res_distance_0,"data_COVID/1_CD4T_CDR3_distance_0.csv")

res_distance_1 <- res_distance_2%>% filter(Distance_cdr3 <=1)
nrow(res_distance_1)
length(unique(res_distance_1$tcr_cdr3))
write.csv(res_distance_1,"data_COVID/1_CD4T_CDR3_distance_1.csv")

nrow(res_distance_2)
length(unique(res_distance_2$tcr_cdr3))
write.csv(res_distance_2,"data_COVID/1_CD4T_CDR3_distance_2.csv")

table(res_distance_0$Epitope.gene,res_distance_0$ct_level_4)
table(res_distance_1$Epitope.gene,res_distance_1$ct_level_4)
table(res_distance_2$Epitope.gene,res_distance_2$ct_level_4)

# 应用函数提取所有数据中的第一条TRB链
CD8T_data$tcr_TRB_cdr3<- sapply(CD8T_data$tcr_cdr3s_aa, extract_first_trb)

# 查看结果
head(CD8T_data$tcr_TRB_cdr3)

# 检查是否存在 NA 值
any_na <- any(is.na(CD8T_data$tcr_TRB_cdr3))
cat("是否存在 NA 值：", any_na, "\n")


# 统计 NA 的总数
na_count <- sum(is.na(CD8T_data$tcr_TRB_cdr3))
cat("总共有", na_count, "个 NA 值\n")
total <- length(CD8T_data$tcr_TRB_cdr3)
cat("总行数：", total, "\n")
cat("NA数量：", na_count, "\n")
cat("NA比例：", round(na_count / total * 100, 2), "%\n")

# 应用函数提取所有数据中的第一条TRB链
CD8T_data$tcr_TRB_or_TRA_cdr3<- sapply(CD8T_data$tcr_cdr3s_aa, extract_first_trb_or_tra)

# 查看结果
head(CD8T_data$tcr_TRB_or_TRA_cdr3)

# 检查是否存在 NA 值
any_na <- any(is.na(CD8T_data$tcr_TRB_or_TRA_cdr3))
cat("是否存在 NA 值：", any_na, "\n")

# 若有 NA 值，可以查看具体哪些是 NA
if (any_na) {
  cat("以下行存在 NA 值：\n")
  print(which(is.na(CD8T_data$tcr_TRB_or_TRA_cdr3)))
}

library(tidyr)

# 假设你的列名是 tcr_TRB_or_TRA_cdr3
CD8T_data <- CD8T_data |>
  separate(tcr_TRB_or_TRA_cdr3, into = c("tcr_TRB_or_TRA_chain", "tcr_TRB_or_TRA_cdr3"), sep = ":", remove = FALSE)
head(CD8T_data,n=1)


#每一种序列的统计数
clonotype_num <- as.data.frame(table(CD8T_data$tcr_TRB_or_TRA_cdr3))
table(clonotype_num$Freq)
CD8T_data$Clonotype_num_b <- clonotype_num$Freq[match(CD8T_data$tcr_TRB_or_TRA_cdr3, clonotype_num$Var1)]


#每个 Condition 下每一种序列的统计数
clonotype_num_Condition <- CD8T_data[ ,c("tcr_TRB_or_TRA_cdr3","Condition")] %>% 
                  group_by(tcr_TRB_or_TRA_cdr3,Condition) %>%
                  mutate(Clonotype_num_b_Condition = n())
nrow(CD8T_data)
nrow(clonotype_num_Condition)
head(clonotype_num_Condition)

## 筛选
nrow(CD8T_data)
CD8T_data$Clontype_num_b_Condition <- clonotype_num_Condition$Clonotype_num_b_Condition
CD8T_data <-  CD8T_data%>%  filter (!(Condition =="A0" & Clontype_num_b_Condition >1))
nrow(CD8T_data)


#更严格的匹配
CD8T_data$tcr_v_j_cdr3 <- paste(paste(CD8T_data$tcr_v_gene, "01", sep = "*"),
                                paste(CD8T_data$tcr_j_gene, "01", sep = "*"),
                                CD8T_data$tcr_TRB_or_TRA_cdr3, sep = "_")

VDJ_db$tcr_v_j_cdr3 <- paste(VDJ_db$V,VDJ_db$J,VDJ_db$CDR3, sep = "_")

group_index <- ceiling(seq(1, nrow(VDJ_db)) / 10000)
split_datasets <- split(VDJ_db, group_index)# 使用分组索引切割数据
process_data <- function(data) {
  CD8T_data %>% stringdist_join(data, by=c("tcr_v_j_cdr3"="tcr_v_j_cdr3"), method="lv",
                  mode="inner", max_dist=2, distance_col="Distance_cdr3")
}

results <- lapply(split_datasets, process_data)
res_distance_2 <- do.call(rbind, results)

res_distance_0 <- res_distance_2%>% filter(Distance_cdr3 ==0)
nrow(res_distance_0)
length(unique(res_distance_0$tcr_cdr3))
write.csv(res_distance_0,"data_COVID/2_CD8T_VJ_CDR3_distance_0.csv")

res_distance_1 <- res_distance_2%>% filter(Distance_cdr3 <=1)
nrow(res_distance_1)
length(unique(res_distance_1$tcr_cdr3))
write.csv(res_distance_1,"data_COVID/2_CD8T_VJ_CDR3_distance_1.csv")

nrow(res_distance_2)
length(unique(res_distance_2$tcr_cdr3))
write.csv(res_distance_2,"data_COVID/2_CD8T_VJ_CDR3_distance_2.csv")

table(res_distance_0$Epitope.gene,res_distance_0$ct_level_4)
table(res_distance_1$Epitope.gene,res_distance_1$ct_level_4)
table(res_distance_2$Epitope.gene,res_distance_2$ct_level_4)

group_index <- ceiling(seq(1, nrow(VDJ_db)) / 5000)
split_datasets <- split(VDJ_db, group_index)# 使用分组索引切割数据
process_data <- function(data) {
  CD8T_data %>% stringdist_join(data, by=c("tcr_TRB_or_TRA_cdr3"="CDR3"), method="lv",
                  mode="inner", max_dist=2, distance_col="Distance_cdr3")
}

results <- lapply(split_datasets, process_data)
res_distance_2 <- do.call(rbind, results)

res_distance_0 <- res_distance_2%>% filter(Distance_cdr3 ==0)
nrow(res_distance_0)
length(unique(res_distance_0$tcr_cdr3))
write.csv(res_distance_0,"data_COVID/2_CD8T_CDR3_distance_0.csv")

res_distance_1 <- res_distance_2%>% filter(Distance_cdr3 <=1)
nrow(res_distance_1)
length(unique(res_distance_1$tcr_cdr3))
write.csv(res_distance_1,"data_COVID/2_CD8T_CDR3_distance_1.csv")

nrow(res_distance_2)
length(unique(res_distance_2$tcr_cdr3))
write.csv(res_distance_2,"data_COVID/2_CD8T_CDR3_distance_2.csv")

table(res_distance_0$Epitope.gene,res_distance_0$ct_level_4)
table(res_distance_1$Epitope.gene,res_distance_1$ct_level_4)
table(res_distance_2$Epitope.gene,res_distance_2$ct_level_4)

CD8T_res_distance_1 <- read.csv("data_COVID/2_CD8T_CDR3_distance_1.csv")
#CD8T_res_distance_1 <- read.csv("data_COVID/2_CD8T_VJ_CDR3_distance_2.csv")
CD4T_res_distance_1 <- read.csv("data_COVID/1_CD4T_CDR3_distance_1.csv")
nrow(CD4T_res_distance_1)
nrow(CD8T_res_distance_1)

# 第一步：优先保留 Distance_cdr3 == 0 的行（如果存在），否则保留原始行
filtered_by_distance <- CD8T_res_distance_1 %>%
  group_by(tcr_barcode) %>%
  group_modify(~{
    df_group <- .
    if (any(df_group$Distance_cdr3 == 0)) {
      df_group[df_group$Distance_cdr3 == 0, , drop = FALSE]
    } else {
      df_group
    }
  }) %>%
  ungroup()

# 第二步：对每个 barcode，选择 Epitope.gene 出现频率最多的那一个，再只保留一行
CD8T_res_filtered <- filtered_by_distance %>%
  group_by(tcr_barcode) %>%
  summarise(
    # 获取频率最高的 Epitope.gene
    Epitope.gene = names(sort(table(Epitope.gene), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  # 与原数据做 join，以获得其他列
  left_join(filtered_by_distance, by = c("tcr_barcode", "Epitope.gene")) %>%
  group_by(tcr_barcode) %>%
  slice(1) %>%  # 保底操作：如果 join 后仍有重复，只保留第一行
  ungroup()

# 去除 Epitope.gene 字段的首尾空格
CD8T_res_filtered$Epitope.gene <- CD8T_res_filtered$Epitope.gene %>%
  as.character() %>%
  trimws()

# 验证结果是否符合预期
cat("唯一 barcode 个数: ", length(unique(CD8T_res_distance_1$tcr_barcode)), "\n")
cat("最终筛选行数: ", nrow(CD8T_res_filtered), "\n")
cat("是否存在重复 barcode: ", any(duplicated(CD8T_res_filtered$tcr_barcode)), "\n")

# 查看结果
head(CD8T_res_filtered, 1)
nrow(CD8T_res_filtered)
table(CD8T_res_filtered$Epitope.gene)

library(dplyr)
# 第一步：优先保留 Distance_cdr3 == 0 > 1 > 2 的行（分组处理）
filtered_by_distance <- CD8T_res_distance_1 %>%
  group_by(tcr_barcode) %>%
  group_modify(~{
    df_group <- .
    if (any(df_group$Distance_cdr3 == 0)) {
      df_group[df_group$Distance_cdr3 == 0, , drop = FALSE]
    } else if (any(df_group$Distance_cdr3 == 1)) {
      df_group[df_group$Distance_cdr3 == 1, , drop = FALSE]
    } else {
      df_group
    }
  }) %>%
  ungroup()

# 第二步：只保留那些 Epitope.gene 唯一的 barcode
filtered_unique_gene <- filtered_by_distance %>%
  group_by(tcr_barcode) %>%
  filter(n_distinct(Epitope.gene) == 1) %>%
  ungroup()

# 第三步：每个 barcode 保留 Distance_cdr3 最小的那一行（即唯一行）
CD8T_res_filtered <- filtered_unique_gene %>%
  group_by(tcr_barcode) %>%
  slice_min(order_by = Distance_cdr3, n = 1, with_ties = FALSE) %>%
  ungroup()

# 去除 Epitope.gene 首尾空格
CD8T_res_filtered$Epitope.gene <- CD8T_res_filtered$Epitope.gene %>%
  as.character() %>%
  trimws()

# 验证结果是否符合预期
cat("唯一barcode 个数: ", length(unique(CD8T_res_distance_1$tcr_barcode)), "\n")
cat("最终筛选行数: ", nrow(CD8T_res_filtered), "\n")
cat("是否存在重复 barcode: ", any(duplicated(CD8T_res_filtered$tcr_barcode)), "\n")

# 查看结果
head(CD8T_res_filtered, 1)
nrow(CD8T_res_filtered)
table(CD8T_res_filtered$Epitope.gene)

library(dplyr)

filtered_by_distance_and_score <- CD8T_res_distance_1 %>%
  # 第一步：去除 Score == 0 的行
  filter(Score > 0) %>%
  group_by(tcr_barcode) %>%
  group_modify(~{
    df_group <- .
    # 优先选择 Distance_cdr3 从小到大
    for (d in 0:2) {
      df_d <- df_group %>% filter(Distance_cdr3 == d)
      if (nrow(df_d) > 0) {
        # 在当前 Distance_cdr3 下，再按 Score 3 > 2 > 1 进行优先筛选
        for (s in 3:1) {
          df_ds <- df_d %>% filter(Score == s)
          if (nrow(df_ds) > 0) {
            return(df_ds)
          }
        }
      }
    }
    # 如果 Distance_cdr3 > 2，或者没有匹配上，就返回空表
    df_group[0, ]
  }) %>%
  ungroup()

# 第二步：只保留那些 Epitope.gene 唯一的 barcode
filtered_unique_gene <- filtered_by_distance_and_score %>%
  group_by(tcr_barcode) %>%
  filter(n_distinct(Epitope.gene) == 1) %>%
  ungroup()

# 第三步：每个 barcode 保留 Distance_cdr3 最小的那一行（唯一行）
CD8T_res_filtered <- filtered_unique_gene %>%
  group_by(tcr_barcode) %>%
  slice_min(order_by = Distance_cdr3, n = 1, with_ties = FALSE) %>%
  ungroup()

# 去除 Epitope.gene 首尾空格
CD8T_res_filtered$Epitope.gene <- CD8T_res_filtered$Epitope.gene %>%
  as.character() %>%
  trimws()

# 验证结果
cat("唯一 barcode 个数: ", length(unique(CD8T_res_distance_1$tcr_barcode)), "\n")
cat("最终筛选行数: ", nrow(CD8T_res_filtered), "\n")
cat("是否存在重复 barcode: ", any(duplicated(CD8T_res_filtered$tcr_barcode)), "\n")

# 查看结果
#head(CD8T_res_filtered, 1)
nrow(CD8T_res_filtered)
table(CD8T_res_filtered$Epitope.gene)

table(CD8T_res_filtered$ct_level_4,CD8T_res_filtered$Epitope.gene)

table(CD8T_res_filtered$Epitope.gene,CD8T_res_filtered$Condition)

table(CD8T_res_filtered$ct_level_4,CD8T_res_filtered$Score)

library(dplyr)

filtered_by_distance_and_score <- CD4T_res_distance_1 %>%
  # 第一步：去除 Score == 0 的行
  filter(Score > 0) %>%
  group_by(tcr_barcode) %>%
  group_modify(~{
    df_group <- .
    # 优先选择 Distance_cdr3 从小到大
    for (d in 0:2) {
      df_d <- df_group %>% filter(Distance_cdr3 == d)
      if (nrow(df_d) > 0) {
        # 在当前 Distance_cdr3 下，再按 Score 3 > 2 > 1 进行优先筛选
        for (s in 3:1) {
          df_ds <- df_d %>% filter(Score == s)
          if (nrow(df_ds) > 0) {
            return(df_ds)
          }
        }
      }
    }
    # 如果 Distance_cdr3 > 2，或者没有匹配上，就返回空表
    df_group[0, ]
  }) %>%
  ungroup()

# 第二步：只保留那些 Epitope.gene 唯一的 barcode
filtered_unique_gene <- filtered_by_distance_and_score %>%
  group_by(tcr_barcode) %>%
  filter(n_distinct(Epitope.gene) == 1) %>%
  ungroup()

# 第三步：每个 barcode 保留 Distance_cdr3 最小的那一行（唯一行）
CD4T_res_filtered <- filtered_unique_gene %>%
  group_by(tcr_barcode) %>%
  slice_min(order_by = Distance_cdr3, n = 1, with_ties = FALSE) %>%
  ungroup()

# 去除 Epitope.gene 首尾空格
CD4T_res_filtered$Epitope.gene <- CD4T_res_filtered$Epitope.gene %>%
  as.character() %>%
  trimws()

# 验证结果
cat("唯一 barcode 个数: ", length(unique(CD4T_res_distance_1$tcr_barcode)), "\n")
cat("最终筛选行数: ", nrow(CD4T_res_filtered), "\n")
cat("是否存在重复 barcode: ", any(duplicated(CD4T_res_filtered$tcr_barcode)), "\n")

# 查看结果
#head(CD8T_res_filtered, 1)
nrow(CD4T_res_filtered)
table(CD4T_res_filtered$Epitope.gene)

table(CD4T_res_filtered$ct_level_4,CD4T_res_filtered$Epitope.gene)

table(CD4T_res_filtered$Epitope.gene,CD4T_res_filtered$Condition)

table(CD4T_res_filtered$ct_level_4,CD4T_res_filtered$Score)

library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

# 计算频数和百分比
species_table <- CD8T_res_filtered %>%
  count(Epitope.gene, name = "Count") %>%
  mutate(
    Percentage = Count / sum(Count),
    Epitope.gene = factor(Epitope.gene, levels = Epitope.gene[order(-Count)])
  )

# 分配颜色（统一顺序）
colors <- setNames(
  colorRampPalette(rev(brewer.pal(12, "Set3")))(nlevels(species_table$Epitope.gene)),
  levels(species_table$Epitope.gene)
)

pdf("figure/2_CD8T_Epitope_Gene_Pie_and_Bar.pdf", 5, 4)
# ==== 饼图 ====
ggplot(species_table, aes(x = "", y = Percentage, fill = Epitope.gene)) +
  geom_bar(width = 0.9, stat = "identity", alpha = 0.95, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colors) +
  theme_void() +
  labs(fill = "Epitope Gene", title = "Epitope Gene Distribution in CD8T") +
  theme(
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.position = 'right',
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))

# ==== 柱状图 ====
ggplot(species_table, aes(y = reorder(Epitope.gene, Count), x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity", alpha = 0.95, color = "white") +
  geom_text(aes(label = Count), hjust = -0.2, color = "black", size = 3.5) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(fill = "Epitope Gene", title = "Epitope Gene Frequency in CD8T", y = "", x = "") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),
    legend.position = 'none',
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
  ) +
  coord_cartesian(clip = "off")
dev.off()

library(ggplot2)
library(dplyr)

# ===== 准备数据（Epitope.gene 总排序）=====
gene_order <- CD8T_res_filtered %>%
  count(Epitope.gene, name = "Count") %>%
  arrange(desc(Count)) %>%
  pull(Epitope.gene)

# 统一颜色（按排序）
colors <- setNames(
  colorRampPalette(rev(RColorBrewer::brewer.pal(12, "Set3")))(length(gene_order)),
  gene_order
)

df1 <- as.data.frame(table(CD8T_res_filtered$Condition, CD8T_res_filtered$Epitope.gene))
colnames(df1) <- c("Condition", "Epitope.gene", "Count")
df1$Epitope.gene <- factor(df1$Epitope.gene, levels = gene_order)
df2 <- as.data.frame(table(CD8T_res_filtered$ct_level_4, CD8T_res_filtered$Epitope.gene))
colnames(df2) <- c("ct_level_4", "Epitope.gene", "Count")
df2$Epitope.gene <- factor(df2$Epitope.gene, levels = gene_order)

# 图1：Epitope.gene 堆叠在 Y = Condition 上（水平堆叠条形图）
pdf("figure/2_CD8T_Condition_and_celltype_EpitopeGene.pdf", 8, 4.2)
ggplot(df1, aes(y = Condition, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "fill",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by Condition", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )

# 图2：Epitope.gene 堆叠在 Y = ct_level_4 上（水平堆叠条形图）
ggplot(df2, aes(y = ct_level_4, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "fill",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by ct_level_4", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )
# 图1：Epitope.gene 堆叠在 Y = Condition 上（水平堆叠条形图）
ggplot(df1, aes(y = Condition, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "stack",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by Condition", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )
# 图3：Epitope.gene 堆叠在 Y = ct_level_4 上（水平堆叠条形图）
ggplot(df2, aes(y = ct_level_4, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "stack",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by ct_level_4", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )
dev.off()





library(dplyr)

# 先找出只对应一种 Epitope.gene 的 tcr_barcode
unique_epitope_tcrs <-CD4T_res_distance_1 %>%
  group_by(tcr_barcode) %>%
  filter(n_distinct(Epitope.gene) == 1) %>%
  ungroup()

# 如果还需要去重，每个 barcode 保留一行（比如 Distance_cdr3 最小的）
CD4T_res_filtered <- unique_epitope_tcrs %>%
  group_by(tcr_barcode) %>%
  slice_min(order_by = Distance_cdr3, n = 1, with_ties = FALSE) %>%
  ungroup()

# 清理列
CD4T_res_filtered$Epitope.gene <-CD4T_res_filtered$Epitope.gene %>%
  as.character() %>%
  trimws() ## 去除首尾空格

# 查看结果
head(CD4T_res_filtered, 1)
nrow(CD4T_res_filtered)
table(CD4T_res_filtered$Epitope.gene)

library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

# 计算频数和百分比
species_table <- CD4T_res_filtered %>%
  count(Epitope.gene, name = "Count") %>%
  mutate(
    Percentage = Count / sum(Count),
    Epitope.gene = factor(Epitope.gene, levels = Epitope.gene[order(-Count)])
  )

# 分配颜色（统一顺序）
colors <- setNames(
  colorRampPalette(rev(brewer.pal(12, "Set3")))(nlevels(species_table$Epitope.gene)),
  levels(species_table$Epitope.gene)
)

pdf("figure/1_CD4T_Epitope_Gene_Pie_and_Bar.pdf", 5, 4)
# ==== 饼图 ====
ggplot(species_table, aes(x = "", y = Percentage, fill = Epitope.gene)) +
  geom_bar(width = 0.9, stat = "identity", alpha = 0.95, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colors) +
  theme_void() +
  labs(fill = "Epitope Gene", title = "Epitope Gene Distribution in CD4T") +
  theme(
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.position = 'right',
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))

# ==== 柱状图 ====
ggplot(species_table, aes(y = reorder(Epitope.gene, Count), x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity", alpha = 0.95, color = "white") +
  geom_text(aes(label = Count), hjust = -0.2, color = "black", size = 3.5) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(fill = "Epitope Gene", title = "Epitope Gene Frequency in CD4T", y = "", x = "") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),
    legend.position = 'none',
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
  ) +
  coord_cartesian(clip = "off")
dev.off()

table(CD4T_res_filtered$Condition,CD4T_res_filtered$Epitope.gene)
table(CD4T_res_filtered$ct_level_4,CD4T_res_filtered$Epitope.gene)

library(ggplot2)
library(dplyr)

# ===== 准备数据（Epitope.gene 总排序）=====
gene_order <- CD4T_res_filtered %>%
  count(Epitope.gene, name = "Count") %>%
  arrange(desc(Count)) %>%
  pull(Epitope.gene)

# 统一颜色（按排序）
colors <- setNames(
  colorRampPalette(rev(RColorBrewer::brewer.pal(12, "Set3")))(length(gene_order)),
  gene_order
)

df1 <- as.data.frame(table(CD4T_res_filtered$Condition, CD4T_res_filtered$Epitope.gene))
colnames(df1) <- c("Condition", "Epitope.gene", "Count")
df1$Epitope.gene <- factor(df1$Epitope.gene, levels = gene_order)
df2 <- as.data.frame(table(CD4T_res_filtered$ct_level_4, CD4T_res_filtered$Epitope.gene))
colnames(df2) <- c("ct_level_4", "Epitope.gene", "Count")
df2$Epitope.gene <- factor(df2$Epitope.gene, levels = gene_order)

# 图1：Epitope.gene 堆叠在 Y = Condition 上（水平堆叠条形图）
pdf("figure/1_CD4T_Condition_and_celltype_EpitopeGene.pdf", 8, 4.2)
ggplot(df1, aes(y = Condition, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "fill",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by Condition", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )

# 图2：Epitope.gene 堆叠在 Y = ct_level_4 上（水平堆叠条形图）
ggplot(df2, aes(y = ct_level_4, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "fill",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by ct_level_4", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )
# 图1：Epitope.gene 堆叠在 Y = Condition 上（水平堆叠条形图）
ggplot(df1, aes(y = Condition, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "stack",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by Condition", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )
# 图3：Epitope.gene 堆叠在 Y = ct_level_4 上（水平堆叠条形图）
ggplot(df2, aes(y = ct_level_4, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  position = "stack",  alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition by ct_level_4", y = "", x = "Count", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "char"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
       legend.key.size = unit(0.2, "inch")
  )
dev.off()

# 1. 筛选出 ct_level_4 为 "CD4_Treg" 的数据
treg_data <- CD4T_res_filtered %>% 
  filter(ct_level_4 == "CD4_Treg")

# 2. 统计 Condition vs Epitope.gene
treg_table <- treg_data %>%
  count(Condition, Epitope.gene, name = "Count")

# 3. 设定颜色（保持和之前一致）
epitope_levels <- treg_table %>%
  group_by(Epitope.gene) %>%
  summarise(Total = sum(Count)) %>%
  arrange(-Total) %>%
  pull(Epitope.gene)
treg_table$Epitope.gene <- factor(treg_table$Epitope.gene, levels = epitope_levels)

# 4. 画图（横向相对比例堆叠条形图）
pdf("figure/1_CD4T_Treg_Condition_vs_EpitopeGene_relative.pdf", 6, 4)
ggplot(treg_table, aes(y = Condition, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity", position = "fill", alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition in CD4_Treg (by Condition)", 
       x = "Proportion", y = "", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(2, 2, 2, 2), "char")
  )
dev.off()

# 1. 筛选出 ct_level_4 为 "CD4_Treg" 的数据
treg_data <- CD4T_res_filtered %>% 
  filter(ct_level_4 == "CD4_Treg")

# 2. 统计 Condition vs Epitope.gene
treg_table <- treg_data %>%
  count(Condition, Epitope.gene, name = "Count")

# 3. 设定颜色（保持和之前一致）
epitope_levels <- treg_table %>%
  group_by(Epitope.gene) %>%
  summarise(Total = sum(Count)) %>%
  arrange(-Total) %>%
  pull(Epitope.gene)
treg_table$Epitope.gene <- factor(treg_table$Epitope.gene, levels = epitope_levels)

# 4. 画图（横向相对比例堆叠条形图）
pdf("figure/1_CD4T_Treg_Condition_vs_EpitopeGene_relative_2.pdf", 6, 4)
ggplot(treg_table, aes(y = Condition, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity",  alpha = 0.95, color = "white") + #position = "fill",
  scale_fill_manual(values = colors) +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition in CD4_Treg (by Condition)", 
       x = "Proportion", y = "", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(2, 2, 2, 2), "char")
  )
dev.off()

# 4. 画图（横向堆叠条形图，显示数量）
pdf("figure/1_CD4T_Treg_Condition_vs_EpitopeGene_absolute.pdf", 6, 4)
ggplot(treg_table, aes(y = Condition, x = Count, fill = Epitope.gene)) +
  geom_bar(stat = "identity", alpha = 0.95, color = "white") +
  scale_fill_manual(values = colors) +
  # 删除百分比格式
  # scale_x_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  labs(title = "Epitope Gene Composition in CD4_Treg (by Condition)", 
       x = "Count", y = "", fill = "Epitope Gene") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(2, 2, 2, 2), "char")
  )
dev.off()



# Epitope Species ---------------------------------------------------------

library("tidyverse")
library("data.table")
library(ggplot2)
library("Matrix")
library(Seurat)
library(viridis)
library(RColorBrewer)
library("ggpubr")
library(pacman)
library(fuzzyjoin)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(dplyr)

getwd()

VDJ_db<- read.table("data/SearchTable-2024-09-30.tsv", header = TRUE, sep = "\t", fill = TRUE)
VDJ_db <- VDJ_db%>% filter(Species == "HomoSapiens")
VDJ_db <- VDJ_db[, 1:(ncol(VDJ_db)-5)] # 删掉后面5列数据
VDJ_db <- VDJ_db[, c(2:5, 11:12)]
nrow(VDJ_db)
head(VDJ_db)

CD4T<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD4T_TCR.rds")
CD8T<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD8T_TCR.rds")

length_distribution <- table(nchar(as.character(CD4T@meta.data$tcr_cdr3)))
barplot(length_distribution, 
        main = "Distribution of sequence length in our data", 
        xlab = "Sequence length", 
        ylab = "Frequency")
length_distribution <- table(nchar(as.character(CD8T@meta.data$tcr_cdr3)))
barplot(length_distribution, 
        main = "Distribution of sequence length in our data", 
        xlab = "Sequence length", 
        ylab = "Frequency")
length_distribution <- table(nchar(as.character(VDJ_db$CDR3)))
barplot(length_distribution, 
        main = "Distribution of sequence length in DB", 
        xlab = "Sequence length", 
        ylab = "Frequency")

CD4T_data <- CD4T@meta.data %>%
  filter(! is.na(tcr_cdr3))%>% 
  select(tcr_barcode,ct_level_4,Donor,Condition,
         tcr_chain,tcr_v_gene,tcr_d_gene,tcr_j_gene,
         tcr_cdr3s_aa,tcr_cdr3_nt,tcr_cdr3)
CD8T_data <- CD8T@meta.data %>%
  filter(! is.na(tcr_cdr3))%>% 
  select(tcr_barcode,ct_level_4,Donor,Condition,
         tcr_chain,tcr_v_gene,tcr_d_gene,tcr_j_gene,
         tcr_cdr3s_aa,tcr_cdr3_nt,tcr_cdr3)

extract_first_trb_or_tra <- function(cdr3_string) {
  # 查找所有TRB链
  trb_matches <- unlist(regmatches(cdr3_string, gregexpr("TRB:[^;]+", cdr3_string)))
  
  if (length(trb_matches) > 0) {
    return(trb_matches[1])
  } else {
    # 若没有TRB，查找第一条TRA链
    tra_matches <- unlist(regmatches(cdr3_string, gregexpr("TRA:[^;]+", cdr3_string)))
    if (length(tra_matches) > 0) {
      return(tra_matches[1])
    } else {
      return(NA)
    }
  }
}

# 应用函数提取所有数据中的第一条TRB链
CD4T_data$tcr_TRB_or_TRA_cdr3<- sapply(CD4T_data$tcr_cdr3s_aa, extract_first_trb_or_tra)

# 查看结果
head(CD4T_data$tcr_TRB_or_TRA_cdr3)

# 检查是否存在 NA 值
any_na <- any(is.na(CD4T_data$tcr_TRB_or_TRA_cdr3))
cat("是否存在 NA 值：", any_na, "\n")

# 若有 NA 值，可以查看具体哪些是 NA
if (any_na) {
  cat("以下行存在 NA 值：\n")
  print(which(is.na(CD4T_data$tcr_TRB_or_TRA_cdr3)))
}

library(tidyr)

# 假设你的列名是 tcr_TRB_or_TRA_cdr3
CD4T_data <- CD4T_data |>
  separate(tcr_TRB_or_TRA_cdr3, into = c("tcr_TRB_or_TRA_chain", "tcr_TRB_or_TRA_cdr3"), sep = ":", remove = FALSE)
head(CD4T_data,n=1)

# 应用函数提取所有数据中的第一条TRB链
CD8T_data$tcr_TRB_or_TRA_cdr3<- sapply(CD8T_data$tcr_cdr3s_aa, extract_first_trb_or_tra)

# 查看结果
head(CD8T_data$tcr_TRB_or_TRA_cdr3)

# 检查是否存在 NA 值
any_na <- any(is.na(CD8T_data$tcr_TRB_or_TRA_cdr3))
cat("是否存在 NA 值：", any_na, "\n")

# 若有 NA 值，可以查看具体哪些是 NA
if (any_na) {
  cat("以下行存在 NA 值：\n")
  print(which(is.na(CD8T_data$tcr_TRB_or_TRA_cdr3)))
}

library(tidyr)

# 假设你的列名是 tcr_TRB_or_TRA_cdr3
CD8T_data <- CD8T_data |>
  separate(tcr_TRB_or_TRA_cdr3, into = c("tcr_TRB_or_TRA_chain", "tcr_TRB_or_TRA_cdr3"), sep = ":", remove = FALSE)
head(CD8T_data,n=1)

num_rows <- nrow(VDJ_db)# 计算数据的行数
group_index <- ceiling(seq(1, num_rows) / 10000)# 创建一个序列, 每10000个数据为一组
split_datasets <- split(VDJ_db, group_index)# 使用分组索引切割数据
process_data <- function(data) {
  CD4T_data %>% stringdist_join(data, by=c("tcr_TRB_or_TRA_cdr3"="CDR3"), method="lv",
                                mode="inner", max_dist=2, distance_col="Distance_cdr3")
}

results <- lapply(split_datasets, process_data)
res_distance_2 <- do.call(rbind, results)

res_distance_0 <- res_distance_2%>% filter(Distance_cdr3 ==0)
nrow(res_distance_0)
length(unique(res_distance_0$tcr_cdr3))
write.csv(res_distance_0,"data/1_CD4T_virus_tcr_distance_0.csv")

res_distance_1 <- res_distance_2%>% filter(Distance_cdr3 <=1)
nrow(res_distance_1)
length(unique(res_distance_1$tcr_cdr3))
write.csv(res_distance_1,"data/1_CD4T_virus_tcr_distance_1.csv")

nrow(res_distance_2)
length(unique(res_distance_2$tcr_cdr3))
write.csv(res_distance_2,"data/1_CD4T_virus_tcr_distance_2.csv")

table(res_distance_0$Epitope.species)

num_rows <- nrow(VDJ_db)# 计算数据的行数
group_index <- ceiling(seq(1, num_rows) / 10000)# 创建一个序列, 每10000个数据为一组
split_datasets <- split(VDJ_db, group_index)# 使用分组索引切割数据
process_data <- function(data) {
  CD8T_data %>% stringdist_join(data, by=c("tcr_TRB_or_TRA_cdr3"="CDR3"), method="lv",
                                mode="inner", max_dist=0, distance_col="Distance_cdr3")
}

results <- lapply(split_datasets, process_data)
res_distance_0 <- do.call(rbind, results)
write.csv(res_distance_0,"data/3_CD8T_virus_tcr_distance_0.csv")
table(res_distance_0$Epitope.species)

CD4T_data$tcr_v_j_cdr3 <- paste(paste(CD4T_data$tcr_v_gene, "01", sep = "*"),
                                paste(CD4T_data$tcr_j_gene, "01", sep = "*"),
                                CD4T_data$tcr_TRB_or_TRA_cdr3, sep = "_")

VDJ_db$tcr_v_j_cdr3 <- paste(VDJ_db$V, VDJ_db$J,VDJ_db$CDR3, sep = "_")

group_index <- ceiling(seq(1, nrow(VDJ_db)) / 10000)
split_datasets <- split(VDJ_db, group_index)# 使用分组索引切割数据
process_data <- function(data) {
  CD4T_data %>% stringdist_join(data, by=c("tcr_v_j_cdr3"="tcr_v_j_cdr3"), method="lv",
                                mode="inner", max_dist=2, distance_col="Distance_cdr3")
}
results <- lapply(split_datasets, process_data)
res_distance_2 <- do.call(rbind, results)

res_distance_0 <- res_distance_2%>% filter(Distance_cdr3 ==0)
nrow(res_distance_0)
length(unique(res_distance_0$tcr_cdr3))
write.csv(res_distance_0,"data/3_CD4T_virus_tcr_distance_tcr_v_j_0.csv")

res_distance_1 <- res_distance_2%>% filter(Distance_cdr3 <=1)
nrow(res_distance_1)
length(unique(res_distance_1$tcr_cdr3))
write.csv(res_distance_1,"data/3_CD4T_virus_tcr_distance_tcr_v_j_1.csv")

nrow(res_distance_2)
length(unique(res_distance_2$tcr_cdr3))
write.csv(res_distance_2,"data/3_CD4T_virus_tcr_distance_tcr_v_j_2.csv")

CD8T_data$tcr_v_j_cdr3 <- paste(paste(CD8T_data$tcr_v_gene, "01", sep = "*"),
                                paste(CD8T_data$tcr_j_gene, "01", sep = "*"),
                                CD8T_data$tcr_TRB_or_TRA_cdr3, sep = "_")

VDJ_db$tcr_v_j_cdr3 <- paste(VDJ_db$V,VDJ_db$J,VDJ_db$CDR3, sep = "_")

group_index <- ceiling(seq(1, nrow(VDJ_db)) / 10000)
split_datasets <- split(VDJ_db, group_index)# 使用分组索引切割数据
process_data <- function(data) {
  CD8T_data %>% stringdist_join(data, by=c("tcr_v_j_cdr3"="tcr_v_j_cdr3"), method="lv",
                                mode="inner", max_dist=2, distance_col="Distance_cdr3")
}

results <- lapply(split_datasets, process_data)
res_distance_2 <- do.call(rbind, results)

res_distance_0 <- res_distance_2%>% filter(Distance_cdr3 ==0)
nrow(res_distance_0)
length(unique(res_distance_0$tcr_cdr3))
write.csv(res_distance_0,"data/4_CD8T_virus_tcr_distance_tcr_v_j_0.csv")

res_distance_1 <- res_distance_2%>% filter(Distance_cdr3 <=1)
nrow(res_distance_1)
length(unique(res_distance_1$tcr_cdr3))
write.csv(res_distance_1,"data/4_CD8T_virus_tcr_distance_tcr_v_j_1.csv")

nrow(res_distance_2)
length(unique(res_distance_2$tcr_cdr3))
write.csv(res_distance_2,"data/4_CD8T_virus_tcr_distance_tcr_v_j_2.csv")

process_data <- function(file_path) {
  # 读取数据
  res_distance <- read.csv(file_path)
  # 统计并标记 TCR 序列的一对多情况
  res_distance <- res_distance %>%
    group_by(tcr_cdr3) %>%
    mutate(cdr3_dis_virus_kind = length(unique(Epitope.species)))
  print("-----------------------------------------------------------------------------------")
  print("原序列的一对多情况统计：")
  print(table(res_distance$cdr3_dis_virus_kind))
  # 过滤掉一对多情况的序列
  res_distance <- res_distance %>% filter(cdr3_dis_virus_kind == 1)
  # 打印清理后的病毒分布
  print("清理后的病毒分布：（全为一一对应的序列）")
  print(table(res_distance$Epitope.species))
  # 过滤出 SARS-CoV-2 的 TCR 序列
  res_distance_COVID <- res_distance %>% filter(Epitope.species == "SARS-CoV-2")
  print(paste("新冠TCR序列数：", nrow(res_distance_COVID)))
  
  # 打印新冠序列的 ORF.Coverage 分布
  print("新冠TCR序列抗原分布：")
  print(table(res_distance_COVID$ORF.Coverage))
  
  # 打印新冠TCR的不同序列数
  print(paste("其中新冠TCR不同序列统计数：", length(unique(res_distance_COVID$tcr_v_j_cdr3.x))))
  return(res_distance)
}

# 处理三个不同的文件
res_distance_0 <- process_data("data/1_CD4T_virus_tcr_distance_0.csv")
res_distance_1 <- process_data("data/1_CD4T_virus_tcr_distance_1.csv")
res_distance_2 <- process_data("data/1_CD4T_virus_tcr_distance_2.csv")
res_distance_0 <- process_data("data/3_CD4T_virus_tcr_distance_tcr_v_j_0.csv")
res_distance_1 <- process_data("data/3_CD4T_virus_tcr_distance_tcr_v_j_1.csv")
res_distance_2 <- process_data("data/3_CD4T_virus_tcr_distance_tcr_v_j_2.csv")

# 处理三个不同的文件
res_distance_0 <- process_data("data/3_CD8T_virus_tcr_distance_0.csv")
res_distance_0 <- process_data("data/4_CD8T_virus_tcr_distance_tcr_v_j_0.csv")
res_distance_1 <- process_data("data/4_CD8T_virus_tcr_distance_tcr_v_j_1.csv")
res_distance_2 <- process_data("data/4_CD8T_virus_tcr_distance_tcr_v_j_2.csv")

res_distance_2 <- process_data("data/2_CD4T_virus_tcr_distance_tcr_v_j_2.csv")
CD4T_data$tcr_v_j_cdr3 <- paste(paste(CD4T_data$tcr_v_gene, "01", sep = "*"),
                                paste(CD4T_data$tcr_j_gene, "01", sep = "*"),
                                CD4T_data$tcr_cdr3, sep = "_")
CD4T_data$ORF.Coverage <- res_distance_2$ORF.Coverage[match(CD4T_data$tcr_v_j_cdr3,res_distance_2$tcr_v_j_cdr3.x)]
CD4T_data$Epitope.species <- res_distance_2$Epitope.species[match(CD4T_data$tcr_v_j_cdr3,res_distance_2$tcr_v_j_cdr3.x)]
CD4T_data$Chain_distance <- res_distance_2$Distance_cdr3[match(CD4T_data$tcr_v_j_cdr3,res_distance_2$tcr_v_j_cdr3.x)]
table(CD4T_data$ORF.Coverage)
table(CD4T_data$Epitope.species)
table(CD4T_data$Chain_distance)
write.csv(CD4T_data,"CD4T_data_res_distance_2.csv")


species_table <- as.data.frame(table(CD4T_data$Epitope.species))
colnames(species_table) <- c("Epitope_species", "Count")

# 计算每种物种的百分比
species_table <- species_table %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  arrange(desc(Percentage))  # 按百分比从高到低排序

species_table$Epitope_species <- factor(species_table$Epitope_species,
                                        levels = unique(species_table$Epitope_species))

colors <- colorRampPalette(brewer.pal(12, "Paired"))(25)# 设置颜色（使用 RColorBrewer 中的调色板）

# 生成饼图
pdf("2_CD4T_Epitope_Species_Pie_distance_2.pdf", 6, 5)
ggplot(species_table, aes(x = "", y = Percentage, fill = Epitope_species)) +
  geom_bar(width = 0.9, stat = "identity", alpha = 0.95, color = "white") +  # 设置透明度
  coord_polar(theta = "y") +  # 使用极坐标绘制饼图
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_void() +  # 移除背景元素
  labs(fill = "Epitope Species", title = "Epitope Species in CD4T") +  # 添加标签
  theme(
    legend.key.size = unit(0.35, "cm"),  # 图例键的大小
    legend.text = element_text(size = 10),  # 图例文本大小
    legend.title = element_text(size = 12,face="bold"),  # 图例标题大小
    legend.spacing.x = unit(0.1, "cm"),  # 图例项之间的水平间距
    legend.position = 'right',  # 图例位置
    plot.margin = unit(c(1, 1, 1, 1), "char"),  # 图的外边距
    plot.title = element_text(size = 12, hjust = 0.5, vjust = 0, face = "bold")  # 总图标题样式
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))  # 图例列数为1

dev.off()

CD4T_data_COVID <- CD4T_data%>% filter(Epitope.species == "SARS-CoV-2")
# 生成数据表并排序
coverage_table <- as.data.frame(table(CD4T_data_COVID$ORF.Coverage))
colnames(coverage_table) <- c("ORF_Coverage", "Count")

# 计算每种覆盖率的百分比
coverage_table <- coverage_table %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  arrange(desc(Percentage))  # 按百分比从高到低排序
coverage_table$ORF_Coverage <- factor(coverage_table$ORF_Coverage,
                                      levels = unique(coverage_table$ORF_Coverage))

colors <- colorRampPalette(brewer.pal(9, "Set1"))(14)# 设置颜色（使用 RColorBrewer 中的调色板）
pdf("2_CD4T_COVID_Coverage_Pie_distance_2.pdf", 6, 4)

ggplot(coverage_table, aes(x = "", y = Percentage, fill = ORF_Coverage)) +
  geom_bar(width = 1, stat = "identity", alpha = 0.95, color = "white") +  # 设置透明度
  coord_polar(theta = "y") +  # 使用极坐标绘制饼图
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_void() +  # 移除背景元素
  labs(fill = "ORF Coverage", title = "ORF Coverage Distribution") +  # 添加标签
  theme(
    legend.key.size = unit(0.35, "cm"),  # 图例键的大小
    legend.text = element_text(size = 10),  # 图例文本大小
    legend.title = element_text(size = 12,face="bold"),  # 图例标题大小
    legend.spacing.x = unit(0.1, "cm"),  # 图例项之间的水平间距
    legend.position = 'right',  # 图例位置
    plot.margin = unit(c(1, 1, 1, 1), "char"),  # 图的外边距
    plot.title = element_text(size = 12, hjust = 0.5, vjust = 0, face = "bold")  # 总图标题样式
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))  # 图例列数为1

dev.off()



# 生成横向条形图，按 Count 对 ORF_Coverage 进行排序，并显示数字
pdf("2_CD4T_ORF_Coverage_Bar_Plot.pdf", 5, 4)  # 调整图形大小为6x4

ggplot(coverage_table, aes(y = reorder(ORF_Coverage, Count), x = Count, fill = ORF_Coverage)) +  # 对 ORF_Coverage 进行排序
  geom_bar(stat = "identity", alpha = 0.95, color = "white") +  # 绘制条形图
  geom_text(aes(label = Count), hjust = -0.2, color = "black", size = 3.5) +  # 在条形图右侧显示数字
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_minimal() +  # 使用简洁主题
  labs(fill = "ORF Coverage", title = "ORF Coverage in CD4T",
       y = "", x = "") +  # 去掉 X 轴标签
  theme(
    panel.background = element_blank(),  # 去掉背景
    panel.grid = element_blank(),  # 去除背景格线
    strip.background = element_blank(),  # 去掉标题的背景框
    axis.text.x = element_blank(),  # 去掉 X 轴文本
    axis.ticks.x = element_blank(),  # 去掉 X 轴刻度
    axis.text.y = element_text(size = 10, color='black'),  # 调整 Y 轴文本样式
    legend.position = 'none',  # 不显示图例
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    plot.title = element_text(size = 15, hjust = 0.5, vjust = 1, face = "bold")  # 总图标题样式
  ) +
  coord_cartesian(clip = "off")  # 允许文字显示在绘图区域外

dev.off()

# 修改的函数，显示每个 celltype 总数
plot_func <- function(data, var_name) {
  # 计算 count_data 为每个 celltype 在不同 ORF.Coverage 下的数量
  count_data <- table(data[[var_name]], data$ORF.Coverage)
  df_counts <- as.data.frame(count_data)
  colnames(df_counts) <- c("celltype", "ORF.Coverage", "count")
  
  # 计算每个 celltype 的总数，并将其添加到数据框
  total_counts <- aggregate(count ~ celltype, df_counts, sum)
  df_counts <- merge(df_counts, total_counts, by = "celltype", suffixes = c("", "_total"))
  
  # 绘制条形图，并在条形图右侧显示总数
  ggplot(df_counts, aes(x = reorder(celltype, count_total), y = count, fill = ORF.Coverage)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count_total, y = count_total), 
              hjust = -0.2, size = 1.2, color = "black") +  # 在右侧显示总数
    coord_flip() +
    theme_minimal() +
    xlab(var_name) +
    scale_fill_manual(values = colors) + 
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 6),
          strip.background = element_blank(),
          panel.background = element_blank(),  # 去掉背景
          panel.grid = element_blank(),  # 去除背景格线
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 4, color = 'black'),
          legend.title = element_text(size = 4,face='bold'),
          legend.text = element_text(size = 4),
          legend.position = 'right',
          legend.key.size = unit(0.08, "inch"),
          plot.margin = unit(c(1, 1, 1, 1), "char")
    ) 
}

# 将要生成图表的变量存入一个向量
variables <- c("ct_level_4", "Condition", "Donor")

# 生成PDF
pdf("2_CD4T_Covid_ORF_Spilt.pdf", 4,1.4)
CD4T_data_COVID$ORF.Coverage <- factor(CD4T_data_COVID$ORF.Coverage,
                                       level=c('Spike','ORF1ab','Nucleocapsid','ORF3','Matrix','ORF7a','RNP','NSP3','Envelope','ORF9b','ORF14'))
lapply(variables, plot_func, data = CD4T_data_COVID) # 使用lapply 循环生成图表

dev.off()

plot_func <- function(data, var_name) {
  count_data <- table(data[[var_name]], data$ORF.Coverage)
  df_counts <- as.data.frame(count_data)
  colnames(df_counts) <- c("celltype", "ORF.Coverage", "count")
  ggplot(df_counts, aes(x = reorder(celltype, count), y = count, fill = ORF.Coverage)) +
    geom_bar(stat = "identity", position = "fill") +
    coord_flip() +
    theme_minimal() +
    xlab(var_name) +
    scale_fill_manual(values =  colors) + 
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 6),
          strip.background = element_blank(),
          panel.background = element_blank(),  # 去掉背景
          panel.grid = element_blank(),  # 去除背景格线
          axis.text.x = element_text(size = 4,color='black'),
          axis.text.y = element_text(size = 4,color='black'),
          # axis.ticks.y = element_line(linewidth = 0.3), # 调整 x 轴刻度线的宽度
          # axis.ticks.x = element_line(linewidth = 0.3), # 调整 x 轴刻度线的宽度
          legend.title = element_text(size = 4),
          legend.text = element_text(size = 4),
          legend.position = 'right',
          legend.key.size = unit(0.08, "inch"),
          plot.margin = unit(c(1, 1, 1, 1), "char")
    )
}
pdf("2_CD4T_Covid_ORF_Spilt_2.pdf", 4,1.4)
CD4T_data_COVID$ORF.Coverage <- factor(CD4T_data_COVID$ORF.Coverage,
                                       level=c('Spike','ORF1ab','Nucleocapsid','ORF3','Matrix','ORF7a','RNP','NSP3','Envelope','ORF9b','ORF14'))
lapply(variables, plot_func, data = CD4T_data_COVID) # 使用lapply 循环生成图表

dev.off()

res_distance_2 <- process_data("data/3_CD8T_virus_tcr_distance_tcr_v_j_2.csv")
CD8T_data$tcr_v_j_cdr3 <- paste(paste(CD8T_data$tcr_v_gene, "01", sep = "*"),
                                paste(CD8T_data$tcr_j_gene, "01", sep = "*"),
                                CD8T_data$tcr_cdr3, sep = "_")
CD8T_data$ORF.Coverage <- res_distance_2$ORF.Coverage[match(CD8T_data$tcr_v_j_cdr3,res_distance_2$tcr_v_j_cdr3.x)]
CD8T_data$Epitope.species <- res_distance_2$Epitope.species[match(CD8T_data$tcr_v_j_cdr3,res_distance_2$tcr_v_j_cdr3.x)]
CD8T_data$Chain_distance <- res_distance_2$Distance_cdr3[match(CD8T_data$tcr_v_j_cdr3,res_distance_2$tcr_v_j_cdr3.x)]
table(CD8T_data$ORF.Coverage)
table(CD8T_data$Epitope.species)
table(CD8T_data$Chain_distance)
write.csv(CD8T_data,"CD8T_data_res_distance_2.csv")

species_table <- as.data.frame(table(CD8T_data$Epitope.species))
colnames(species_table) <- c("Epitope_species", "Count")

# 计算每种物种的百分比
species_table <- species_table %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  arrange(desc(Percentage))  # 按百分比从高到低排序

species_table$Epitope_species <- factor(species_table$Epitope_species,
                                        levels = unique(species_table$Epitope_species))

colors <- colorRampPalette(brewer.pal(12, "Paired"))(25)# 设置颜色（使用 RColorBrewer 中的调色板）

# 生成饼图
pdf("3_CD8T_Epitope_Species_Pie_distance_2.pdf", 6, 5)
ggplot(species_table, aes(x = "", y = Percentage, fill = Epitope_species)) +
  geom_bar(width = 0.9, stat = "identity", alpha = 0.95, color = "white") +  # 设置透明度
  coord_polar(theta = "y") +  # 使用极坐标绘制饼图
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_void() +  # 移除背景元素
  labs(fill = "Epitope Species", title = "Epitope Species in CD8T") +  # 添加标签
  theme(
    legend.key.size = unit(0.35, "cm"),  # 图例键的大小
    legend.text = element_text(size = 10),  # 图例文本大小
    legend.title = element_text(size = 12,face="bold"),  # 图例标题大小
    legend.spacing.x = unit(0.1, "cm"),  # 图例项之间的水平间距
    legend.position = 'right',  # 图例位置
    plot.margin = unit(c(1, 1, 1, 1), "char"),  # 图的外边距
    plot.title = element_text(size = 12, hjust = 0.5, vjust = 0, face = "bold")  # 总图标题样式
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))  # 图例列数为1

dev.off()

CD8T_data_COVID <- CD8T_data%>% filter(Epitope.species == "SARS-CoV-2")
# 生成数据表并排序
coverage_table <- as.data.frame(table(CD8T_data_COVID$ORF.Coverage))
colnames(coverage_table) <- c("ORF_Coverage", "Count")

# 计算每种覆盖率的百分比
coverage_table <- coverage_table %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  arrange(desc(Percentage))  # 按百分比从高到低排序
coverage_table$ORF_Coverage <- factor(coverage_table$ORF_Coverage,
                                      levels = unique(coverage_table$ORF_Coverage))

colors <- colorRampPalette(brewer.pal(9, "Set1"))(14)# 设置颜色（使用 RColorBrewer 中的调色板）
pdf("3_CD8T_COVID_Coverage_Pie_distance_2.pdf", 6, 4)

ggplot(coverage_table, aes(x = "", y = Percentage, fill = ORF_Coverage)) +
  geom_bar(width = 1, stat = "identity", alpha = 0.95, color = "white") +  # 设置透明度
  coord_polar(theta = "y") +  # 使用极坐标绘制饼图
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_void() +  # 移除背景元素
  labs(fill = "ORF Coverage", title = "ORF Coverage Distribution") +  # 添加标签
  theme(
    legend.key.size = unit(0.35, "cm"),  # 图例键的大小
    legend.text = element_text(size = 10),  # 图例文本大小
    legend.title = element_text(size = 12,face="bold"),  # 图例标题大小
    legend.spacing.x = unit(0.1, "cm"),  # 图例项之间的水平间距
    legend.position = 'right',  # 图例位置
    plot.margin = unit(c(1, 1, 1, 1), "char"),  # 图的外边距
    plot.title = element_text(size = 12, hjust = 0.5, vjust = 0, face = "bold")  # 总图标题样式
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))  # 图例列数为1

dev.off()


# 生成横向条形图，按 Count 对 ORF_Coverage 进行排序，并显示数字
pdf("3_CD8T_ORF_Coverage_Bar_Plot.pdf", 5, 4)  # 调整图形大小为6x4

ggplot(coverage_table, aes(y = reorder(ORF_Coverage, Count), x = Count, fill = ORF_Coverage)) +  # 对 ORF_Coverage 进行排序
  geom_bar(stat = "identity", alpha = 0.95, color = "white") +  # 绘制条形图
  geom_text(aes(label = Count), hjust = -0.2, color = "black", size = 3.5) +  # 在条形图右侧显示数字
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_minimal() +  # 使用简洁主题
  labs(fill = "ORF Coverage", title = "ORF Coverage in CD8T",
       y = "", x = "") +  # 去掉 X 轴标签
  theme(
    panel.background = element_blank(),  # 去掉背景
    panel.grid = element_blank(),  # 去除背景格线
    strip.background = element_blank(),  # 去掉标题的背景框
    axis.text.x = element_blank(),  # 去掉 X 轴文本
    axis.ticks.x = element_blank(),  # 去掉 X 轴刻度
    axis.text.y = element_text(size = 10, color='black'),  # 调整 Y 轴文本样式
    legend.position = 'none',  # 不显示图例
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    plot.title = element_text(size = 15, hjust = 0.5, vjust = 1, face = "bold")  # 总图标题样式
  ) +
  coord_cartesian(clip = "off")  # 允许文字显示在绘图区域外

dev.off()

# 定义函数
plot_func <- function(data, var_name) {
  count_data <- table(data[[var_name]], data$ORF.Coverage)
  df_counts <- as.data.frame(count_data)
  colnames(df_counts) <- c("celltype", "ORF.Coverage", "count")
  ggplot(df_counts, aes(x = reorder(celltype, count), y = count, fill = ORF.Coverage)) +
    geom_bar(stat = "identity") +#, position = "fill"
    coord_flip() +
    theme_minimal() +
    xlab(var_name) +
    scale_fill_manual(values =  colors) + 
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 6),
          strip.background = element_blank(),
          panel.background = element_blank(),  # 去掉背景
          panel.grid = element_blank(),  # 去除背景格线
          axis.text.x = element_text(size = 4,color='black'),
          axis.text.y = element_text(size = 4,color='black'),
          # axis.ticks.y = element_line(linewidth = 0.3), # 调整 x 轴刻度线的宽度
          # axis.ticks.x = element_line(linewidth = 0.3), # 调整 x 轴刻度线的宽度
          legend.title = element_text(size = 4),
          legend.text = element_text(size = 4),
          legend.position = 'right',
          legend.key.size = unit(0.08, "inch"),
          plot.margin = unit(c(1, 1, 1, 1), "char")
    )
}

# 将要生成图表的变量存入一个向量
variables <- c("ct_level_4", "Condition", "Donor")

# 生成PDF
pdf("3_CD8T_Covid_ORF_Spilt_2.pdf", 4,1.4)
CD8T_data_COVID$ORF.Coverage <- factor(CD8T_data_COVID$ORF.Coverage,
                                       level=c('Spike','ORF1ab','Nucleocapsid','ORF3','NSP3','Matrix','ORF9b','ORF8','ORF7a','Envelope'))
lapply(variables, plot_func, data = CD8T_data_COVID) # 使用lapply 循环生成图表

dev.off()

# 修改的函数，显示每个 celltype 总数
plot_func <- function(data, var_name) {
  # 计算 count_data 为每个 celltype 在不同 ORF.Coverage 下的数量
  count_data <- table(data[[var_name]], data$ORF.Coverage)
  df_counts <- as.data.frame(count_data)
  colnames(df_counts) <- c("celltype", "ORF.Coverage", "count")
  
  # 计算每个 celltype 的总数，并将其添加到数据框
  total_counts <- aggregate(count ~ celltype, df_counts, sum)
  df_counts <- merge(df_counts, total_counts, by = "celltype", suffixes = c("", "_total"))
  
  # 绘制条形图，并在条形图右侧显示总数
  ggplot(df_counts, aes(x = reorder(celltype, count_total), y = count, fill = ORF.Coverage)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count_total, y = count_total), 
              hjust = -0.2, size = 1.3, color = "black") +  # 在右侧显示总数
    coord_flip() +
    theme_minimal() +
    xlab(var_name) +
    scale_fill_manual(values = colors) + 
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 6),
          strip.background = element_blank(),
          panel.background = element_blank(),  # 去掉背景
          panel.grid = element_blank(),  # 去除背景格线
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 4, color = 'black'),
          legend.title = element_text(size = 4,face='bold'),
          legend.text = element_text(size = 4),
          legend.position = 'right',
          legend.key.size = unit(0.08, "inch"),
          plot.margin = unit(c(1, 1, 1, 1), "char")
    ) 
}

# 将要生成图表的变量存入一个向量
variables <- c("ct_level_4", "Condition", "Donor")

# 生成PDF
pdf("3_CD8T_Covid_ORF_Spilt_2.pdf", 4, 1.4)
CD8T_data_COVID$ORF.Coverage <- factor(CD8T_data_COVID$ORF.Coverage,
                                       levels = c('Spike', 'ORF1ab', 'Nucleocapsid', 'ORF3', 'NSP3', 'Matrix', 'ORF9b', 'ORF8', 'ORF7a', 'Envelope'))
lapply(variables, plot_func, data = CD8T_data_COVID)  # 使用lapply循环生成图表

dev.off()






