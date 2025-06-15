# 0.环境准备 --------------------------------------------------------------------
library("tidyverse")
library("ggpubr")
library("hrbrthemes")
library("gridExtra")
options(tibble.width = Inf, tibble.print_max = 10000)
library("Matrix")
library("doParallel")
library(BiocParallel)
library(Seurat)
library(future)
library(scales)
library(ggplot2)
library(data.table)

getwd()

# 文件路径列表
tcr_files <- c(
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX10tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX11tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX12tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX50tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX51tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX52tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX53tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX5tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX6tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX7tcr/outs/filtered_contig_annotations.csv',
  '/share/home/qlab/project_10x/3_data/project_cvv/sc5rMIX8tcr/outs/filtered_contig_annotations.csv'
)

# add function
add_clonotype <- function(file_dir, type) {
  file_dir <- gsub("/outs/filtered_contig_annotations.csv", "", file_dir) #字符串替换为空字符串
  mix_name <- gsub("/share/home/qlab/project_10x/3_data/project_cvv/", "", file_dir)
  mix_name <- gsub("sc5r", "", mix_name)
  mix_name <- gsub(type, "", mix_name)

  tcr <- read.csv(paste0(file_dir, "/outs/filtered_contig_annotations.csv"), sep = ",")
  tcr$barcode <- paste(mix_name, tcr$barcode, sep = "_") #合并为新列
  print(paste0(mix_name,
               "  unique:  ",length(unique(tcr$barcode)),
               "  all:  ",length(tcr$barcode),
               "  duplicated:  ",length(tcr$barcode)-length(unique(tcr$barcode)),
               "  pro:  ",length(unique(tcr$barcode))/length(tcr$barcode)
              ))
  tcr <- tcr[!duplicated(tcr$barcode), ] #筛选出不重复的barcode
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

  # Clonotype-centric info
  clono <- read.csv(paste0(file_dir, "/outs/clonotypes.csv"))
    
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono) 

  # Reorder so barcodes are first column and set them as rownames.
  tcr$clonotype_mix <- paste(mix_name, tcr$clonotype_id, sep = "_")
  colnames(tcr) <- paste(type, colnames(tcr), sep = "_")
  return(tcr)
}


# ---- add TCR metadata
tcr_metadata <- lapply(1:length(tcr_files), function(x) {
  add_clonotype(tcr_files[x], type = "tcr")
})
tcr_metadata <- rbindlist(tcr_metadata)
rownames_barcode <- gsub("/", "", tcr_metadata$tcr_barcode)
tcr_metadata$tcr_clonotype_mix <- gsub("/", "", tcr_metadata$tcr_clonotype_mix)
rownames(tcr_metadata) <- rownames_barcode

colnames(tcr_metadata)

data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3"
CD8T  <- readRDS(file.path(data_save_path, "4_Tcell_CD8T/2_Reumap/Trajectory_CD8T.rds"))
CD4T <- readRDS(file.path(data_save_path, "4_Tcell_CD4T/3_ReCCA_with_Treg/CD4T_ReCCA_with_Treg.rds"))

# ---- add TCR clonotype, rename tcr_clonotype_id based on cdr3s_aa
CD4T <- AddMetaData(object = CD4T, metadata = tcr_metadata)
table(!is.na(CD4T$tcr_clonotype_id))
change_id <- data.frame(raw = na.omit(unique(CD4T$tcr_cdr3s_aa)))
change_id$new <- paste0("clonotype", 1:nrow(change_id))
CD4T$tcr_clonotype_id_new <- change_id$new[match(CD4T$tcr_cdr3s_aa, change_id$raw)]
table(CD4T$ct_level_4, CD4T$tcr_clonotype_id_new > 0)
head(change_id)

# ---- add TCR clonotype, rename tcr_clonotype_id based on cdr3s_aa
CD8T <- AddMetaData(object = CD8T, metadata = tcr_metadata)
table(!is.na(CD8T$tcr_clonotype_id))
change_id <- data.frame(raw = na.omit(unique(CD8T$tcr_cdr3s_aa)))
change_id$new <- paste0("clonotype", 1:nrow(change_id))
CD8T$tcr_clonotype_id_new <- change_id$new[match(CD8T$tcr_cdr3s_aa, change_id$raw)]
table(CD8T$ct_level_4, CD8T$tcr_clonotype_id_new > 0)
head(change_id)

nrow(tcr_metadata)
nrow(CD4T@meta.data)
length(unique(CD4T$tcr_clonotype_id_new))
length(unique(CD8T$tcr_clonotype_id_new))

saveRDS(tcr_metadata, "data/1_tcr_metadata.rds")
write.csv(tcr_metadata, "data/1_tcr_metadata.csv")
saveRDS(CD4T,  "data/1_CD4T_TCR.rds")
saveRDS(CD8T,  "data/1_CD8T_TCR.rds")

CD4T <- readRDS("data/1_CD4T_TCR.rds")
CD8T <- readRDS("data/1_CD8T_TCR.rds")

# 定义函数来处理细胞对象的克隆型数据
nclos_clonotype <-  c("lightgrey", "darkgrey", "#E0C5BEFF",
                      "#F39B7FB2", "#E57E7EFF", "#DE5C00FF", "#B22C2CFF", "#573333FF")
process_clonotype_data <- function(cell_type_obj,ct_seq) {
  # 设置 ct_level_4 因子顺序
  cell_type_obj$ct_level_4 <- factor(
    cell_type_obj$ct_level_4,
    levels = ct_seq
  )
  
  # 设置 Donor 因子顺序
  cell_type_obj$Donor <- factor(
    cell_type_obj$Donor,
    levels = c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5", "Donor6", "Donor7", "Donor8",
               "Donor9", "Donor10", "Donor11", "Donor12")
  )
  
  # 设置 Batch 因子顺序
  cell_type_obj$Batch <- factor(
    cell_type_obj$Batch,
    levels = c("MIX5", "MIX6", "MIX7", "MIX8", "MIX10", "MIX11", "MIX12", "MIX50", "MIX51",
               "MIX52", "MIX53")
  )
  
  # 初始化 Clonotype_num_ld 列
  cell_type_obj$Clonotype_num_ld <- "Not detected"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num == 1] <- "n = 1"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 1 & cell_type_obj$Clonotype_num <= 5] <- "1 < n <= 5"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 5 & cell_type_obj$Clonotype_num <= 10] <- "5 < n <= 10"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 10 & cell_type_obj$Clonotype_num <= 20] <- "10 < n <= 20"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 20 & cell_type_obj$Clonotype_num <= 50] <- "20 < n <= 50"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 50 & cell_type_obj$Clonotype_num <= 100] <- "50 < n <= 100"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 100] <- "n > 100"
  
  # 将 Clonotype_num_ld 转换为因子，并设定顺序
  cell_type_obj$Clonotype_num_ld <- factor(
    cell_type_obj$Clonotype_num_ld,
    levels = c("Not detected", "n = 1", "1 < n <= 5", "5 < n <= 10", "10 < n <= 20", 
               "20 < n <= 50", "50 < n <= 100", "n > 100")
  )
  
  return(cell_type_obj)
}


# 定义函数来处理细胞对象的克隆型数据
process_clonotype_data <- function(cell_type_obj,ct_seq) {
  # 设置 ct_level_4 因子顺序
  cell_type_obj$ct_level_4 <- factor(
    cell_type_obj$ct_level_4,
    levels = ct_seq
  )
  
  # 设置 Donor 因子顺序
  cell_type_obj$Donor <- factor(
    cell_type_obj$Donor,
    levels = c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5", "Donor6", "Donor7", "Donor8",
               "Donor9", "Donor10", "Donor11", "Donor12")
  )
  
  # 设置 Batch 因子顺序
  cell_type_obj$Batch <- factor(
    cell_type_obj$Batch,
    levels = c("MIX5", "MIX6", "MIX7", "MIX8", "MIX10", "MIX11", "MIX12", "MIX50", "MIX51",
               "MIX52", "MIX53")
  )
  
  # 初始化 Clonotype_num_ld 列
  cell_type_obj$Clonotype_num_ld <- "Not detected"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num == 1] <- "n = 1"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 1 & cell_type_obj$Clonotype_num <= 3] <- "1 < n <= 3"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 3 & cell_type_obj$Clonotype_num <= 5] <- "3 < n <= 5"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 5] <- "n > 5"
  
  # 将 Clonotype_num_ld 转换为因子，并设定顺序
  cell_type_obj$Clonotype_num_ld <- factor(
    cell_type_obj$Clonotype_num_ld,
    levels = c("Not detected", "n = 1", "1 < n <= 3", "3 < n <= 5",  "n > 5")
  )
  
  return(cell_type_obj)
}


# 定义通用绘图函数
plot_clonotype_umap <- function(cell_type_obj, file_name, nclos_clonotype,w,h,pt) {
  pdf(file_name, w,h)
  print(
    DimPlot(
      object = cell_type_obj, reduction = "umap", label = FALSE, pt.size = pt, order = TRUE,
      group.by = "Clonotype_num_ld", cols = nclos_clonotype
    )
  )
  dev.off()
}
# 定义通用绘图函数
plot_dim <- function(data, split_by) {
  print(
    DimPlot(
      object = data, 
      reduction = "umap", 
      label = FALSE, 
      pt.size = 0.01, 
      order = TRUE,
      group.by = "Clonotype_num_ld", 
      split.by = split_by,
      cols = nclos_clonotype
    )
  )
}


unique(CD4T$Condition)
unique(CD8T$Condition)

# CD4T 调用：
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
ct_seq <- c("CD4_Tn", "CD4_Tcm", "CD4_Th2", "CD4_T_CCR6", "CD4_Tem", "CD4_Temra", "CD4_Treg")
clonotype_num <- as.data.frame(table(CD4T$tcr_clonotype_id_new))
CD4T$Clonotype_num <- clonotype_num$Freq[match(CD4T$tcr_clonotype_id_new, clonotype_num$Var1)]
CD4T <- process_clonotype_data(CD4T,ct_seq)

plot_clonotype_umap(CD4T, "2_CD4T_clonotype_umap.pdf", nclos_clonotype,5,3.5,0.0001)

pdf("2_CD4T_clonotype_umap_split.pdf", 18, 4)
split_conditions <- c("Clonotype_num_ld", "Condition", "ct_level_4", "Donor", "Batch")
lapply(split_conditions, function(split_by) plot_dim(CD4T, split_by))
dev.off()

# CD8T 调用：
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
ct_seq <- c("CD8_Tn", "CD8_Tem", "CD8_Temra", "CD8_NELL2")
clonotype_num <- as.data.frame(table(CD8T$tcr_clonotype_id_new))
CD8T$Clonotype_num <- clonotype_num$Freq[match(CD8T$tcr_clonotype_id_new, clonotype_num$Var1)]
CD8T <- process_clonotype_data(CD8T,ct_seq)

plot_clonotype_umap(CD8T, "2_CD8T_clonotype_umap.pdf", nclos_clonotype,5,3.5,0.0001)

pdf("2_CD8T_clonotype_umap_split.pdf", 18, 4)
split_conditions <- c("Clonotype_num_ld", "Condition", "ct_level_4", "Donor", "Batch")
lapply(split_conditions, function(split_by) plot_dim(CD8T, split_by))
dev.off()



CD4T <- readRDS("data/1_CD4T_TCR.rds")

CD8T <- readRDS("data/1_CD8T_TCR.rds")

clonotype_num <- as.data.frame(table(CD8T$tcr_clonotype_id_new))
CD8T$Clonotype_num <- clonotype_num$Freq[match(CD8T$tcr_clonotype_id_new, clonotype_num$Var1)]

df_l1 <- CD8T@meta.data %>%
  filter(!is.na(Pseudotime1)) %>%
  rownames_to_column(var = "Cell") %>%
  mutate(lineage = "Lineage1_CD8_Temra")
df_l2 <- CD8T@meta.data %>%
  filter(!is.na(Pseudotime2)) %>%
  rownames_to_column(var = "Cell") %>%
  mutate(lineage = "Lineage2_CD8_NELL2")
df <- rbind(df_l1, df_l2)
head(df)

# ---- lineage
p <- ggplot(df, aes(x = Pseudotime, y = Clonotype_num, color = lineage)) +
  geom_smooth(method = "gam", se = F) +
  scale_color_manual(values = c("Lineage1_CD8_Temra" = "#9b3a74", 
                                "Lineage2_CD8_NELL2" = "#4da0a0")) +
  ylab("Number of Clonotypes") +
  theme_bw()
ggsave("3_CD8T_density_all.pdf", p, width = 6, height = 4)

# ---- upset plot
library(UpSetR)
clonotypes <- data.frame(
  id = unique(CD4T$tcr_clonotype_id_new),
  CD4_Tn = 0,
  CD4_Tcm = 0,
     CD4_Th2 = 0,
  CD4_T_CCR6 = 0,
  CD4_Tem = 0,
  CD4_Temra = 0,
    CD4_Treg = 0
)
clonotypes <- clonotypes[!is.na(clonotypes$id), ]
for (i in c(1:ncol(CD4T))) {
  #print(i)
  id <- CD4T$tcr_clonotype_id_new[i]
  if (is.na(id)) {
    next
  }
  CellType <- as.character(CD4T$ct_level_4[i])
  clonotypes[clonotypes$id == id, CellType] <- clonotypes[clonotypes$id == id, CellType] + 1
}
clonotypes[clonotypes > 0] <- 1
head(clonotypes)
print(i)

pdf("4_CD4T_upset_ct_level_4.pdf", width = 8, height = 6)
upset(clonotypes,sets = c( "CD4_Tn", "CD4_Tcm",
                          "CD4_Th2", "CD4_T_CCR6","CD4_Tem",
                           "CD4_Temra","CD4_Treg"
                          ),
            keep.order = TRUE, text.scale = 1,
            mb.ratio = c(0.7, 0.3),
            matrix.color = "#B22C2CFF", main.bar.color = "#B22C2CFF",
           shade.color = "#F3B1A0",sets.bar.color = "#B22C2CFF", 
            point.size = 2, line.size = 0.8)
dev.off()

# ---- upset plot
library(UpSetR)
clonotypes <- data.frame(
  id = unique(CD4T$tcr_clonotype_id_new),
  Donor1 = 0,
  Donor2 = 0,
  Donor3 = 0,
  Donor4 = 0,
  Donor5 = 0,
  Donor6 = 0,
  Donor7 = 0,
  Donor8 = 0,
  Donor9 = 0,
  Donor10 = 0,
  Donor11 = 0,
  Donor12 = 0
)
clonotypes <- clonotypes[!is.na(clonotypes$id), ]

for (i in c(1:ncol(CD4T))) {
  #print(i)
  id <- CD4T$tcr_clonotype_id_new[i]
  if (is.na(id)) {
    next
  }
  donor_n <- as.character(CD4T$Donor[i])
  clonotypes[clonotypes$id == id, donor_n] <- clonotypes[clonotypes$id == id, donor_n] + 1
}
clonotypes[clonotypes > 0] <- 1

pdf("4_CD4T_upset_donor.pdf", width = 10, height = 8)
upset(clonotypes,sets = c("Donor1" ,"Donor2","Donor3" ,"Donor4","Donor5" ,"Donor6","Donor7" ,"Donor8",
                          "Donor9" ,"Donor10","Donor11" ,"Donor12"),
            keep.order = TRUE, text.scale = 1,
            matrix.color = "#B22C2CFF", main.bar.color = "#B22C2CFF",
           shade.color = "#F3B1A0",sets.bar.color = "#B22C2CFF", 
            mb.ratio = c(0.6, 0.4),
            point.size = 2, line.size = 1)
dev.off()

# ---- upset plot
library(UpSetR)
clonotypes <- data.frame(
  id = unique(CD8T$tcr_clonotype_id_new),
  CD8_Tn = 0,
  CD8_Temra = 0,
  CD8_NELL2 = 0,
  CD8_Tem = 0
)
clonotypes <- clonotypes[!is.na(clonotypes$id), ]
for (i in c(1:ncol(CD8T))) {
  #print(i)
  id <- CD8T$tcr_clonotype_id_new[i]
  if (is.na(id)) {
    next
  }
  CellType <- as.character(CD8T$ct_level_4[i])
  clonotypes[clonotypes$id == id, CellType] <- clonotypes[clonotypes$id == id, CellType] + 1
}
clonotypes[clonotypes > 0] <- 1
print(i)
head(clonotypes)

pdf("4_CD8T_upset_ct_level_4.pdf", width = 6, height = 4.5)
upset(clonotypes,sets = c("CD8_Tn","CD8_NELL2", "CD8_Tem","CD8_Temra"),
            keep.order = TRUE, text.scale = 1,
            mb.ratio = c(0.75, 0.25),
            matrix.color = "#B22C2CFF", main.bar.color = "#B22C2CFF",
           shade.color = "#F3B1A0",sets.bar.color = "#B22C2CFF", 
            point.size = 2, line.size = 0.6)
dev.off()

# ---- upset plot
library(UpSetR)
clonotypes <- data.frame(
  id = unique(CD8T$tcr_clonotype_id_new),
  Donor1 = 0,
  Donor2 = 0,
  Donor3 = 0,
  Donor4 = 0,
  Donor5 = 0,
  Donor6 = 0,
  Donor7 = 0,
  Donor8 = 0,
  Donor9 = 0,
  Donor10 = 0,
  Donor11 = 0,
  Donor12 = 0
)
clonotypes <- clonotypes[!is.na(clonotypes$id), ]
for (i in c(1:ncol(CD8T))) {
  #print(i)
  id <- CD8T$tcr_clonotype_id_new[i]
  if (is.na(id)) {
    next
  }
  donor_n <- as.character(CD8T$Donor[i])
  clonotypes[clonotypes$id == id, donor_n] <- clonotypes[clonotypes$id == id, donor_n] + 1
}
clonotypes[clonotypes > 0] <- 1
print(i)
head(clonotypes)

pdf("4_CD8T_upset_donor.pdf", width = 10, height = 8)
upset(clonotypes,sets = c("Donor1" ,"Donor2","Donor3" ,"Donor4","Donor5" ,"Donor6","Donor7" ,"Donor8",
                          "Donor9" ,"Donor10","Donor11" ,"Donor12"),
            keep.order = TRUE, text.scale = 1,
            mb.ratio = c(0.6, 0.4),
            matrix.color = "#B22C2CFF", main.bar.color = "#B22C2CFF",
           shade.color = "#F3B1A0",sets.bar.color = "#B22C2CFF", 
            point.size = 2, line.size = 1)
dev.off()

CD8T <- readRDS("data/1_CD8T_TCR.rds")
#colnames(CD4T@meta.data)
table(CD8T$ct_level_4)

#估算下Treg细胞的TCR richness在每个节点的变化，richness = unique（TCR）/all TCR
library(dplyr)

tem_temra_df <- subset(CD8T, subset = ct_level_4 %in% c("CD8_Tem", "CD8_Temra"))
# 提取相关元数据
richness_df <- as.data.frame(tem_temra_df@meta.data[ ,c("tcr_clonotype_id_new","Condition","Donor")]) %>%
  na.omit()
# 计算 unique TCR 数量
richness_unique <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())

# 计算总的 TCR 数量
richness_total <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])

# 合并并计算 richness
richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(richness$Var2)
richness$richness <- richness$unique_TCR / richness$total_TCR
# 整理列名
names(richness) <- c("Condition", "Donor", "unique_TCR", "total_TCR", "richness")
# 按 Condition 分组计算均值/中位数
richness$Condition <- factor(richness$Condition, level = unique(richness$Condition))
richness <- as.data.frame(richness) %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),
         median_richness = median(richness, na.rm = TRUE))
# 查看结果
head(richness)
# 导出为 CSV
write.csv(richness, "5_richness_CD8_Tem_Temra.csv", row.names = FALSE)


library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

pdf("5_CD8T_Temra_Tem_line_richness.pdf", width = 3, height = 2.2)

# 去除 Condition 为 "C3" 的样本
filtered_richness <- richness 
# 计算平均值和标准误
summary_df <- filtered_richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

# 绘图
ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#53A85F", size = 0.5) +
  geom_point(color = "#53A85F", size = 1.8) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.15, linetype = "dashed", alpha = 0.6, color = "#53A85F") +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  theme_bw(base_size = 10) +
  labs(title = "CD8+ Temra and Tem Richness",
       y = "Richness", x = "Condition")

dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

# 1. Mean ± SEM 折线图
summary_df <- richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

p1 <- ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#2c7fb8", size = 1) +
  geom_point(color = "#2c7fb8", size = 2) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.2, color = "#2c7fb8") +
  theme_bw() +
  labs(title = "CD8T_Temra_Tem TCR Richness Mean ± SEM",
       y = "Richness", x = "Condition")

# 2. 热图处理：保存为 grob
richness_mat <- richness %>%
  select(Condition, Donor, richness) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = richness) %>%
  tibble::column_to_rownames("Donor")

pheat <- pheatmap(as.matrix(richness_mat),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  color = viridis::viridis(100),
                  main = "CD8T_Temra_Tem TCR Richness Heatmap",
                  silent = TRUE)  # 返回 grob 对象

# 3. 每个供体的轨迹图
richness$Condition <- factor(richness$Condition, levels = unique(richness$Condition))  # 确保顺序
p3 <- ggplot(richness, aes(x = Condition, y = richness, group = Donor, color = Donor)) +
  geom_line(alpha = 0.6, size = 0.8) +
  geom_point(size = 1.2) +
  theme_bw() +
  scale_color_viridis_d(option = "D") +
  labs(title = "CD8T_Temra_Tem TCR Richness Dynamics per Donor",
       y = "Richness", x = "Condition") +
  theme(legend.position = "right")

# 4. 合并三图
pdf("5_CD8T_Temra_Tem_richness_combined.pdf", width = 6, height = 10)
grid.arrange(
  p1,
  pheat[[4]],  # pheatmap 的 grob
  p3,
  ncol = 1,
  heights = c(1, 1, 1)  # 可微调每幅图的高度比例
)
dev.off()





CD4T <- readRDS("data/1_CD4T_TCR.rds")
#colnames(CD4T@meta.data)
table(CD4T$ct_level_4)

#估算下Treg细胞的TCR richness在每个节点的变化，richness = unique（TCR）/all TCR
library(dplyr)

# 筛选 CD4_Treg 子集
treg_df <- subset(CD4T, subset = ct_level_4 == "CD4_Treg")

# 提取相关元数据
richness_df <- as.data.frame(treg_df@meta.data[ ,c("tcr_clonotype_id_new","Condition","Donor")]) %>%
  na.omit()

# 计算 unique TCR 数量
richness_unique <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())

# 计算总的 TCR 数量
richness_total <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])

# 合并并计算 richness
richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(richness$Var2)
richness$richness <- richness$unique_TCR / richness$total_TCR

# 整理列名
names(richness) <- c("Condition", "Donor", "unique_TCR", "total_TCR", "richness")
# 按 Condition 分组计算均值/中位数
richness$Condition <- factor(richness$Condition, level = unique(richness$Condition))
richness <- as.data.frame(richness) %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),
         median_richness = median(richness, na.rm = TRUE))

# 查看结果
head(richness)

# 导出为 CSV
write.csv(richness, "5_richness_CD4_Treg.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)
library(viridis)

# Donor 排序（根据实际供体名称可调整）
richness$Donor <- factor(richness$Donor,
                         levels = c("Donor1" ,"Donor2","Donor3" ,"Donor4","Donor5" ,"Donor6",
                                    "Donor7" ,"Donor8","Donor9" ,"Donor10","Donor11" ,"Donor12"))
# --- 图1：剔除离群点（距离均值最远） ---
pdf("5_Treg_richness_Condition_filtered.pdf", width = 4, height = 2.3)

richness_filtered <- richness %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),
         diff = abs(richness - mean_richness)) %>%
  filter(diff != max(diff, na.rm = TRUE)) %>%  # 去掉每组中离均值最远的一个点
  ungroup()
p1 <- ggplot(richness_filtered, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2,
             position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +
  theme_bw(base_size = 10) +
  labs(title = "     Richness in CD4_Treg (filtered)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +
  ylab("Richness")
print(p1)
dev.off()

# --- 图2：全数据保留 ---
pdf("5_Treg_richness_Condition_all.pdf", width = 4, height = 2.3)

p2 <- ggplot(richness, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2,
             position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +
  theme_bw(base_size = 10) +
  labs(title = "     Richness in CD4_Treg (all)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +
  ylab("Richness")

print(p2)
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

# 1. Mean ± SEM 折线图
summary_df <- richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

p1 <- ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#2c7fb8", size = 1) +
  geom_point(color = "#2c7fb8", size = 2) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.2, color = "#2c7fb8") +
  theme_bw() +
  labs(title = "Treg TCR Richness Mean ± SEM",
       y = "Richness", x = "Condition")
# 2. 热图处理：保存为 grob
richness_mat <- richness %>%
  select(Condition, Donor, richness) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = richness) %>%
  tibble::column_to_rownames("Donor")

pheat <- pheatmap(as.matrix(richness_mat),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  color = viridis::viridis(100),
                  main = "CD4_Treg TCR Richness Heatmap",
                  silent = TRUE)  # 返回 grob 对象

# 3. 每个供体的轨迹图
richness$Condition <- factor(richness$Condition, levels = unique(richness$Condition))  # 确保顺序
p3 <- ggplot(richness, aes(x = Condition, y = richness, group = Donor, color = Donor)) +
  geom_line(alpha = 0.6, size = 0.8) +
  geom_point(size = 1.2) +
  theme_bw() +
  scale_color_viridis_d(option = "D") +
  labs(title = "Treg TCR Richness Dynamics per Donor",
       y = "Richness", x = "Condition") +
  theme(legend.position = "right")

# 4. 合并三图
pdf("5_Treg_richness_combined.pdf", width = 6, height = 10)
grid.arrange(
  p1,
  pheat[[4]],  # pheatmap 的 grob
  p3,
  ncol = 1,
  heights = c(1, 1, 1)  # 可微调每幅图的高度比例
)
dev.off()



CD4T <- readRDS("data/1_CD4T_TCR.rds")

#colnames(CD4T@meta.data)
table(CD4T$ct_level_4)

#估算下Treg细胞的TCR richness在每个节点的变化，richness = unique（TCR）/all TCR
library(dplyr)

tem_temra_df <- subset(CD4T, subset = ct_level_4 %in% c("CD4_Tem", "CD4_Temra"))
# 提取相关元数据
richness_df <- as.data.frame(tem_temra_df@meta.data[ ,c("tcr_clonotype_id_new","Condition","Donor")]) %>%
  na.omit()
# 计算 unique TCR 数量
richness_unique <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())

# 计算总的 TCR 数量
richness_total <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])

# 合并并计算 richness
richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(richness$Var2)
richness$richness <- richness$unique_TCR / richness$total_TCR

# 整理列名
names(richness) <- c("Condition", "Donor", "unique_TCR", "total_TCR", "richness")
# 按 Condition 分组计算均值/中位数
richness$Condition <- factor(richness$Condition, level = unique(richness$Condition))
richness <- as.data.frame(richness) %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),
         median_richness = median(richness, na.rm = TRUE))
# 查看结果
head(richness)
# 导出为 CSV
write.csv(richness, "5_richness_CD4_Tem_Temra.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)
library(viridis)

# Donor 排序（根据实际供体名称可调整）
richness$Donor <- factor(richness$Donor,
                         levels = c("Donor1" ,"Donor2","Donor3" ,"Donor4","Donor5" ,"Donor6",
                                    "Donor7" ,"Donor8","Donor9" ,"Donor10","Donor11" ,"Donor12"))
# --- 图1：剔除离群点（距离均值最远） ---
pdf("5_CD4T_Temra_Tem_richness_Condition_filtered.pdf", width = 4, height = 2.3)

richness_filtered <- richness %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),
         diff = abs(richness - mean_richness)) %>%
  filter(diff != max(diff, na.rm = TRUE)) %>%  # 去掉每组中离均值最远的一个点
  ungroup()
p1 <- ggplot(richness_filtered, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2,
             position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +
  theme_bw(base_size = 10) +
  labs(title = "    Richness in CD4T_Temra_Tem (filtered)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +
  ylab("Richness")
print(p1)
dev.off()


# --- 图2：全数据保留 ---
pdf("5_CD4T_Temra_Tem_richness_Condition_all.pdf", width = 4, height = 2.3)

p2 <- ggplot(richness, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2,
             position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +
  theme_bw(base_size = 10) +
  labs(title = "    Richness in CD4T_Temra_Tem (all)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +
  ylab("Richness")

print(p2)
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

# 1. Mean ± SEM 折线图
summary_df <- richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

p1 <- ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#2c7fb8", size = 1) +
  geom_point(color = "#2c7fb8", size = 2) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.2, color = "#2c7fb8") +
  theme_bw() +
  labs(title = "CD4T_Temra_Tem TCR Richness Mean ± SEM",
       y = "Richness", x = "Condition")

# 2. 热图处理：保存为 grob
richness_mat <- richness %>%
  select(Condition, Donor, richness) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = richness) %>%
  tibble::column_to_rownames("Donor")

pheat <- pheatmap(as.matrix(richness_mat),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  color = viridis::viridis(100),
                  main = "CD4T_Temra_Tem TCR Richness Heatmap",
                  silent = TRUE)  # 返回 grob 对象

# 3. 每个供体的轨迹图
richness$Condition <- factor(richness$Condition, levels = unique(richness$Condition))  # 确保顺序
p3 <- ggplot(richness, aes(x = Condition, y = richness, group = Donor, color = Donor)) +
  geom_line(alpha = 0.6, size = 0.8) +
  geom_point(size = 1.2) +
  theme_bw() +
  scale_color_viridis_d(option = "D") +
  labs(title = "CD4T_Temra_Tem TCR Richness Dynamics per Donor",
       y = "Richness", x = "Condition") +
  theme(legend.position = "right")

# 4. 合并三图
pdf("5_CD4T_Temra_Tem_richness_combined.pdf", width = 6, height = 10)
grid.arrange(
  p1,
  pheat[[4]],  # pheatmap 的 grob
  p3,
  ncol = 1,
  heights = c(1, 1, 1)  # 可微调每幅图的高度比例
)
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

pdf("5_CD4T_Temra_Tem_line_richness.pdf", width = 3, height = 2.2)


filtered_richness <- richness 

# 计算平均值和标准误
summary_df <- filtered_richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

# 绘图
ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#E95C39", size = 0.5) +
  geom_point(color = "#E95C39",size = 1.8) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.15, linetype = "dashed", alpha = 0.6, color =  "#E95C39") +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  theme_bw(base_size = 10) +
  labs(title = "CD4+ Temra and Tem Richness",
       y = "Richness", x = "Condition")

dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

pdf("5_CD4T_Temra_Tem_line_richness.pdf", width = 3, height = 2.2)
# 1. Mean ± SEM 折线图
summary_df <- richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#E95C39", size = 0.5)+
  geom_point(color = "#E95C39", size = 1.8) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.15, linetype = "dashed", alpha = 0.6,color = "#E95C39") +
  theme_bw(base_size = 10) +
  labs(title = "CD4+ Temra and Tem Richness",
       y = "Richness", x = "Condition")

dev.off()

richness_df <- as.data.frame(CD4T@meta.data[ ,c("tcr_clonotype_id_new","Condition","Donor")])%>% 
                  na.omit()
# unique
richness_unique <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())
# total
richness_total <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])
# richness 
richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(
  richness$Var2
)
richness$richness <- richness$unique_TCR / richness$total_TCR
#tail(richness)
names(richness) <- c("Condition", "Donor", "unique_TCR", "total_TCR", "richness")
unique(richness$Condition)
richness$Condition <- factor(richness$Condition,level=unique(richness$Condition))
richness <- as.data.frame(richness) %>% 
            group_by(Condition) %>% 
            mutate(mean_richness = mean(richness, na.rm = TRUE),median_richness=median(richness, na.rm = TRUE))
head(richness)
write.csv(richness, "5_richness_CD4T.csv")

richness <- read.csv("5_richness_CD4T.csv")
head(richness)

richness$Donor <- factor(richness$Donor,
                        level=c("Donor1" ,"Donor2","Donor3" ,"Donor4","Donor5" ,"Donor6","Donor7" ,"Donor8",
                          "Donor9" ,"Donor10","Donor11" ,"Donor12"))
pdf("5_CD4T_richness_Condition.pdf", width = 4, height = 2.3)
richness_filtered <- richness %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),  # 计算均值
         diff = abs(richness - mean_richness)) %>%  # 计算与均值的差异
  filter(diff != max(diff, na.rm = TRUE)) %>%  # 去掉距离均值差异最大的数据点
  ungroup()
p2 <- ggplot(richness_filtered, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2, position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +  # 使用 viridis 包中的渐变色
  theme_bw(base_size = 10) +
   labs(title = "     Richness in CD4T") +
  theme(
    panel.grid.major = element_blank(),  # 去掉主网格
    panel.grid.minor = element_blank(),  # 去掉次网格
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +ylab("Richness")
# 添加虚线和点到 p2 图上
# overlay_plot <- p2 +
#   geom_point(data = richness_filtered, aes(x = Condition, y = mean_richness, group = 1), alpha = 0.95, size = 2, color = "#B22C2CFF") +
#   geom_line(data = richness_filtered, aes(x = Condition, y = mean_richness, group = 1), alpha = 0.95, linewidth = 0.8, color = "#B22C2CFF")
#print(overlay_plot)
print(p2)
dev.off()

pdf("5_CD4T_richness_Condition_all.pdf", width = 4, height = 2.3)
p2 <- ggplot(richness, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2, position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +  # 使用 viridis 包中的渐变色
  theme_bw(base_size = 10) + # 去掉背景格子
  labs(title = "     Richness in CD4T(all)") +
  theme(
    panel.grid.major = element_blank(),  # 去掉主网格
    panel.grid.minor = element_blank(),  # 去掉次网格
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +ylab("Richness")
#overlay_plot <- p2 +geom_point(data = richness, aes(x = Condition, y = median_richness, group = 1), alpha = 0.95, size = 2, color = "#B22C2CFF") +
  #  geom_line(data = richness, aes(x = Condition, y = median_richness, group = 1), alpha = 0.95, linewidth = 0.8, color = "#B22C2CFF")
#print(overlay_plot)
print(p2)
# 关闭PDF设备
dev.off()

richness_df <- as.data.frame(CD8T@meta.data[ ,c("tcr_clonotype_id_new","Condition","Donor")])%>% 
                  na.omit()
# unique
richness_unique <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())
# total
richness_total <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])
# richness 
richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(
  richness$Var2
)
richness$richness <- richness$unique_TCR / richness$total_TCR
#tail(richness)
names(richness) <- c("Condition", "Donor", "unique_TCR", "total_TCR", "richness")
unique(richness$Condition)
richness$Condition <- factor(richness$Condition,level=unique(richness$Condition))
richness <- as.data.frame(richness) %>% 
            group_by(Condition) %>% 
            mutate(mean_richness = mean(richness, na.rm = TRUE),median_richness=median(richness, na.rm = TRUE))
head(richness)
write.csv(richness, "5_richness_CD8T.csv")

richness <- read.csv("5_richness_CD8T.csv")
head(richness)

richness$Donor <- factor(richness$Donor,
                        level=c("Donor1" ,"Donor2","Donor3" ,"Donor4","Donor5" ,"Donor6","Donor7" ,"Donor8",
                          "Donor9" ,"Donor10","Donor11" ,"Donor12"))
pdf("5_CD8T_richness_Condition.pdf", width = 4, height = 2.3)
richness_filtered <- richness %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),  # 计算均值
         diff = abs(richness - mean_richness)) %>%  # 计算与均值的差异
  filter(diff != max(diff, na.rm = TRUE)) %>%  # 去掉距离均值差异最大的数据点
  ungroup()
p2 <- ggplot(richness_filtered, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2, position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +  # 使用 viridis 包中的渐变色
  theme_bw(base_size = 10) +
  labs(title = "       Richness in CD8T") +
  theme(
    panel.grid.major = element_blank(),  # 去掉主网格
    panel.grid.minor = element_blank(),  # 去掉次网格
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +ylab("Richness")
# 添加虚线和点到 p2 图上
# overlay_plot <- p2 +
#   geom_point(data = richness_filtered, aes(x = Condition, y = mean_richness, group = 1), alpha = 0.95, size = 2, color = "#B22C2CFF") +
#   geom_line(data = richness_filtered, aes(x = Condition, y = mean_richness, group = 1), alpha = 0.95, linewidth = 0.8, color = "#B22C2CFF")
#print(overlay_plot)
print(p2)
dev.off()

pdf("5_CD8T_richness_Condition_all.pdf", width = 4, height = 2.3)
p2 <- ggplot(richness, aes(x = Condition, y = richness)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 2, position = position_jitter(width = 0.25)) +
  scale_fill_viridis_d(option = "D") +  # 使用 viridis 包中的渐变色
  theme_bw(base_size = 10) + # 去掉背景格子
  labs(title = "     Richness in CD8T(all)") +
  theme(
    panel.grid.major = element_blank(),  # 去掉主网格
    panel.grid.minor = element_blank(),  # 去掉次网格
    axis.text.x = element_text(angle = 0),
    legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    legend.key.size = unit(0.15, "inch"),
    legend.text = element_text(size = 8)
  ) +ylab("Richness")
#overlay_plot <- p2 +geom_point(data = richness, aes(x = Condition, y = median_richness, group = 1), alpha = 0.95, size = 2, color = "#B22C2CFF") +
  #  geom_line(data = richness, aes(x = Condition, y = median_richness, group = 1), alpha = 0.95, linewidth = 0.8, color = "#B22C2CFF")
#print(overlay_plot)
print(p2)
dev.off()






