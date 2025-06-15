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

All_virus_data <- read.csv("2_All_virus_data_select_info.csv")
table(All_virus_data$Virus_name)
head(All_virus_data)

Bcell<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_BCR/1_BCR_info/data/1_Bcell_BCR.rds")
Bcell.df <- as.data.frame(Bcell@meta.data[ ,colnames(Bcell@meta.data)[c(19:25,68:107)]])
write.csv(Bcell.df,  "data/0_Bcell_BCR_df.csv")
colnames(Bcell.df)
table(Bcell.df$bcr_chain)

# 提取第一条TRB链的函数
extract_first_trb <- function(cdr3_string) {
  # 使用正则表达式查找所有TRB链
  trb_matches <- unlist(regmatches(cdr3_string, gregexpr("IGH:[^;]+", cdr3_string)))
  # 如果存在TRB链，返回第一条；否则返回NA
  if (length(trb_matches) > 0) {
    return(trb_matches[1])
  } else {
    return(NA)
  }
}

# 应用函数提取所有数据中的第一条TRB链
Bcell$bcr_IGH_cdr3<- sapply(Bcell$bcr_cdr3s_aa, extract_first_trb)

# 查看结果
head(Bcell$bcr_IGH_cdr3)
Bcell$cdr3_IGH_aa <-sub("IGH:", "", Bcell$bcr_IGH_cdr3)
head(Bcell$cdr3_IGH_aa)

length_distribution <- table(nchar(as.character(Bcell$cdr3_IGH_aa)))
barplot(length_distribution, 
        main = "Distribution of sequence length in our data", 
        xlab = "Sequence length", 
        ylab = "Frequency")
length_distribution <- table(nchar(as.character(All_virus_data$cdr3_aa_heavy)))
barplot(length_distribution, 
        main = "Distribution of sequence length in our data", 
        xlab = "Sequence length", 
        ylab = "Frequency")


tail(All_virus_data$cdr3_aa_heavy,n=20)
head(Bcell$cdr3_IGH_aa)

Bcell_data <- Bcell@meta.data %>%
                  filter(! is.na(cdr3_IGH_aa))%>% 
                  select(bcr_barcode,ct_level_4,Donor,Condition,
                         bcr_chain,bcr_v_gene,bcr_d_gene,bcr_j_gene,
                         bcr_cdr3s_aa,bcr_cdr3_nt,bcr_cdr3,cdr3_IGH_aa)

# 计算数据的行数
num_rows <- nrow(All_virus_data)

# 创建一个序列, 每20000个数据为一组
group_index <- ceiling(seq(1, num_rows) / 100000)

# 使用分组索引切割数据
split_datasets <- split(All_virus_data, group_index)

process_data <- function(data) {
  Bcell_data %>% stringdist_join(data, by=c("cdr3_IGH_aa"="cdr3_aa_heavy"), method="lv",
                  mode="inner", max_dist=3, distance_col="Distance_cdr3")
}

results_3 <- lapply(split_datasets, process_data)

# 合并结果
final_result_3 <- do.call(rbind, results_3)

nrow(final_result_3)
length(unique(final_result_3$cdr3_IGH_aa))
table(final_result_3$Virus_name)

head(final_result_3)

table(final_result_3$Distance_cdr3)

write.csv(final_result_3,"data/2_Bcell_virus_bcr_distance_3_IGH.csv")

res_distance_2 <- final_result_3%>% filter(Distance_cdr3 ==2)
nrow(res_distance_2)
length(unique(res_distance_2$cdr3_IGH_aa))
write.csv(res_distance_2,"data/2_Bcell_virus_bcr_distance_2_IGH.csv")

colnames(final_result_3)

process_data <- function(file_path) {
  # 读取数据
  res_distance <- read.csv(file_path)
  # 统计并标记 TCR 序列的一对多情况
  res_distance <- res_distance %>%
                    group_by(cdr3_IGH_aa) %>%
                    mutate(cdr3_dis_virus_kind = length(unique(Virus_name)))
  print("-----------------------------------------------------------------------------------")
  print("原序列的一对多情况统计：")
  print(table(res_distance$cdr3_dis_virus_kind))
 # 过滤掉一对多情况的序列
  res_distance <- res_distance %>% filter(cdr3_dis_virus_kind == 1)
  # 打印清理后的病毒分布
  print("清理后的病毒分布：（全为一一对应的序列）")
  print(table(res_distance$Virus_name))
    # 过滤出 SARS-CoV-2 的 TCR 序列
  res_distance_COVID <- res_distance %>% filter(Virus_name == "COVID")
  print(paste("新冠TCR序列数：", nrow(res_distance_COVID)))
  
  # 打印新冠TCR的不同序列数
  print(paste("其中新冠TCR不同序列统计数：", length(unique(res_distance_COVID$cdr3_IGH_aa))))
  return(res_distance)
}

# 处理三个不同的文件
res_distance_3 <- process_data("data/2_Bcell_virus_bcr_distance_3_IGH.csv")
res_distance_2 <- process_data("data/2_Bcell_virus_bcr_distance_2_IGH.csv")

getwd()

Bcell_data$Virus_name <- res_distance_3$Virus_name[match(Bcell_data$cdr3_IGH_aa,res_distance_3$cdr3_IGH_aa)]
Bcell_data$Chain_distance <- res_distance_3$Distance_cdr3[match(Bcell_data$cdr3_IGH_aa,res_distance_3$cdr3_IGH_aa)]
table(Bcell_data$Virus_name)
table(Bcell_data$Chain_distance)
write.csv(Bcell_data,"data/Bcell_data_res_distance_3.csv")

species_table <- as.data.frame(table(Bcell_data$Virus_name))
colnames(species_table) <- c("Virus_name", "Count")

# 计算每种物种的百分比
species_table <- species_table %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  arrange(desc(Percentage))  # 按百分比从高到低排序

species_table$Virus_name <- factor(species_table$Virus_name,
                               levels = unique(species_table$Virus_name))

colors <- colorRampPalette(brewer.pal(12, "Paired"))(7)# 设置颜色（使用 RColorBrewer 中的调色板）

# 生成饼图
pdf("2_Bcell_Virus_name_Pie_distance_3.pdf", 5, 5)
ggplot(species_table, aes(x = "", y = Percentage, fill = Virus_name)) +
  geom_bar(width = 0.9, stat = "identity", alpha = 0.95, color = "white") +  # 设置透明度
  coord_polar(theta = "y") +  # 使用极坐标绘制饼图
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_void() +  # 移除背景元素
  labs(fill = "Virus name", title = "Virus BCR in Database") +  # 添加标签
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

# 生成横向条形图，按 Count 对 Virus_name 进行排序，并显示数字
pdf("2_Bcell_Virus_name_Bar_Plot.pdf", 5, 4)  # 调整图形大小为6x4

ggplot(species_table, aes(y = reorder(Virus_name, Count), x = Count, fill = Virus_name)) +  # 对 Virus_name 进行排序
  geom_bar(stat = "identity", alpha = 0.95, color = "white") +  # 绘制条形图
  geom_text(aes(label = Count), hjust = -0.2, color = "black", size = 3.5) +  # 在条形图右侧显示数字
  scale_fill_manual(values = colors) +  # 设置调色板
  theme_minimal() +  # 使用简洁主题
  labs(fill = "Virus name", title = "Virus BCR in Database",
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



library("tidyverse")
library("data.table")
library(ggplot2)
library(dplyr)
library(Biostrings)

dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/6_BCR/3_BCR_SHM'

# SHM 计算函数
# 读取IMGT参考序列（修改路径为你保存的FASTA文件路径）
IGHV_reference_sequences <- readDNAStringSet(paste0(dir_path,"/IMGT/IGHV.fasta.txt"))
IGKV_reference_sequences <- readDNAStringSet(paste0(dir_path,"/IMGT/IGKV.fasta.txt"))
IGLV_reference_sequences <- readDNAStringSet(paste0(dir_path,"/IMGT/IGLV.fasta.txt"))

# 存放到一个列表中便于查找
reference_sequences <- list(
  IGHV = IGHV_reference_sequences,
  IGKV = IGKV_reference_sequences,
  IGLV = IGLV_reference_sequences
)


get_reference_sequence <- function(v_gene) {
  # 假设v_gene是一个字符，例如 "IGLV2-14"
  chain_type <- substr(v_gene, 1, 4)  # 提取链类型
  
  # 检查这个链类型是否存在于参考序列集合中
  if (chain_type %in% names(reference_sequences)) {
    seq_set <- reference_sequences[[chain_type]]
    
    # 在序列集中匹配基因名
    matched_seq <- seq_set[grep(v_gene, names(seq_set))]
    
    # 如果找到匹配的序列，返回；否则返回NULL
    if (length(matched_seq) > 0) {
      return(matched_seq[[1]])
    } else {
      warning(paste("No reference sequence found for", v_gene))
      return(NULL)
    }
  } else {
    warning(paste("Invalid chain type for", v_gene))
    return(NULL)
  }
}

library(stringdist)

calculate_shm_rate <- function(sequence, reference_sequence) {
  # 如果参考序列或序列为空，返回NA
  if (is.null(reference_sequence) || is.null(sequence)) {
    return(NA)
  }

  # 输出调试信息
  # print(paste("Sequence:", sequence, "Length:", nchar(sequence)))
  # print(paste("Reference Sequence:", reference_sequence, "Length:", nchar(reference_sequence)))

  # 计算序列之间的Levenshtein距离
  distance <- stringdist::stringdist(sequence, reference_sequence, method = "lv")

  # 确定参考序列的长度
  ref_length <- nchar(reference_sequence)

  # 计算突变率
  SHM_rate <- distance / ref_length
  return(SHM_rate)
}

Bcell_data <- read.csv("data/Bcell_data_res_distance_3.csv")
colnames(Bcell_data)
rownames(Bcell_data) <- Bcell_data$X
Bcell_data$X <- NULL

Bcell.df <- read.csv("data/0_Bcell_BCR_df.csv")
colnames(Bcell.df)

table(Bcell.df$bcr_c_gene)

table(Bcell_data$Virus_name)

Bcell_data_COVID <- Bcell_data %>% filter(Virus_name == 'COVID') 
plot_df <- Bcell.df  %>%
  filter(!is.na(bcr_c_gene)) %>%
  filter(bcr_c_gene %in% c("IGHM",  "IGHD",  "IGHA1",  "IGHA2", 
                                "IGHG1",  "IGHG2","IGHG3",  "IGHG4" )) %>%
  filter(bcr_barcode %in% Bcell_data_COVID$bcr_barcode) %>%
  count(ct_level_4, bcr_c_gene) %>%
  mutate(ct_level_4 = forcats::fct_rev(as.factor(ct_level_4)))
head(plot_df)
unique(plot_df$bcr_c_gene)
table(plot_df$bcr_c_gene)

plot_df

plot_df$ct_level_4 <- factor(plot_df$ct_level_4,
                               levels=rev(c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                        'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-',
                                        'Plasma')))
unique(plot_df$ct_level_4)

# 设置BCR的链类型顺序
plot_df$bcr_c_gene <- factor(
  plot_df$bcr_c_gene,
    levels = c(
    "IGHM",  # IgM: 第一个产生的抗体
    "IGHD",  # IgD: 参与 B 细胞的成熟
    "IGHA1", # IgA: 主要在粘膜和体液中
    "IGHA2", # IgA: 主要在粘膜和体液中
    "IGHG1", # IgG1: 最常见的 IgG 亚型
    "IGHG2", # IgG2: 与特定病原体的免疫反应相关
    "IGHG3", # IgG3: 强烈的免疫反应，通常与感染相关
    "IGHG4" # IgG4: 与过敏反应和特定免疫反应相关
    )
)

# 图1：BCR链类型在细胞亚群中的分布
pdf("3_COVID_DB_CVV_BCR_IGH_c_gene.pdf",4.2,2.5)
ggplot(plot_df, aes(ct_level_4, n, fill = bcr_c_gene)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(
    values =  c(
      "IGHM"  = "#51A3CCFF",  # IgM
      "IGHD"  = "#f18800",  # IgD
      "IGHA1" = "#85B22CFF",  # IgA1
      "IGHA2" = "#B22C2CFF",  # IgA2
      "IGHG1" = "#ffd401",  # IgG1
      "IGHG2" = "#e20612",  # IgG2
      "IGHG3" = "#2e409a",  # IgG3
      "IGHG4" = "#E57E7EFF"  # IgG4
    )
  ) +
  labs(x = "", y = "Proportion", fill = 'BCR isotype') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.3)#,
    #legend.position = "none"
  ) +
  coord_flip()
dev.off()











