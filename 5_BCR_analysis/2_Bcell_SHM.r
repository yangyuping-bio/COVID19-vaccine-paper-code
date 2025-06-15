library("tidyverse")
library("data.table")
library(ggplot2)
library(dplyr)
library(Biostrings)

getwd()

Bcell<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_BCR/1_BCR_info/data/1_Bcell_BCR.rds")

# SHM 计算函数
# 读取IMGT参考序列（修改路径为你保存的FASTA文件路径）
IGHV_reference_sequences <- readDNAStringSet("IMGT/IGHV.fasta.txt")
IGKV_reference_sequences <- readDNAStringSet("IMGT/IGKV.fasta.txt")
IGLV_reference_sequences <- readDNAStringSet("IMGT/IGLV.fasta.txt")

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

#BCR数据集
Bcell_SHM <- as.data.frame(Bcell@meta.data) %>%
  rowwise() %>%
  mutate(full_bcr_seq = paste(
      bcr_fwr1_nt, bcr_cdr1_nt, bcr_fwr2_nt, bcr_cdr2_nt,
      bcr_fwr3_nt, bcr_cdr3_nt, bcr_fwr4_nt, sep = "")
  ) %>% 
  mutate( shm_rate_ratio = calculate_shm_rate(full_bcr_seq, get_reference_sequence(bcr_v_gene))
      )
range(Bcell_SHM$shm_rate_ratio , na.rm = TRUE)
min_shm <- min(Bcell_SHM$shm_rate_ratio, na.rm = TRUE)
max_shm <- max(Bcell_SHM$shm_rate_ratio, na.rm = TRUE)
Bcell_SHM$SHM_rate <- (Bcell_SHM$shm_rate_ratio - min_shm) / (max_shm - min_shm) * 10
range(Bcell_SHM$SHM_rate , na.rm = TRUE)

write.csv(Bcell_SHM,"Bcell_BCR_SHM.csv")

# 加载必要的库
library(dplyr)
library(ggplot2)
library(forcats)
library(cowplot)

table(Bcell_SHM$bcr_v_gene)
table(Bcell_SHM$bcr_d_gene)
table(Bcell_SHM$bcr_j_gene)
table(Bcell_SHM$bcr_c_gene)

unique(Bcell_SHM$ct_level_4)

#数据准备与SHM率计算
# 过滤出我们需要的数据，
 Bcell_SHM$ct_level_4 <- factor(Bcell_SHM$ct_level_4,
                               levels=c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                        'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-',
                                        'Plasma'))
plot_df <- Bcell_SHM %>%
  filter(!is.na(bcr_c_gene)) %>%
  count(ct_level_4, bcr_c_gene) %>%
  mutate(ct_level_4 = forcats::fct_rev(as.factor(ct_level_4)))

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
    "IGHG4", # IgG4: 与过敏反应和特定免疫反应相关
    "IGKC",  # kappa 轻链
    "IGLC1", # lambda 轻链
    "IGLC2", # lambda 轻链
    "IGLC3", # lambda 轻链
    "IGLC7"  # lambda 轻链
    )
)

# 图1：BCR链类型在细胞亚群中的分布
p1 <- ggplot(plot_df, aes(ct_level_4, n, fill = bcr_c_gene)) +
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
      "IGHG4" = "#E57E7EFF",  # IgG4
      "IGKC"  = "#DDC0DC",  # kappa 轻链
      "IGLC1" = "#CDA2CC",  # lambda 轻链 1
      "IGLC2" = "#BE84BC",  # lambda 轻链 2
      "IGLC3" = "#AF67AC",  # lambda 轻链 3
      "IGLC7" = "#9F499C"   # lambda 轻链 7
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
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  ) +
  coord_flip()
p1

range(Bcell_SHM$SHM_rate , na.rm = TRUE)

#SHM水平分类与可视化
plot_df <- Bcell_SHM %>%
  filter(!is.na(bcr_c_gene)) %>%
  mutate(ct_level_4 = forcats::fct_rev(as.factor(ct_level_4)))

plot_df$SHM_levels <- ifelse(plot_df$SHM_rate < 3, "Low", "Median")
plot_df$SHM_levels <- ifelse(plot_df$SHM_rate > 5, "High", plot_df$SHM_levels)
plot_df$SHM_levels <- factor(plot_df$SHM_levels, levels = c("Low", "Median", "High"))

plot_df <- plot_df %>% 
  group_by(ct_level_4, SHM_levels) %>% 
  summarise(n = n()) %>%
  group_by(ct_level_4) %>% 
  mutate(sum_n = sum(n)) %>% 
  ungroup() %>% 
  mutate(percent = n / sum_n)

p2 <- ggplot(plot_df, aes(x = ct_level_4, y = percent, fill = SHM_levels)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "High" = "#78290f",
    "Median" = "#ff7d00",
    "Low" = "#ffecd1"
  )) +
  labs(x = "", y = "Proportion", fill = 'SHM level') +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  )
p2

#中位数SHM率可视化
plot_df <- Bcell_SHM %>%
  filter(!is.na(bcr_c_gene)) %>%
  mutate(ct_level_4 = forcats::fct_rev(as.factor(ct_level_4)))

plot_df <- plot_df %>% 
  group_by(ct_level_4) %>% 
  summarise(median_SHM = median(SHM_rate))

p3 <- ggplot(plot_df, aes(x = ct_level_4, y = median_SHM)) +
  geom_point(shape = 16, stroke = 0, size = 1) +
  labs(x = "", y = "SHM rate") +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.ticks.y  = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  )
p3

#将三个图表合并成一个布局
cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "h", rel_widths = c(1, 1, 0.4))

# 保存图表
ggsave("Bcell_SHM_with_light_chain.pdf", width = 6, height = 1.6)

Bcell_SHM <- Bcell_SHM %>% 
       filter(bcr_c_gene %in%  c("IGHM",  "IGHD",  "IGHA1",  "IGHA2", 
                                "IGHG1",  "IGHG2","IGHG3",  "IGHG4" 
                      ))

#数据准备与SHM率计算
# 过滤出我们需要的数据，
 Bcell_SHM$ct_level_4 <- factor(Bcell_SHM$ct_level_4,
                               levels=c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                        'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-',
                                        'Plasma'))
plot_df <- Bcell_SHM %>%
  filter(!is.na(bcr_c_gene)) %>%
  count(ct_level_4, bcr_c_gene) %>%
  mutate(ct_level_4 = forcats::fct_rev(as.factor(ct_level_4)))

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
p1 <- ggplot(plot_df, aes(ct_level_4, n, fill = bcr_c_gene)) +
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
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  ) +
  coord_flip()
p1

# 计算排除 NA 后的取值范围
range_value <- range(Bcell_SHM$SHM_rate, na.rm = TRUE)
print(range_value)

# 计算排除 NA 后的取值范围
summary(Bcell_SHM$SHM_rate, na.rm = TRUE)

#SHM水平分类与可视化
plot_df <- Bcell_SHM %>%
  filter(!is.na(bcr_c_gene)) %>%
  mutate(ct_level_4 = forcats::fct_rev(as.factor(ct_level_4)))

plot_df$SHM_levels <- ifelse(plot_df$SHM_rate < 4.688, "Low", "Median")
plot_df$SHM_levels <- ifelse(plot_df$SHM_rate > 5.523, "High", plot_df$SHM_levels)
plot_df$SHM_levels <- factor(plot_df$SHM_levels, levels = c("Low", "Median", "High"))

plot_df <- plot_df %>% 
  group_by(ct_level_4, SHM_levels) %>% 
  summarise(n = n()) %>%
  group_by(ct_level_4) %>% 
  mutate(sum_n = sum(n)) %>% 
  ungroup() %>% 
  mutate(percent = n / sum_n)

p2 <- ggplot(plot_df, aes(x = ct_level_4, y = percent, fill = SHM_levels)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "High" = "#810F7C",
    "Median" = "#FFA8C3",
    "Low" = "#ffecd1"
  )) +
  labs(x = "", y = "Proportion", fill = 'SHM level') +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  )
p2

#中位数SHM率可视化
plot_df <- Bcell_SHM %>%
  filter(!is.na(bcr_c_gene)) %>%
  mutate(ct_level_4 = forcats::fct_rev(as.factor(ct_level_4)))

plot_df <- plot_df %>% 
  group_by(ct_level_4) %>% 
  summarise(median_SHM = median(SHM_rate))

p3 <- ggplot(plot_df, aes(x = ct_level_4, y = median_SHM)) +
  geom_point(shape = 16, stroke = 0, size = 1) +
  labs(x = "", y = "SHM rate") +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.ticks.y  = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  )
p3

#将三个图表合并成一个布局
cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "h", rel_widths = c(1, 1, 0.4))

# 保存图表
ggsave("Bcell_SHM_with_legend.pdf", width = 8, height = 2.8)

#将三个图表合并成一个布局
cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "h", rel_widths = c(1, 1, 0.4))

# 保存图表
ggsave("Bcell_SHM.pdf", width = 6, height = 1.6)




