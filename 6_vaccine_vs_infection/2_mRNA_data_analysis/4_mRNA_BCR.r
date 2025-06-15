library("tidyverse")
library("data.table")
library(ggplot2)
library("Matrix")
library(Seurat)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tidyr)
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

mRNA_obj <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/data/GSE247917_cv19_combined_obj.rds")

Idents(mRNA_obj) <- 'group'
mRNA_obj <- subset(mRNA_obj,idents=c('HC','Vax','booster'))

# 筛选数据
Idents(mRNA_obj) <- 'celltype.l1'
Bcell <- subset(mRNA_obj,  idents = c('B/Plasma'))

colnames(Bcell@meta.data)

nrow(Bcell@meta.data)

Bcell@meta.data$celltype.l2 <- factor(Bcell@meta.data$celltype.l2,
                                     level=unique(Bcell@meta.data$celltype.l2))
unique(Bcell@meta.data$celltype.l2)

table(Bcell@meta.data$celltype.l2)

write.csv(Bcell@meta.data,"data/0_mRNA_bcr_info.csv")

saveRDS(Bcell,  "data/1_mRNA_Bcell_bcr.rds")

Bcell <- readRDS("data/1_mRNA_Bcell_bcr.rds")

nclos_clonotype <-  c("lightgrey", "darkgrey", "#E0C5BEFF",
                      "#F39B7FB2", "#E57E7EFF", "#DE5C00FF", "#B22C2CFF", "#573333FF")

# 定义函数来处理细胞对象的克隆型数据
process_clonotype_data <- function(cell_type_obj,ct_seq) {
  # 设置 celltype.l2 因子顺序
  cell_type_obj$celltype.l2 <- factor(
    cell_type_obj$celltype.l2,
    levels = ct_seq
  )
  
  # 设置 donor_id 因子顺序
  cell_type_obj$donor_id <- factor(
    cell_type_obj$donor_id,
    levels = unique(cell_type_obj$donor_id)
  )
  
  # 设置 experiment 因子顺序
  cell_type_obj$experiment <- factor(
    cell_type_obj$experiment,
    levels = unique(cell_type_obj$experiment)
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
    levels = rev(c("Not detected", "n = 1", "1 < n <= 3", "3 < n <= 5",  "n > 5"))
  )
  
  return(cell_type_obj)
}

# 定义通用绘图函数
plot_clonotype_umap <- function(cell_type_obj, file_name, nclos_clonotype,w,h,pt) {
  pdf(file_name, w,h)
  print(
    DimPlot(
      object = cell_type_obj, reduction = "umap", raster=FALSE,label = FALSE, pt.size = pt, order = TRUE,
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
        raster=FALSE,
      pt.size = 0.01, 
      order = TRUE,
      group.by = "Clonotype_num_ld", 
      split.by = split_by,
      cols = nclos_clonotype
    )
  )
}
    

Bcell@meta.data <- Bcell@meta.data %>%
  mutate(
    bcr_cdr3 = ifelse(
      is.na(BCR_IGH.CDR3.1) & is.na(BCR_IGL.CDR3.1) & is.na(BCR_IGK.CDR3.1),
      NA,
      paste0(
        ifelse(!is.na(BCR_IGH.CDR3.1), paste0("IGH:", BCR_IGH.CDR3.1), ""),
        ifelse(!is.na(BCR_IGL.CDR3.1), paste0("; IGL:", BCR_IGL.CDR3.1), ""),
        ifelse(!is.na(BCR_IGK.CDR3.1), paste0("; IGK:", BCR_IGK.CDR3.1), "")
      ) %>%
        # 清理开头或多余的分号
        gsub("^; ", "", .) %>%
        gsub("; $", "", .)
    )
  )

# 查看结果
head(Bcell@meta.data$bcr_cdr3)

unique(Bcell$celltype.l2)

unique(Bcell$celltype.l2)
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0", "#DE5C00FF", "#B22C2CFF","#573333FF")
ct_seq <-c('Resting B','Activated B','Memory B','B Exhausted','Plasmablasts')
clonotype_num <- as.data.frame(table(Bcell$bcr_cdr3))
Bcell$Clonotype_num <- clonotype_num$Freq[match(Bcell$bcr_cdr3, clonotype_num$Var1)]
Bcell <- process_clonotype_data(Bcell,ct_seq)

plot_clonotype_umap(Bcell, "figures/1_Bcell_clonotype_umap.pdf", nclos_clonotype,5,3.5,0.0001)

pdf("figures/1_Bcell_clonotype_umap_split.pdf", 18, 4)
split_conditions <- c("donor_id", "group", "timepoint", "celltype.l2")
lapply(split_conditions, function(split_by) plot_dim(Bcell, split_by))
dev.off()

table(Bcell$Clonotype_num)

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(Bcell$celltype.l2, Bcell$Clonotype_num_ld))
colnames(data) <- c("celltype.l2", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/1_Bcell_Bar_Clonotype_Levels_all.pdf",3.6,2.5)
ggplot(data, aes(x = celltype.l2, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(
       x = "Cell type", 
       y = "Proportion", 
       fill = "Clonotype size") +
  theme_minimal()+
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 8, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5,color='black'), # 保留 X 轴和 Y 轴的线
      panel.grid = element_blank(),
     axis.text.y = element_text(color='black'),
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1,color='black')# 设置 X 轴文本
  ) 
ggplot(data, aes(x = celltype.l2, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_mapping) +
  labs(
       x = "Cell type", 
       y = "Count", 
       fill = "Clonotype size") +
  theme_minimal()+
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 8, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5,color='black'), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(),
    axis.text.y = element_text(color='black'),
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1,color='black')# 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)
Bcell_f <- Bcell@meta.data%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
Bcell_f$Clonotype_num_ld <- factor(Bcell_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(Bcell_f$celltype.l2, Bcell_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("celltype.l2", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/1_Bcell_Bar_Clonotype_Levels.pdf",3.6,2.5)
ggplot(data, aes(x = celltype.l2, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number", 
       x = "Cell type", 
       y = "Proportion", 
       fill = "Clonotype size") +
  theme_minimal()+
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 8, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5,color='black'), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(),
      axis.text.y = element_text(color='black'),
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1,color='black'), # 设置 X 轴文本
  ) 
dev.off()

# 加载 dplyr 包
library(dplyr)

# 定义映射关系，包含 ILC 类别
time_mapping <- c(
  'HC' = 'A0',
  'Vax_d7' = 'A1',
  'Vax_d10' = 'A2',
    'Vax_d14' = 'A2',
    'Vax_d20' = 'B0',
    'Vax_d21' = 'B0',
     'Vax_d28' = 'B1',
    'Vax_d29' = 'B1',
    'Vax_d35' = 'B2',
    'Vax_d36' = 'B2',
     'booster_d0' = 'C0',
    'booster_d7' = 'C1',
     'booster_d28' = 'C2',
    'booster_d120' = 'C3',           
     'acute' = 'acute',
    'cv19_d11' = 'acute',
    'convalescent' = 'convalescent'
)

# 使用 dplyr::mutate 和 recode 应用映射
Bcell@meta.data <- Bcell@meta.data %>%
  mutate(Condition = recode(timepoint, !!!time_mapping))
summary(Bcell@meta.data$timepoint)
summary(Bcell@meta.data$Condition)
# 查看映射结果
unique(Bcell@meta.data[, c("timepoint", "Condition")])

Bcell$Condition <- factor(Bcell$Condition,
                         levels=c('A0','A1','A2',
                                  'B0','B1','B2',
                                  'C0','C1','C2','C3'))
#'acute','convalescent'
#table(Bcell$Clonotype_num,Bcell$Condition)

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(Bcell$Condition, Bcell$Clonotype_num_ld))
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/1_Bcell_Bar_Clonotype_Levels_Condition_all.pdf",3.6,2.5)
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(
       x = "Condition", 
       y = "Proportion", 
       fill = "Clonotype size") +
  theme_minimal()+
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 8, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5,color='black'), # 保留 X 轴和 Y 轴的线
      panel.grid = element_blank(),
     axis.text.y = element_text(color='black'),
    axis.text.x = element_text(size = 8, color='black')# 设置 X 轴文本 #angle = 45,hjust = 1,
  ) 
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_mapping) +
  labs(
       x = "Condition", 
       y = "Count", 
       fill = "Clonotype size") +
  theme_minimal()+
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 8, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5,color='black'), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(),
    axis.text.y = element_text(color='black'),
    axis.text.x = element_text(size = 8, color='black')# 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)
Bcell_f <- Bcell@meta.data%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
Bcell_f$Clonotype_num_ld <- factor(Bcell_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(Bcell_f$Condition, Bcell_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/1_Bcell_Bar_Clonotype_Levels_Condition.pdf",3.6,2.5)
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number", 
       x = "Condition", 
       y = "Proportion", 
       fill = "Clonotype size") +
  theme_minimal()+
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 8, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5,color='black'), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(),
      axis.text.y = element_text(color='black'),
    axis.text.x = element_text(size = 8,color='black'), # 设置 X 轴文本
  ) 
dev.off()





All_virus_data <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/6_BCR/2_BCR_viral/2_All_virus_data_select_info.csv")
table(All_virus_data$Virus_name)
head(All_virus_data)

tail(All_virus_data$cdr3_aa_heavy,n=2)
head(Bcell$BCR_IGH.CDR3.1)

colnames(Bcell@meta.data)

Bcell@meta.data$bcr_barcode <- rownames(Bcell@meta.data)
head(Bcell@meta.data$bcr_barcode)

Bcell_data <- Bcell@meta.data %>%
                  filter(!is.na(BCR_IGH.CDR3.1))%>% 
                  select( bcr_barcode,celltype.l1,celltype.l2,Condition,
                      timepoint,day,assay,
                         BCR_IGH.CDR3.1,BCR_IGH.V.1,BCR_IGH.D.1,
                         BCR_IGH.J.1,BCR_IGH.C.1)
head(Bcell_data,n=1)

# 计算数据的行数
num_rows <- nrow(All_virus_data)

# 创建一个序列, 每20000个数据为一组
group_index <- ceiling(seq(1, num_rows) / 100000)

# 使用分组索引切割数据
split_datasets <- split(All_virus_data, group_index)

process_data <- function(data) {
  Bcell_data %>% stringdist_join(data, by=c("BCR_IGH.CDR3.1"="cdr3_aa_heavy"), method="lv",
                  mode="inner", max_dist=3, distance_col="Distance_cdr3")
}

results_3 <- lapply(split_datasets, process_data)

# 合并结果
final_result_3 <- do.call(rbind, results_3)

nrow(final_result_3)
length(unique(final_result_3$BCR_IGH.CDR3.1))
table(final_result_3$Virus_name)

write.csv(final_result_3,"data/2_mRNA_Bcell_virus_bcr_distance_3_IGH.csv")

res_distance_2 <- final_result_3%>% filter(Distance_cdr3 ==2)
nrow(res_distance_2)
length(unique(res_distance_2$BCR_IGH.CDR3.1))
write.csv(res_distance_2,"data/2_Bcell_virus_bcr_distance_2_IGH.csv")

## 删除一对多的数据
process_data <- function(file_path) {
  # 读取数据
  res_distance <- read.csv(file_path)
  # 统计并标记 TCR 序列的一对多情况
  res_distance <- res_distance %>%
                    group_by(BCR_IGH.CDR3.1) %>%
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
  print(paste("其中新冠TCR不同序列统计数：", length(unique(res_distance_COVID$BCR_IGH.CDR3.1))))
  return(res_distance)
}

# 处理三个不同的文件
res_distance_3 <- process_data("data/2_mRNA_Bcell_virus_bcr_distance_3_IGH.csv")
res_distance_2 <- process_data("data/2_Bcell_virus_bcr_distance_2_IGH.csv")

getwd()

Bcell_data$Virus_name <- res_distance_3$Virus_name[match(Bcell_data$BCR_IGH.CDR3.1,res_distance_3$BCR_IGH.CDR3.1)]
Bcell_data$Chain_distance <- res_distance_3$Distance_cdr3[match(Bcell_data$BCR_IGH.CDR3.1,res_distance_3$BCR_IGH.CDR3.1)]
table(Bcell_data$Virus_name)
table(Bcell_data$Chain_distance)
write.csv(Bcell_data,"data/2_Bcell_data_res_distance_3.csv")

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
pdf("figures/2_Bcell_Virus_name_Pie_distance_3.pdf", 5, 5)
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
pdf("figures/2_Bcell_Virus_name_Bar_Plot.pdf", 5, 4)  # 调整图形大小为6x4

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

Bcell <- readRDS("data/1_mRNA_Bcell_bcr.rds")

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

table(Bcell$BCR_IGH.C.1)
plot_df <- data.frame(Bcell@meta.data) %>%
  filter(!is.na(BCR_IGH.C.1)) %>%
  filter(BCR_IGH.C.1 %in% c("IGHM",  "IGHD",  "IGHA1",  "IGHA2", 
                                "IGHG1",  "IGHG2","IGHG3",  "IGHG4" )) %>%
  count(celltype.l2, BCR_IGH.C.1) %>%
  mutate(celltype.l2 = forcats::fct_rev(as.factor(celltype.l2)))
head(plot_df)
unique(plot_df$BCR_IGH.C.1)
table(plot_df$BCR_IGH.C.1)

# 设置BCR的链类型顺序
plot_df$BCR_IGH.C.1 <- factor(
  plot_df$BCR_IGH.C.1,
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
unique(plot_df$celltype.l2)
plot_df$celltype.l2 <- factor(plot_df$celltype.l2,
                              c('Plasmablasts',
                                    'B Exhausted',
                                    'Memory B',
                                'Activated B',
                                    'Resting B'))

write.csv(plot_df,"data/3_BCR_IGH_plot_df.csv")

# 图1：BCR链类型在细胞亚群中的分布
pdf("figures/3_mRNA_BCR_IGH_c_gene.pdf",2.7,1.5)
ggplot(plot_df, aes(celltype.l2, n, fill = BCR_IGH.C.1)) +
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
    text = element_text(size = 6),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
       legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.1, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
    #legend.position = "none"
  ) +
  coord_flip()
dev.off()

Bcell_data <- read.csv("data/2_Bcell_data_res_distance_3.csv")
colnames(Bcell_data)
rownames(Bcell_data) <- Bcell_data$X
Bcell_data$X <- NULL
table(Bcell_data$BCR_IGH.C.1)

table(Bcell_data$Virus_name)

plot_df <- Bcell_data %>%
  filter(!is.na(BCR_IGH.C.1)) %>%
  filter(BCR_IGH.C.1 %in% c("IGHM",  "IGHD",  "IGHA1",  "IGHA2", 
                                "IGHG1",  "IGHG2","IGHG3",  "IGHG4" )) %>%
  filter(Virus_name == 'COVID') %>%
  count(celltype.l2, BCR_IGH.C.1) %>%
  mutate(celltype.l2 = forcats::fct_rev(as.factor(celltype.l2)))
head(plot_df)
unique(plot_df$BCR_IGH.C.1)
table(plot_df$BCR_IGH.C.1)

plot_df

# 设置BCR的链类型顺序
plot_df$BCR_IGH.C.1 <- factor(
  plot_df$BCR_IGH.C.1,
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
unique(plot_df$celltype.l2)
plot_df$celltype.l2 <- factor(plot_df$celltype.l2,
                              c('Plasmablasts',
                                    'B Exhausted',
                                    'Memory B',
                                'Activated B',
                                    'Resting B'))

# 图1：BCR链类型在细胞亚群中的分布
pdf("figures/3_COVID_DB_mRNA_BCR_IGH_c_gene.pdf",4,3)
ggplot(plot_df, aes(celltype.l2, n, fill = BCR_IGH.C.1)) +
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







Idents(Bcell) <- 'group'
Bcell <- subset(Bcell,idents=c('HC','Vax','booster'))
Bcell$new_ct <- 'Plasmablasts'
Bcell$new_ct[Bcell$celltype.l2 %in% c('Resting B', 'Memory B',
                                 'B Exhausted','Activated B')] <- 'B cell'
nrow(Bcell@meta.data)

table(Bcell$new_ct)

source("/share/home/qlab/projects/project_cvv/yyp_results_3/yyp_function.R") #加载到当前的R环境中

dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/2_Bcell/4_diff_gene_in_Bcell'

GO_result <- read.csv(paste0(dir_path,"/data/3_Bm_CD27+_IGHM+_vs_Others_GO_results.csv"))

#展示GO某个的基因Moudule_scor时期变化
Bcell <- analyze_and_visualize_GO_genes (GO_result, Bcell, 
                                        Description_name="regulation of B cell activation",
                                        ctrl_genes = NULL, 
                                        Celltype=  'new_ct',
                                         Condition = 'Condition',
                                   module_score_name = 'IGHM_act', module_score_real_name = 'IGHM_act1',
                                   output_heatmap_file = "figures/4_mRNA_Bm_CD27+_IGHM+_Module_gene.pdf",
                                         col_fun_values = c(0, 0.5, 1, 1.8, 2.4, 3),
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 5, heatmap_height = 3)

colnames(GO_result)

GO_result <- read.csv(paste0(dir_path,"/data/1_Plasma_vs_Others_GO_results.csv"))


#a <- "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains"
#GO_result$Description[GO_result$Description == a] <- "adaptive immune response"

#展示GO某个的基因Moudule_scor时期变化
Bcell <- analyze_and_visualize_GO_genes (GO_result, Bcell, 
                                        Description_name="B cell receptor signaling pathway",
                                        ctrl_genes = NULL, 
                                        Celltype=  'new_ct',
                                         Condition = 'Condition',
                                   module_score_name = 'IGHM_act', module_score_real_name = 'IGHM_act1',
                                   output_heatmap_file = "figures/4_mRNA_Plasma_Module_gene.pdf",
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 5, heatmap_height = 3)

GO_result <- read.csv(paste0(dir_path,"/data/2_Bm_CD27+_IGHM+_SOX5_vs_Others_GO_results.csv"))


#展示GO某个的基因Moudule_scor时期变化
Bcell <- analyze_and_visualize_GO_genes (GO_result, Bcell, 
                                        Description_name="immune response-regulating signaling pathway",
                                        ctrl_genes = NULL, 
                                        Celltype=  'new_ct',
                                         Condition = 'Condition',
                                   module_score_name = 'IGHM_act', module_score_real_name = 'IGHM_act1',
                                   output_heatmap_file = "figures/4_mRNA_Bm_CD27+_IGHM+_SOX5_Module_gene.pdf",
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 5, heatmap_height = 3)

GO_result <- read.csv(paste0(dir_path,"/data/5_Bm_CD27+_IGHM-_vs_Others_GO_results.csv"))


#展示GO某个的基因Moudule_scor时期变化
Bcell <- analyze_and_visualize_GO_genes (GO_result, Bcell, 
                                        Description_name="immune response-activating signaling pathway",
                                        ctrl_genes = NULL, 
                                        Celltype=  'new_ct',
                                         Condition = 'Condition',
                                   module_score_name = 'IGHM_act', module_score_real_name = 'IGHM_act1',
                                   output_heatmap_file = "figures/4_mRNA_Bm_CD27+_IGHM-_Module_gene.pdf",
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 5, heatmap_height = 3)

GO_result <- read.csv(paste0(dir_path,"/data/4_Bmn_IGHD_vs_Others_GO_results.csv"))


#展示GO某个的基因Moudule_scor时期变化
Bcell <- analyze_and_visualize_GO_genes (GO_result, Bcell, 
                                        Description_name="MHC class II protein complex assembly",
                                        ctrl_genes = NULL, 
                                        Celltype=  'new_ct',
                                         Condition = 'Condition',
                                   module_score_name = 'IGHM_act', module_score_real_name = 'IGHM_act1',
                                   output_heatmap_file = "figures/4_mRNA_Bmn_IGHD_Module_gene.pdf",
                                    #col_fun_values = c(3, 3.5, 4, 4.5, 4.8, 4.9),  
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 5, heatmap_height = 3)










