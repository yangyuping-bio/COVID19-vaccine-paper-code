library("tidyverse")
library(Seurat)
library(future)
# 加载所需的包
library(readxl)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(dplyr)

#ct_merge
colorl2 <- c(
  Bcell="#AB3282",  # B细胞类型
  CD4T="#6778AE", # CD4 T细胞类型
  CD8T="#53A85F", # CD8 T细胞类型
  cDC="#9FA3A8",  # CDC细胞类型
  gdT="#23452F",  #gdT
  HPSC="#585658",#HPSC
  ILC="#F1BC62", #ILC
  MAIT="#3A6963",#MAIT
  Monocyte="#E95C39",# Mono细胞类型 
  NK="#E1A111", # NK细胞类型
  pDC="#938175", 
  #Platelet="#BD956A",
  Tdn="#495608"# 其他细胞类型
)

dir_path <- '/share/home/qlab/projects/project_mRNA/yyp_results_3/0_paper_data/6_mRNA_data/4_new_celltype'

sample.integrated<- readRDS(paste0(dir_path,"/1_mRNA_final_CCA_integrated.rds"))

colnames(sample.integrated@meta.data)

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
sample.integrated@meta.data <- sample.integrated@meta.data %>%
  mutate(Condition = recode(timepoint, !!!time_mapping))
summary(sample.integrated@meta.data$timepoint)
summary(sample.integrated@meta.data$Condition)
# 查看映射结果
unique(sample.integrated@meta.data[, c("timepoint", "Condition")])

table(sample.integrated@meta.data$Condition)

unique(sample.integrated@meta.data$group)

Idents(sample.integrated) <- 'group'
mRNA_part <- subset(sample.integrated,
                    idents=c('Vax','HC','booster'))
table(mRNA_part$Condition)

sample.integrated <- mRNA_part



unique(sample.integrated$Condition)
unique(sample.integrated$ct_merge)

DefaultAssay(sample.integrated) <- "RNA"
sample.integrated$ct_merge_condition <- paste(sample.integrated$ct_merge, sample.integrated$Condition, sep = "_")
Idents(sample.integrated) <- "ct_merge_condition"
length(table(sample.integrated$ct_merge_condition)) # 119
length(table(sample.integrated$ct_merge)) # 12

table(sample.integrated$Condition,sample.integrated$ct_merge)

library(future)
plan(multisession, workers = 10) # 使用多线程以加速计算

# 计算细胞类型和条件的数量矩阵
number_matrix <- as.data.frame.matrix(table(sample.integrated$ct_merge, sample.integrated$Condition))
celltypes_filter <- row.names(number_matrix[which(number_matrix$A0 > 2 & number_matrix$A1 > 2 & number_matrix$B1 > 2 & number_matrix$C1 > 2), ])

# 比较信息
VS <- data.frame(
  compare = c('A0vsA1', 'A0vsB1', 'A0vsC1','B0vsB1','C0vsC1'),
  ident1 = c('_A0', '_A0', '_A0','_B0','_C0'),
  ident2 = c('_A1', '_B1', '_C1','_B1','_C1')
)


DEG_mRNA <- data.frame() # 初始化一个空数据框来存储结果

# 对于每一个比较 (VS) 和细胞类型 (celltypes_filter)，找出差异表达的基因
for (j in (1:dim(VS)[1])) {
  eval(parse(text = paste0(format(VS[j,1]), " <- data.frame()"))) # 动态创建数据框
  for (i in (1:length(celltypes_filter))) {
    # 创建有效的变量名称
    valid_name <- make.names(paste0("deg_markers_", celltypes_filter[i], "_", format(VS[j,1]), "_mRNA"))
    
    # 找到两个条件间的差异表达基因
    assign(valid_name, 
           FindMarkers(sample.integrated, ident.1 = paste0(celltypes_filter[i], format(VS[j,2])),
                       ident.2 = paste0(celltypes_filter[i], format(VS[j,3])),
                       logfc.threshold = 0.25, verbose = FALSE))
    
    # 检查数据框是否有结果，并且子集不为空
    eval(parse(text = paste0("if(dim(", valid_name, ")[1] > 0) {
                             ", valid_name, " <- subset(", valid_name, ", p_val_adj < 0.05);
                             if (nrow(", valid_name, ") > 0) {
                               ", valid_name, "$Gene_symbol <- row.names(", valid_name, ");
                               ", valid_name, "$Cell_type <- '", format(celltypes_filter[i]), "';
                               ", valid_name, "$Change <- ifelse(", valid_name, "$avg_log2FC > 0, 'Up', 'Down');
                               ", valid_name, "$Tissue <- 'sample.integrated';
                               ", valid_name, "$Comparison <- '", format(VS[j,1]), "';
                               ", format(VS[j,1]), " <- bind_rows(", format(VS[j,1]), ", ", valid_name, ")
                             }}; 
                             rm(", valid_name, ")")))
  }
  gc() # 手动垃圾收集以节省内存
  eval(parse(text = paste0("DEG_mRNA <- bind_rows(DEG_mRNA, ", format(VS[j,1]), "); rm(", format(VS[j,1]), ")")))
}
write.csv(DEG_mRNA, 'data/1_DEG_mRNA_ct_merge_Condition_Innate.csv')

DEG <- DEG_mRNA

head(DEG)

library(future)
plan(multisession, workers = 10) # 使用多线程以加速计算

# 计算细胞类型和条件的数量矩阵
number_matrix <- as.data.frame.matrix(table(sample.integrated$ct_merge, sample.integrated$Condition))
celltypes_filter <- row.names(number_matrix[which(number_matrix$A0 > 2), ])

# 比较信息
VS <- data.frame(
  compare = c('A0vsA2', 'B0vsB2', 'C0vsC2'),
  ident1 = c('_A0', '_B0', '_C0'),
  ident2 = c('_A2', '_B2', '_C2')
)

DEG_mRNA <- data.frame() # 初始化一个空数据框来存储结果

# 对于每一个比较 (VS) 和细胞类型 (celltypes_filter)，找出差异表达的基因
for (j in 1:nrow(VS)) {
  eval(parse(text = paste0(format(VS[j, 1]), " <- data.frame()"))) # 动态创建数据框
  for (i in 1:length(celltypes_filter)) {
    # 生成组的细胞名称
    cells_1 <- WhichCells(sample.integrated, idents = paste0(celltypes_filter[i], format(VS[j, 2])))
    cells_2 <- WhichCells(sample.integrated, idents = paste0(celltypes_filter[i], format(VS[j, 3])))

    # 检查两组的细胞数量是否都大于等于 3
    if (length(cells_1) >= 3 & length(cells_2) >= 3) {
      # 创建有效的变量名称
      valid_name <- make.names(paste0("deg_markers_", celltypes_filter[i], "_", format(VS[j, 1]), "_mRNA"))
      
      # 找到两个条件间的差异表达基因
      assign(valid_name,
             FindMarkers(
               sample.integrated,
               ident.1 = paste0(celltypes_filter[i], format(VS[j, 2])),
               ident.2 = paste0(celltypes_filter[i], format(VS[j, 3])),
               logfc.threshold = 0.25, verbose = FALSE
             ))
      
      # 检查数据框是否有结果，并且子集不为空
      eval(parse(text = paste0(
        "if(dim(", valid_name, ")[1] > 0) {
           ", valid_name, " <- subset(", valid_name, ", p_val_adj < 0.05);
           if (nrow(", valid_name, ") > 0) {
             ", valid_name, "$Gene_symbol <- row.names(", valid_name, ");
             ", valid_name, "$Cell_type <- '", format(celltypes_filter[i]), "';
             ", valid_name, "$Change <- ifelse(", valid_name, "$avg_log2FC > 0, 'Up', 'Down');
             ", valid_name, "$Tissue <- 'sample.integrated';
             ", valid_name, "$Comparison <- '", format(VS[j, 1]), "';
             ", format(VS[j, 1]), " <- bind_rows(", format(VS[j, 1]), ", ", valid_name, ")
           }}; 
           rm(", valid_name, ")"
      )))
    } else {
      # 如果细胞数量不足，记录空结果
      message("Skipping comparison for ", celltypes_filter[i], " in ", VS[j, 1], " due to insufficient cells.")
    }
  }
  gc() # 手动垃圾收集以节省内存
  eval(parse(text = paste0("DEG_mRNA <- bind_rows(DEG_mRNA, ", format(VS[j, 1]), "); rm(", format(VS[j, 1]), ")")))
}

write.csv(DEG_mRNA, 'data/2_DEG_mRNA_ct_merge_Condition_adaptive.csv')

#up
DEG <- read.csv('data/1_DEG_mRNA_ct_merge_Condition_Innate.csv')
DEG$Comparison[DEG$Comparison == 'A0vsA1'] <- 'A1 vs A0'
DEG$Comparison[DEG$Comparison == 'A0vsB1'] <- 'B1 vs A0'
DEG$Comparison[DEG$Comparison == 'A0vsC1'] <- 'C1 vs A0'
DEG$Comparison[DEG$Comparison == 'B0vsB1'] <- 'B1 vs B0'
DEG$Comparison[DEG$Comparison == 'C0vsC1'] <- 'C1 vs C0'
# 设置因子顺序
DEG$Comparison <- factor(
  DEG$Comparison,
  levels = unique(DEG$Comparison))
DEG <- DEG%>% filter(p_val<= 0.05 & abs(avg_log2FC) > 0.5)%>% filter(Cell_type != 'Platelet')#avg_log2FC<(-0.5)
nrow(DEG)
unique(DEG$Comparison)
table(DEG$Cell_type,DEG$Comparison)
# 计算计数并将其转换为数据框
count_table <- as.data.frame(table(DEG$Cell_type, DEG$Comparison))
colnames(count_table) <- c("Cell_type", "Comparison", "count")
head(count_table)
count_table$Cell_type <- factor(count_table$Cell_type,
                               level=rev(c('Monocyte','NK','CD4T','CD8T','gdT','cDC','ILC','Bcell',
                                       'MAIT','pDC','Tdn','HPSC')))#'Platelet',
unique(count_table$Cell_type)

head(count_table)

# 用于保存所有 count_plot 数据的空数据框
all_count_plot_data_innate <- data.frame()

# 数据处理逻辑
for (i in unique(count_table$Comparison)) {
    # 筛选特定比较的数据
    count_plot <- count_table %>% filter(Comparison == i)
    
    # 添加新的列 "covid_status"，并设置为 "mRNA Vaccine"
    count_plot <- count_plot %>% mutate(covid_status = "mRNA Vaccine")
    
    # 将当前的 count_plot 数据合并到总数据框
    all_count_plot_data_innate <- bind_rows(all_count_plot_data_innate, count_plot)
}

# 保存处理后的数据到 CSV 文件
write.csv(all_count_plot_data_innate, "3_mRNA_all_count_plot_data_innate.csv", row.names = FALSE)
unique(all_count_plot_data_innate$Comparison)
head(all_count_plot_data_innate)

library(grid)
library(RColorBrewer)
pdf("1_mRNA_ct_merge_lollipop_0.05_0.5_innate.pdf", width = 3, height = 2.4)
j = 0
Num <- c('First', 'Second', 'Third','Second', 'Third')

# 获取 "YlOrRd" 调色板的前6个颜色
colors <- RColorBrewer::brewer.pal(4, "YlOrRd")
for (i in unique(count_table$Comparison)) {
    j = j + 1
    count_plot <- count_table %>% filter(Comparison == i)
    p <- ggplot(count_plot, aes(x = Cell_type, y = count)) +
        geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = count), 
                     color = "grey", size = 2) +
        geom_point(aes(color = count), size = 4.5) +
        scale_color_gradientn(name = "Number of DEGs", colors = colors) + # 使用指定的颜色
        geom_text(aes(label = count), color = "Black", size = 4 * 0.4) +
        theme_light() +
        theme(
            text = element_text(size = 9),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 9, colour = "black"),
            plot.title = element_text(size = 9,hjust = 0.5),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            legend.key.width = unit(0.15, "inch"),
            legend.key.height = unit(0.1, "inch"),
            legend.text = element_text(size = 6), 
            legend.title = element_text(size = 6),
            legend.position = c(0.8, 0.2) # 将图例移动到左下角
        ) +
        coord_flip() +
    labs(title=paste0('After ', Num[j], ' Dose (', i, ')     '),
         x="",y='Number of DEGs                    ')
    print(p)
}
dev.off()

#up
DEG <- read.csv('data/2_DEG_mRNA_ct_merge_Condition_adaptive.csv')
DEG$Comparison[DEG$Comparison == 'A0vsA2'] <- 'A2 vs A0'
DEG$Comparison[DEG$Comparison == 'B0vsB2'] <- 'B2 vs B0'
DEG$Comparison[DEG$Comparison == 'C0vsC2'] <- 'C2 vs C0'
# 设置因子顺序
DEG$Comparison <- factor(
  DEG$Comparison,
  levels = unique(DEG$Comparison))
DEG <- DEG%>% filter(p_val<= 0.05 & abs(avg_log2FC) > 0.5)%>% filter(Cell_type != 'Platelet')
nrow(DEG)
unique(DEG$Comparison)
table(DEG$Cell_type,DEG$Comparison)
# 计算计数并将其转换为数据框
count_table <- as.data.frame(table(DEG$Cell_type, DEG$Comparison))
colnames(count_table) <- c("Cell_type", "Comparison", "count")
head(count_table)
count_table$Cell_type <- factor(count_table$Cell_type,
                               level=rev(c('Monocyte','NK','CD4T','CD8T','gdT','cDC','ILC','Bcell',
                                       'MAIT','pDC','Tdn','HPSC')))#'Platelet',
unique(count_table$Cell_type)

# 用于保存所有 count_plot 数据的空数据框
all_count_plot_data <- data.frame()

# 数据处理逻辑
for (i in unique(count_table$Comparison)) {
    # 筛选特定比较的数据
    count_plot <- count_table %>% filter(Comparison == i)
    
    # 添加新的列 "covid_status"，并设置为 "mRNA Vaccine"
    count_plot <- count_plot %>% mutate(covid_status = "mRNA Vaccine")
    
    # 将当前的 count_plot 数据合并到总数据框
    all_count_plot_data <- bind_rows(all_count_plot_data, count_plot)
}

# 保存处理后的数据到 CSV 文件
write.csv(all_count_plot_data, "3_mRNA_all_count_plot_data_adaptive.csv", row.names = FALSE)
unique(all_count_plot_data$Comparison)
head(all_count_plot_data)

library(grid)
library(RColorBrewer)
pdf("2_mRNA_Adaptive_ct_merge_0.05_0.5.pdf", width = 3, height = 2.4)
j = 0
Num <- c('First','Second', 'Third')

# 获取 "YlOrRd" 调色板的前6个颜色
colors <- RColorBrewer::brewer.pal(4, "YlOrRd")
for (i in unique(count_table$Comparison)) {
    j = j + 1
    count_plot <- count_table %>% filter(Comparison == i)
    p <- ggplot(count_plot, aes(x = Cell_type, y = count)) +
        geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = count), 
                     color = "grey", size = 2) +
        geom_point(aes(color = count), size = 4.5) +
        scale_color_gradientn(name = "Number of DEGs", colors = colors) + # 使用指定的颜色
        geom_text(aes(label = count), color = "Black", size = 4 * 0.4) +
        theme_light() +
        theme(
            text = element_text(size = 9),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 9, colour = "black"),
            plot.title = element_text(size = 9,hjust = 0.5),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            legend.key.width = unit(0.15, "inch"),
            legend.key.height = unit(0.1, "inch"),
            legend.text = element_text(size = 6), 
            legend.title = element_text(size = 6),
            legend.position = c(0.9, 0.1) # 将图例移动到左下角
        ) +
        coord_flip() +
    labs(title=paste0('After ', Num[j], ' Dose (', i, ')     '),
         x="",y='Number of DEGs                    ')
    print(p)
}
dev.off()

# pdf("3_ct_merge_lollipop_ggdotchart_UP.pdf", width = 5, height = 18)
# library(ggpubr)
# # 绘制点状图
# ggdotchart(count_table, x='Cell_type', y = 'count',
#            color = "Comparison",                                 # 按照cyl填充颜色
#            palette = c("#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF"), # 修改颜色
#            sorting = "descending",                       # 按降序排序
#            add = "segments",                             # 添加棒子
#            add.params = list(color = "lightgray", size = 1.5), # 改变棒子参数
#            rotate = TRUE,                                # 方向转为垂直
#            group = "Comparison",                                # 按cyl分组
#            dot.size = 6.5,                                 # 改变点的大小
#            label = round(count_table$count),                       # 添加label
#            font.label = list(color = "white", size = 6, vjust = 0.5), # 设置label参数
#            ggtheme = theme_pubr(),                       # 改变主题
#            xlab = "")+                                    # 设置x轴标签
#            facet_wrap(~Comparison, ncol = 1, scales = "free_x")
# dev.off()



# 假设 `all_count_plot_data` 和 `all_count_plot_data_innate` 是两个需要合并的数据框
# 检查列名是否一致
if (all(names(all_count_plot_data) == names(all_count_plot_data_innate))) {
  # 合并两个数据框
  combined_data <- bind_rows(all_count_plot_data, all_count_plot_data_innate)
} else {
  stop("列名不一致，无法合并。请检查两个数据框的列名。")
}

# 查看合并后的数据框
unique(combined_data$Comparison)
table(combined_data$Comparison)
# 保存合并后的数据框到 CSV 文件
write.csv(combined_data, "mRNA_combined_all_count_plot_data.csv", row.names = FALSE)














