library(tidyverse)
library(grid)
library(ggplot2)
library(scales)
require(RColorBrewer)

#celltypel2
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
  Platelet="#BD956A",
  Tdn="#495608"# 其他细胞类型
)

#up
DEG <- read.csv('data/2_DEG_challenge_ct_merge_Condition_covid_status.csv')

head(DEG,n=2)
unique(DEG$Cell_type)

unique(DEG$Comparison)

unique(DEG$covid_status)

library(dplyr)
library(ggplot2)
library(grid)
library(RColorBrewer)

# 初始化一个空的数据框，用于存储所有子集的 count_plot 数据
all_count_plots <- data.frame()

# 按照 'covid_status' 切割 DEG 数据并对每个子集进行计算和绘图
for (status in unique(DEG$covid_status)) {
  # 过滤出特定 covid_status 的数据
  DEG_subset <- DEG %>% filter(covid_status == status)
  
  # 替换 'Comparison' 中的 'D3vsD-1' 为 'D3 vs D-1'
  DEG_subset$Comparison[DEG_subset$Comparison == 'D3vsD-1'] <- 'D3 vs D-1'
  
  # 设置因子顺序
  DEG_subset$Cell_type <- factor(DEG_subset$Cell_type, levels = unique(DEG_subset$Cell_type))
  DEG_subset$Comparison <- factor(DEG_subset$Comparison, levels = unique(DEG_subset$Comparison))
  
  # 筛选条件 p_val_adj <= 0.05 且 avg_log2FC < -0.5
  DEG_filtered <- DEG_subset %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) > 0.5)
  
  # 查看筛选后的数据行数
  print(paste("covid_status:", status, "Number of rows:", nrow(DEG_filtered)))
  
  # 创建 Cell_type 和 Comparison 的频次表，并将其转换为数据框
  count_table <- as.data.frame(table(DEG_filtered$Cell_type, DEG_filtered$Comparison))
  colnames(count_table) <- c("Cell_type", "Comparison", "count")
  
  # 增加 covid_status 列
  count_table$covid_status <- status
  
  # 设置 Cell_type 的顺序
  count_table$Cell_type <- factor(count_table$Cell_type, 
                                  levels = rev(c('Monocyte', 'NK', 'CD4T', 'CD8T', 'gdT', 'cDC', 'ILC', 'Bcell',
                                                 'MAIT', 'pDC', 'Tdn', 'HPSC')))
  
  # 将当前子集的 count_table 添加到所有结果的汇总表中
  all_count_plots <- rbind(all_count_plots, count_table)
  
  # 设置调色板
  colors <- brewer.pal(4, "YlOrRd")
  Num <- c('First', 'Second', 'Third', 'Second', 'Third')
  
  # 绘图并保存为 PDF
  pdf(paste0(status, "_lollipop_chart_UP_0.05_0.5.pdf"), width = 3, height = 3)
  j <- 0
  
  for (i in unique(count_table$Comparison)) {
    j <- j + 1
    count_plot <- count_table %>% filter(Comparison == i)
    
    p <- ggplot(count_plot, aes(x = Cell_type, y = count)) +
      geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = count), 
                   color = "grey", size = 2) +
      geom_point(aes(color = count), size = 4.5) +
      scale_color_gradientn(name = "Number of DEGs", colors = colors) + # 使用指定的颜色
      geom_text(aes(label = count), color = "Black", size = 4 * 0.4) +
      theme_light() +
      theme(
        text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        legend.key.width = unit(0.04, "inch"),
        legend.key.height = unit(0.05, "inch"),
        legend.text = element_text(size = 4), 
        legend.title = element_text(size = 4),
        legend.position = c(0.9, 0.4) # 将图例移动到左下角
      ) +
      coord_flip() +
      xlab("") +
      labs(title = paste0('After ', Num[j], ' Dose (', i, ', ', status, ')            '),
           x = "", y = 'Number of DEGs                          ')
    
    print(p)
  }
  
  dev.off()
}

# 保存所有子集的合并结果为 CSV 文件
write.csv(all_count_plots, "data/all_count_plots_with_covid_status.csv", row.names = FALSE)




