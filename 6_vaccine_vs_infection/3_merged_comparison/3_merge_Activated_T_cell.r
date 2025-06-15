library(dplyr)
library(tidyr)
library("tidyverse")



#CVV
dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR/8_Cell2TCR_Activated_cell/1_Tcell'
cvv_res <- read.csv(paste0(dir_path,'/data/All_data_activated_cells_data.csv'))
colnames(cvv_res)

#mRNA
dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/2_All_cell2tcr'
mRNA_res <- read.csv(paste0(dir_path,'/data/Tcell_data_activated_cells_data.csv'))
colnames(mRNA_res)

## challenge
dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/2_All_cell2tcr'
COVID_res <- read.csv(paste0(dir_path,'/data/Tcell_data_activated_cells_Challenge.csv'))
colnames(COVID_res)

#以下处理均删除baseline样本

unique(COVID_res$time_point)

#删除baseline
nrow(COVID_res)
COVID_res <- COVID_res %>% filter(time_point != 'D-1') 
nrow(COVID_res)

COVID_activated_res <- COVID_res %>% filter(activated == 'True')
table(COVID_activated_res$predicted_labels)
table(COVID_res$ct_merge,COVID_res$covid_status)

#CD4T 和CD8T
33756+15776
35332+15266
19125+11648

unique(COVID_activated_res$covid_status)

Sustained_activated_res <- COVID_activated_res %>% filter(covid_status == 'Sustained infection')
Transient_activated_res <- COVID_activated_res %>% filter(covid_status == 'Transient infection')
Abortive_activated_res <- COVID_activated_res %>% filter(covid_status == 'Abortive infection')

sum(Sustained_activated_res$ct_merge %in% c("CD4T", "CD8T"))
sum(Transient_activated_res$ct_merge %in% c("CD4T", "CD8T"))
sum(Abortive_activated_res$ct_merge %in% c("CD4T", "CD8T"))

#删除baseline
nrow(cvv_res)
cvv_res <- cvv_res %>% filter(Condition != 'A0') 
nrow(cvv_res)

cvv_activated_res <- cvv_res %>% filter(activated == 'True')
table(cvv_activated_res$predicted_labels)
table(cvv_activated_res$ct_level_4,cvv_activated_res$predicted_labels)

unique(mRNA_res$group)

#删除baseline
nrow(mRNA_res)
mRNA_res <- mRNA_res %>% filter(group %in% c('Vax','booster'))
nrow(mRNA_res)

table(mRNA_res$celltype.l2)

mRNA_activated_res <- mRNA_res %>% 
      filter(activated == 'True')%>% 
      filter(group %in% c('Vax','booster') )
nrow(mRNA_activated_res)
table(mRNA_activated_res$predicted_labels)
table(mRNA_activated_res$celltype.l2,mRNA_activated_res$predicted_labels)

# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 数据转换函数
process_data <- function(data, group_name) {
  species_table <- as.data.frame(table(data$predicted_labels))
  colnames(species_table) <- c("Predicted_Label", "Count")
  
  species_table <- species_table %>%
    mutate(Percentage = Count / sum(Count) * 100,  # 计算百分比
           Group = group_name)  # 添加分组信息
  return(species_table)
}

# 分别处理三个数据
cvv_data <- process_data(cvv_activated_res, "Inactivated Vaccine (55918 T cells)")
mRNA_data <- process_data(mRNA_activated_res, "mRNA Vaccine (61083 T cells)")
Sustained_data <- process_data(Sustained_activated_res, "Sustained infection (49532 T cells)")
Transient_data <- process_data(Transient_activated_res, "Transient infection (50598 T cells)")
Abortive_data <- process_data(Abortive_activated_res, "Abortive infection (30773 T cells)")
#cv_data <- process_data(cv_activated_res, "CV_19 infection(15983 T cells)")


# 合并三个数据
combined_data <- bind_rows(cvv_data, mRNA_data,Sustained_data, Transient_data,Abortive_data)

# 设置分类因子顺序（确保不同组的分类一致）
combined_data$Predicted_Label <- factor(combined_data$Predicted_Label,
                                        levels = unique(combined_data$Predicted_Label))
combined_data$Group <- factor(combined_data$Group, levels = c("Inactivated Vaccine (55918 T cells)", 
                                                              "mRNA Vaccine (61083 T cells)",
                                                              "Sustained infection (49532 T cells)",
                                                              "Transient infection (50598 T cells)",
                                                              "Abortive infection (30773 T cells)"))
combined_data$Predicted_Label <- factor(combined_data$Predicted_Label,
                                       level=c('T Reg Activated','T Reg Activated Cycling',
                                               'T CD4 Activated CTL','T CD4 Activated CTL Cycling',
                                               'T CD4 Activated Helper 1','T CD4 Activated Helper 1 Cycling',
                                               'T CD8 Activated CTL','T CD8 Activated CTL Cycling'
                                               ))
unique(combined_data$Predicted_Label)
# 设置颜色（每个组一个颜色）
colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(combined_data$Predicted_Label)))

write.csv(combined_data, "data/1_combined_5_data_activated_Tcell.csv")

# 加载RColorBrewer包
library(RColorBrewer)

combined_data <- read.csv("data/1_combined_5_data_activated_Tcell.csv")

# 设置分类因子顺序（确保不同组的分类一致）
combined_data$Predicted_Label <- factor(combined_data$Predicted_Label,
                                        levels = unique(combined_data$Predicted_Label))
combined_data$Group <- factor(combined_data$Group, levels = c("Inactivated Vaccine (55918 T cells)", 
                                                              "mRNA Vaccine (61083 T cells)",
                                                              "Sustained infection (49532 T cells)",
                                                              "Transient infection (50598 T cells)",
                                                              "Abortive infection (30773 T cells)"))
combined_data$Predicted_Label <- factor(combined_data$Predicted_Label,
                                       level=c('T Reg Activated','T Reg Activated Cycling',
                                               'T CD4 Activated CTL','T CD4 Activated CTL Cycling',
                                               'T CD4 Activated Helper 1','T CD4 Activated Helper 1 Cycling',
                                               'T CD8 Activated CTL','T CD8 Activated CTL Cycling'
                                               ))

# 绘制分组条形图
pdf("figures/3_Combined_Predicted_Label_Bar_Plot_5_data.pdf", 6.5, 3)  

ggplot(combined_data, aes(x = Predicted_Label, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.95, color = "white") +  # 分组条形图
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.2, size = 1.5) +  # 显示数字
  scale_fill_manual(values = brewer.pal(5, "Set2")) +  # 使用调色板设置颜色
  theme_minimal() +
  labs(fill = "Group", title = "Predicted Labels Distribution Across Groups",
       x = "Predicted Label", y = "Count") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 8 ,angle = 45, hjust = 1,  color = "black"),  # 调整 X 轴文本
    axis.text.y = element_text(size = 8 ,color = "black"),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
       legend.key.size = unit(0.2, "inch"),
    plot.title = element_text(size = 8, hjust = 0.5, vjust = 1, face = "bold")  # 总图标题样式
  )

dev.off()

# 分组总数量（括号中的值）

total_counts <-c(
    "Inactivated Vaccine (55918 T cells)"= 55918,
  "mRNA Vaccine (61083 T cells)"= 61083, 
  "Sustained infection (49532 T cells)"= 49532, 
  "Transient infection (50598 T cells)"= 50598, 
  "Abortive infection (30773 T cells)"= 30773
)
# 添加总数量列并计算占比
combined_data <- combined_data %>%
  mutate(
    Total_Count = total_counts[Group],  # 对应总数量
    Percentage = Count / Total_Count * 100  # 计算占比（百分比）
  )
head(combined_data )

# 绘制分组占比条形图
pdf("figures/3_Combined_Predicted_Label_Percentage_Plot.pdf",  8, 3.8)

ggplot(combined_data, aes(x = Predicted_Label, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.95, color = "white") +  # 分组条形图
  geom_text(
    aes(label = sprintf("%.1f", Percentage)),  # 显示百分比（保留1位小数）
    position = position_dodge(width = 0.8),
    vjust = -0.2,
    size = 2.3
  ) +
  scale_fill_manual(values = brewer.pal(5, "Set2")) +  # 使用调色板设置颜色
  theme_minimal() +
  labs(
    fill = "Group",
    title = "Predicted Labels Proportion Across Groups",
    x = "Predicted Label",
    y = "Percentage (%)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 8 ,angle = 45, hjust = 1,  color = "black"),  # 调整 X 轴文本
    axis.text.y = element_text(size = 8 ,color = "black"),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
       legend.key.size = unit(0.2, "inch"),
    plot.title = element_text(size = 8, hjust = 0.5, vjust = 1, face = "bold")  # 总图标题样式
  )

dev.off()

# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 设置固定颜色映射，确保所有饼图中颜色一致
all_labels <- c('T Reg Activated','T Reg Activated Cycling',
                                               'T CD4 Activated CTL','T CD4 Activated CTL Cycling',
                                               'T CD4 Activated Helper 1','T CD4 Activated Helper 1 Cycling',
                                               'T CD8 Activated CTL','T CD8 Activated CTL Cycling'
                                               )
colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(all_labels)), all_labels)
colors
head(combined_data)

combined_data$Group <- factor(combined_data$Group, levels = c("Sustained infection (49532 T cells)",
                                                              "Transient infection (50598 T cells)",
                                                          "Abortive infection (30773 T cells)",
                                                          "Inactivated Vaccine (55918 T cells)", 
                                                              "mRNA Vaccine (61083 T cells)"))

# 加载必要的库
library(ggplot2)
library(dplyr)

# 数据预处理和饼图绘制
plot_faceted_pie_chart <- function(data, output_file) {
  # 计算百分比
  data <- data %>%
    group_by(Group) %>%
    mutate(Percentage = Count / sum(Count) * 100)
  
  # 绘制带分面的小饼图
  p <- ggplot(data, aes(x = "", y = Percentage, fill = Predicted_Label)) +
    geom_bar(stat = "identity", width = 0.9, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors) +  # 设置颜色映射
    theme_void() +  # 移除背景
    labs(fill = "Predicted Label") +  # 设置图例标题
    facet_wrap(~ Group, ncol = 3) +  # 按分组分面，指定每行最多显示3个
    theme(
      strip.text = element_text(size = 14, face = "bold"),  # 分面标签样式
      legend.position = "right",  # 图例放在右侧
        plot.margin = unit(c(2, 2, 2, 2), "char"),  # 图的外边距
      legend.title = element_text(size = 20),  # 图例标题样式
      legend.text = element_text(size = 16)  # 图例文字大小
    )
  
  # 输出为 PDF
  ggsave(output_file, plot = p, width = 14, height = 10)
}

# 调用函数
plot_faceted_pie_chart(
  combined_data,
  "figures/4_Faceted_Predicted_Labels_Pie_Charts.pdf"
)


