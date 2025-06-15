library(tidyverse)
library(ggplot2)
require(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

DEG_infection <- read.csv('data/all_count_plots_with_covid_status.csv')
unique(DEG_infection$Comparison)
head(DEG_infection)

DEG_mRNA_VIA <- read.csv('data/all_count_plot_data_with_mRNA_Vaccine.csv')
unique(DEG_mRNA_VIA$Comparison)
head(DEG_mRNA_VIA)

DEG_cvv <- read.csv('data/CVV_combined_all_count_plot_data.csv')
unique(DEG_cvv$Comparison)
head(DEG_cvv)

DEG_mRNA<- read.csv('data/mRNA_combined_all_count_plot_data.csv')
unique(DEG_mRNA$Comparison)
head(DEG_mRNA)





DEG_mRNA_inn <- DEG_mRNA  %>% filter(Comparison %in% c('A1 vs A0','B1 vs B0','C1 vs C0'))
DEG_cvv_inn <- DEG_cvv  %>% filter(Comparison %in% c('A1 vs A0','B1 vs B0','C1 vs C0'))
DEG_merge_inn <- bind_rows(DEG_infection, DEG_mRNA_inn,DEG_cvv_inn)
DEG_merge_inn <- DEG_merge_inn %>%
  mutate(Comparison_covid_status = paste(Comparison, covid_status, sep = "_"))
unique(DEG_merge_inn$Comparison_covid_status)
unique(DEG_merge_inn$Cell_type)

DEG_merge_inn$abbreviation <- DEG_merge_inn$Comparison_covid_status
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="D3 vs D-1_Sustained infection"] <- "Sustained"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="D3 vs D-1_Transient infection"] <- "Transient"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="D3 vs D-1_Abortive infection"] <- "Abortive"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="A1 vs A0_mRNA Vaccine"] <- "First"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="B1 vs B0_mRNA Vaccine"] <- "Second"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="C1 vs C0_mRNA Vaccine"] <- "Third"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="A1 vs A0_Inactivated Vaccine"] <- "I_First"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="B1 vs B0_Inactivated Vaccine"] <- "I_Second"
DEG_merge_inn$abbreviation[DEG_merge_inn$abbreviation=="C1 vs C0_Inactivated Vaccine"] <- "I_Third"
unique(DEG_merge_inn$abbreviation)
head(DEG_merge_inn)
write.csv(DEG_merge_inn,'data/Innate_Immunity_DEG_merge.csv')

# 计算计数并将其转换为数据框
DEG_merge_complete <- DEG_merge_inn
DEG_merge_complete$Cell_type <- factor(DEG_merge_complete$Cell_type,
                               level=rev(c('Monocyte','NK','CD4T','CD8T','gdT','cDC','ILC','Bcell',
                                       'MAIT','pDC','Tdn','HPSC')))
DEG_merge_complete$abbreviation <- factor(DEG_merge_complete$abbreviation,
                               level=c('I_First','I_Second','I_Third',
                                       'First','Second','Third',
                                       'Sustained','Transient','Abortive'))
# 确保数据完整，填补缺失细胞的 count 为 0
DEG_merge_complete <- DEG_merge_complete%>% 
  complete(Cell_type, abbreviation, fill = list(count = 0))
table(DEG_merge_complete$Cell_type,DEG_merge_complete$abbreviation)

pdf("figures/2_Innate_Immunity_ct_merge_lollipop.pdf", width = 5, height = 4)
# 获取渐变色
colors <- RColorBrewer::brewer.pal(4, "YlOrRd")
# 绘制图形
ggplot(DEG_merge_complete, aes(x = abbreviation, y = Cell_type)) +
  geom_point(aes(color = count), size=5,show.legend = TRUE) +  # 点的大小和颜色根据 count 变化
  geom_text(aes(label = count), color = "black", size = 1.8, vjust = 0.5) +  # 显示 count 值
  scale_color_gradientn(name = "Number of DEGs", colors = colors) +  # 渐变色
  theme_light() +
  theme(
    text = element_text(size = 6,hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45 , hjust = 1,colour = "black"),  # X 轴标签角度调整
    axis.text.y = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 12,hjust = 0.5),#,vjust = -75
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    legend.key.width = unit(0.2, "inch"),
    legend.key.height = unit(0.15, "inch"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "right"  # 图例位置
  ) +
  xlab("") +
  ylab("") +
  labs(title = 'Number of DEGs ')
dev.off()

library(dplyr)
library(ggplot2)

# 确保每个 abbreviation 内的 count 按照数量大小计算颜色
DEG_merge_complete <- DEG_merge_complete %>%
  group_by(abbreviation) %>% 
  mutate(color_group = scales::rescale(count, to = c(0, 1))) %>% # 归一化 count 值到 [0, 1]
  ungroup()

# 创建颜色梯度
colors <- RColorBrewer::brewer.pal(3, "YlOrRd")

# 绘制图形
pdf("figures/2_Innate_Immunity_ct_merge_lollipop_Normalized_2.pdf", width = 5, height = 4)
ggplot(DEG_merge_complete, aes(x = abbreviation, y = Cell_type)) +
  geom_point(aes(color = color_group), size = 5, show.legend = TRUE) +  # 点的大小和颜色根据归一化 count 变化
  geom_text(aes(label = count), color = "black", size = 1.8, vjust = 0.5) +  # 显示 count 值
  #scale_color_gradientn(name = "Normalized Count", colors = colors) +  # 渐变色
  scale_color_gradientn(
        name = "Normalized Count", 
        colors = colors, 
        limits = c(0.4, 1),  # 设置渐变色范围
        oob = scales::squish  # 将超出范围的值压缩
      ) +  # 渐变色
  theme_light() +
  theme(
    text = element_text(size = 6, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, colour = "black"),  # X 轴标签角度调整
    axis.text.y = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    legend.key.width = unit(0.2, "inch"),
    legend.key.height = unit(0.15, "inch"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "right"  # 图例位置
  ) +
  xlab("") +
  ylab("") +
  labs(title = 'Number of DEGs ')
dev.off()

unique(DEG_mRNA$Comparison)

DEG_mRNA_ada <- DEG_mRNA %>% filter(Comparison %in% c('A2 vs A0','B2 vs B0','C2 vs C0'))
DEG_cvv_ada <- DEG_cvv  %>% filter(Comparison %in% c('A2 vs A0','B2 vs B0','C2 vs C0'))
DEG_merge_ada <- bind_rows(DEG_cvv_ada, DEG_mRNA_ada)
DEG_merge_ada <- DEG_merge_ada %>%
  mutate(Comparison_covid_status = paste(Comparison, covid_status, sep = "_"))
unique(DEG_merge_ada$Comparison_covid_status)
unique(DEG_merge_ada$Cell_type)

DEG_merge_ada$abbreviation <- DEG_merge_ada$Comparison_covid_status
DEG_merge_ada$abbreviation[DEG_merge_ada$abbreviation=="A2 vs A0_Inactivated Vaccine"] <- "I_First"
DEG_merge_ada$abbreviation[DEG_merge_ada$abbreviation=="B2 vs B0_Inactivated Vaccine"] <- "I_Second"
DEG_merge_ada$abbreviation[DEG_merge_ada$abbreviation=="C2 vs C0_Inactivated Vaccine"] <- "I_Third"
DEG_merge_ada$abbreviation[DEG_merge_ada$abbreviation=="A2 vs A0_mRNA Vaccine"] <- "First"
DEG_merge_ada$abbreviation[DEG_merge_ada$abbreviation=="B2 vs B0_mRNA Vaccine"] <- "Second"
DEG_merge_ada$abbreviation[DEG_merge_ada$abbreviation=="C2 vs C0_mRNA Vaccine"] <- "Third"
unique(DEG_merge_ada$abbreviation)
head(DEG_merge_ada)
write.csv(DEG_merge_ada,'data/Adaptive_Immunity_DEG_merge.csv')

# 计算计数并将其转换为数据框
DEG_merge_complete <- DEG_merge_ada
DEG_merge_complete$Cell_type <- factor(DEG_merge_complete$Cell_type,
                               level=rev(c('Monocyte','NK','CD4T','CD8T','gdT','cDC','ILC','Bcell',
                                       'MAIT','pDC','Tdn','HPSC')))
DEG_merge_complete$abbreviation <- factor(DEG_merge_complete$abbreviation,
                               level=c('I_First','I_Second','I_Third',
                                       'First','Second','Third'))
# 确保数据完整，填补缺失细胞的 count 为 0
DEG_merge_complete <- DEG_merge_complete%>% 
  complete(Cell_type, abbreviation, fill = list(count = 0))
table(DEG_merge_complete$Cell_type,DEG_merge_complete$abbreviation)

pdf("figures/1_Adaptive_Immunity_ct_merge_lollipop.pdf", width = 4.2, height = 3.6)
# 获取渐变色
colors <- RColorBrewer::brewer.pal(5, "YlOrRd")
# 绘制图形
ggplot(DEG_merge_complete, aes(x = abbreviation, y = Cell_type)) +
  geom_point(aes(color = count), size=5,show.legend = TRUE) +  # 点的大小和颜色根据 count 变化
  geom_text(aes(label = count), color = "black", size = 1.8, vjust = 0.5) +  # 显示 count 值
  scale_color_gradientn(name = "Number of DEGs", colors = colors) +  # 渐变色
  theme_light() +
  theme(
    text = element_text(size = 6,hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45 , hjust = 1,colour = "black"),  # X 轴标签角度调整
    axis.text.y = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 12,hjust = 0.5),#,vjust = -75
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    legend.key.width = unit(0.2, "inch"),
    legend.key.height = unit(0.15, "inch"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "right"  # 图例位置
  ) +
  xlab("") +
  ylab("") +
  labs(title = 'Number of DEGs ')
dev.off()

library(dplyr)
library(ggplot2)

# 确保每个 abbreviation 内的 count 按照数量大小计算颜色
DEG_merge_complete <- DEG_merge_complete %>%
  group_by(abbreviation) %>% 
  mutate(color_group = scales::rescale(count, to = c(0, 1))) %>% # 归一化 count 值到 [0, 1]
  ungroup()


# 创建颜色梯度
colors <- RColorBrewer::brewer.pal(3, "YlOrRd")

# 绘制图形
pdf("figures/1_Adaptive_Immunity_ct_merge_lollipop_Normalized_2.pdf", width = 4.2, height = 3.6)
ggplot(DEG_merge_complete, aes(x = abbreviation, y = Cell_type)) +
  geom_point(aes(color = color_group), size = 5, show.legend = TRUE) +  # 点的大小和颜色根据归一化 count 变化
  geom_text(aes(label = count), color = "black", size = 1.8, vjust = 0.5) +  # 显示 count 值
  #scale_color_gradientn(name = "Normalized Count", colors = colors) +  # 渐变色
    scale_color_gradientn(
        name = "Normalized Count", 
        colors = colors, 
        limits = c(0.4, 1),  # 设置渐变色范围
        oob = scales::squish  # 将超出范围的值压缩
      ) +  # 渐变色
  theme_light() +
  theme(
    text = element_text(size = 6, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, colour = "black"),  # X 轴标签角度调整
    axis.text.y = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    legend.key.width = unit(0.2, "inch"),
    legend.key.height = unit(0.15, "inch"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "right"  # 图例位置
  ) +
  xlab("") +
  ylab("") +
  labs(title = 'Number of DEGs ')
dev.off()


