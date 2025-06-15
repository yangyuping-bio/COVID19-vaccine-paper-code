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
  #Platelet="#BD956A",
  Tdn="#495608"# 其他细胞类型
)

RdBu_new <- c('#2166AC','#4393C3','#92C5DE','#F4A582','#D6604D','#B2182B')

challenge_data <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/3_Challenge/data/challenge_adata_obs.csv")
colnames(challenge_data)
cvv_data <-read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/4_CVV/1_IFN_HLA_DQA2/data/1_meta_data.df.csv")
colnames(cvv_data)
VIA_data <-read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/2_VIA_paper_data/data/1_VIA_meta_data.df.csv")
colnames(VIA_data)
mRNA_data <-read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/4_new_celltype/data/1_mRNA_meta_data_df.csv")
colnames(mRNA_data)

unique(VIA_data$celltypel2)
unique(challenge_data$cell_type)
unique(cvv_data$celltypel3)
unique(mRNA_data$celltype.l2)

mRNA_data$ct_merge.l2 <- mRNA_data$celltype.l2
cvv_data$ct_merge.l2 <- cvv_data$celltypel3
challenge_data$ct_merge.l2 <- challenge_data$cell_type
VIA_data$ct_merge.l2 <- VIA_data$celltypel2

challenge_data$ct_merge.l2[challenge_data$ct_merge.l2 %in% c('T Reg')] <- "CD4_Treg"
mRNA_data$ct_merge.l2[mRNA_data$ct_merge.l2 %in% c('Treg')] <- "CD4_Treg"
VIA_data$ct_merge.l2[VIA_data$ct_merge.l2 %in% c('Treg')] <- "CD4_Treg"

challenge_data$ct_merge.l2[challenge_data$ct_merge.l2 %in% c('Plasma Cell')] <- "Plasma"
mRNA_data$ct_merge.l2[mRNA_data$ct_merge.l2 %in% c('Plasmablasts')] <- "Plasma"
VIA_data$ct_merge.l2[VIA_data$ct_merge.l2 %in% c('Plasmablast')] <- "Plasma"

mRNA_data$HLA_score <- mRNA_data$HLA_DQA21
#mRNA_data <- mRNA_data[ ,colnames(mRNA_data)[c(15:23,100:110)]]
colnames(mRNA_data)
unique(mRNA_data$ct_merge)
mRNA_split_data <- mRNA_data%>% filter(state =="mRNA")
cv_19_split_data <- mRNA_data%>% filter(state =="cv_19")
mRNA_split_data$Condition <- factor(
  mRNA_split_data$Condition,
  levels = c("A0" ,"A1","A2" ,"B0" ,"B1","B2" ,"C0" ,"C1","C2","C3" ))
mRNA_split_data$covid_status <- "mRNA Vaccine"
cv_19_split_data$Condition <- factor(
  cv_19_split_data$Condition,
  levels = c('convalescent','acute'))
cv_19_split_data$covid_status <- "COVID_19"

challenge_data$time_point <- factor(challenge_data$time_point,
                                     level=c('D-1','D3','D7','D10','D14','D28'))
challenge_data$covid_status <- factor(challenge_data$covid_status,
                                     level=c('Sustained infection','Transient infection','Abortive infection'))
challenge_data$ct_merge <- challenge_data$cell_compartment
challenge_data$ct_merge[challenge_data$ct_merge %in% c('T CD4','T Reg')] <- "CD4T"
challenge_data$ct_merge[challenge_data$ct_merge %in% c('T CD8')] <- "CD8T"
challenge_data$ct_merge[challenge_data$ct_merge %in% c('T G/D')] <- "gdT"
challenge_data$ct_merge[challenge_data$ct_merge %in% c('T MAI')] <- "MAIT"
challenge_data$ct_merge[challenge_data$ct_merge %in% c('T DN')] <- "Tdn"
challenge_data$ct_merge[challenge_data$ct_merge %in% c('NK','NK CD56+')] <- "NK"
challenge_data$ct_merge[challenge_data$ct_merge %in% c('B','B ABS')] <- "Bcell"
challenge_data$ct_merge[challenge_data$ct_merge %in% c('HPC')] <- "HPSC"
challenge_data$ct_merge[challenge_data$cell_type %in% c('pDC')] <- "pDC"
challenge_data$ct_merge[challenge_data$cell_type %in% c('cDC1','cDC2','cDC3','AS-DC')] <- "cDC"
unique(challenge_data$ct_merge)
challenge_data$Condition <- challenge_data$time_point
#challenge_data$ifn_stim1 <- challenge_data$ifn_stim
challenge_data$HLA_score <- challenge_data$HLA_DQA21

cvv_data$Condition <- factor(
  cvv_data$Condition,
  levels = c("A0" ,"A1","A2" ,"B0" ,"B1","B2" ,"C0" ,"C1","C2" ))
cvv_data$ct_merge <- cvv_data$celltypel2
#cvv_data$ct_merge[cvv_data$ct_merge %in% c('Platelet')] <- "Other"
unique(cvv_data$ct_merge)
cvv_data$covid_status <- "Inactivated Vaccine"

VIA_data$timepoint <- factor(VIA_data$timepoint,
                                     level=c('Day0','Day2','Day11','Day28'))
VIA_data$ct_merge <- VIA_data$celltypel1
VIA_data$ct_merge[VIA_data$ct_merge %in% c('B')] <- "Bcell"
VIA_data$ct_merge[VIA_data$ct_merge %in% c('CD4 T')] <- "CD4T"
VIA_data$ct_merge[VIA_data$ct_merge %in% c('CD8 T')] <- "CD8T"
VIA_data$ct_merge[VIA_data$ct_merge %in% c('Mono')] <- "Monocyte"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c("group_A","group_B")] <- "CD8T"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('MAIT')] <- "MAIT"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('dnT')] <- "Tdn"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('gdT')] <- "gdT"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('ILC')] <- "ILC"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('HSPC')] <- 'HPSC'
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('cDC1','cDC2','ASDC')] <- "cDC"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('pDC')] <- "pDC"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('Platelet')] <- "Platelet"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('NK Proliferating')] <- "NK"
VIA_data$ct_merge[VIA_data$celltypel2 %in% c('CD4 Proliferating')] <- "CD4T"
unique(VIA_data$ct_merge)
VIA_data$Condition <- VIA_data$timepoint
VIA_data$covid_status <- "mRNA Vaccine(VIA)"

# colnames(VIA_data)
# colnames(cvv_data)
colnames(challenge_data)
# colnames(mRNA_split_data)
# colnames(cv_19_split_data)

VIA_data_merge <- VIA_data[,c('ct_merge','Condition','ifn_stim1','HLA_score','covid_status','ct_merge.l2','VIA1','Apoptosis1','ATG71','mac_ATG71')]
cvv_data_merge <- cvv_data[,c('ct_merge','Condition','ifn_stim1','HLA_score','covid_status','ct_merge.l2','VIA1','Apoptosis1','ATG71','mac_ATG71')]
challenge_data_merge <- challenge_data[,c('ct_merge','Condition','ifn_stim1','HLA_score','covid_status','ct_merge.l2','VIA1','Apoptosis1','ATG71','mac_ATG71')]
mRNA_split_data <- mRNA_split_data[,c('ct_merge','Condition','ifn_stim1','HLA_score','covid_status','ct_merge.l2','VIA1','Apoptosis1','ATG71','mac_ATG71')]
cv_19_split_data <- cv_19_split_data[,c('ct_merge','Condition','ifn_stim1','HLA_score','covid_status','ct_merge.l2','VIA1','Apoptosis1','ATG71','mac_ATG71')]
merged_data <- rbind(mRNA_split_data, cvv_data_merge, challenge_data_merge)
unique(merged_data$covid_status)
head(merged_data)

unique(VIA_data_merge$ct_merge)
unique(cvv_data_merge$ct_merge)
unique(mRNA_split_data$ct_merge)

unique(merged_data$Condition)
# Condition_seq <- c('Day0','Day2','Day11','Day28','A0','A1','A2','B0','B1','B2','C0','C1','C2','D-1','D3','D7','D10','D14','D28')
# merged_data$Condition	<- factor(merged_data$Condition,level=Condition_seq)
#head(merged_data)

write.csv(merged_data,"data/1_Merged_data_CVV_VIA_Challenge_with_Other.csv")

#merged_data <- read.csv("data/1_Merged_data_CVV_VIA_Challenge_with_Other.csv")

merged_data <- merged_data %>%
  filter(! ct_merge %in% c('Platelet')) #三个数据统一有的细胞类型
merged_data$ct_merge <- factor(merged_data$ct_merge,
                               levels = rev(c('CD4T','CD8T','MAIT','gdT','Tdn',
                                              'NK','Bcell','Monocyte','cDC','pDC',
                                              'ILC','HPSC')))
unique(merged_data$ct_merge)
merged_data$covid_status <- factor(merged_data$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine','Sustained infection','Transient infection','Abortive infection'))
unique(merged_data$covid_status)
write.csv(merged_data,"data/1_Merged_data_CVV_VIA_Challenge.csv")

head(merged_data)







# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(ifn_stim1, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(ifn_stim1 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )
average_scores$covid_status1 <- as.character(average_scores$covid_status)
average_scores$covid_status1[average_scores$covid_status1 %in% c('Sustained infection',
                                                                 'Transient infection',
                                                                 'Abortive infection')] <- 'Infection'
average_scores <- average_scores %>%
  group_by(covid_status1) %>% # 按照 covid_status1 分组
  mutate(mean_score_normalized = rescale(mean_score, to = c(-1, 1))) # 在每个 covid_status1 内归一化
head(average_scores)

# 绘制热图
pdf("figures/2_Merge_ct_merge_IFN_module_Dotplot.pdf", width = 4.5, height = 1.8)
ggplot(average_scores, aes(x = Condition, y = ct_merge)) +
  geom_point(aes(size = percentage_expressing, color = mean_score_normalized)) +
  scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
   #scale_color_gradientn(colors = brewer.pal(9, "Greens")) +
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    panel.border = element_blank(), # 去掉边框
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5), # 分图标题字体设置为5，居中
    text = element_text(size = 0),
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    #axis.text.x = element_text(size = 4), # 设置 X 轴文本
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.05, vjust = 0.5),
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个legend之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放
dev.off()

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
average_scores$covid_status <- factor(average_scores$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/2_Merge_ct_merge_IFN_Dotplot_with_background.pdf", width = 4.5, height = 1.68)

ggplot(average_scores, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = mean_score_normalized)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new) +
  scale_fill_manual(values = c("red", "blue", "#3A6963")) + # 设置矩形颜色
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(ifn_stim1, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(ifn_stim1 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- average_scores %>%
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_score = mean(mean_score, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>% 
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = mean_score / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
final_merged_data$covid_status <- factor(final_merged_data$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/2_Merge_ct_merge_IFN_Dotplot_with_background_FoldChange_baseline.pdf", width = 5, height = 1.68)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new, name = "Fold change over baseline",
                       limits = c(0, 2),#limits = c(0.5, 1.5),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  scale_fill_manual(values = c("red", "blue", "#3A6963"))+
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression(Fold change over baseline)")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(HLA_score, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(HLA_score > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )
average_scores$covid_status1 <- as.character(average_scores$covid_status)
average_scores$covid_status1[average_scores$covid_status1 %in% c('Sustained infection',
                                                                 'Transient infection',
                                                                 'Abortive infection')] <- 'Infection'
average_scores <- average_scores %>%
  group_by(covid_status1) %>% # 按照 covid_status1 分组
  mutate(mean_score_normalized = rescale(mean_score, to = c(-1, 1))) # 在每个 covid_status1 标准化，Challenge是0-3
head(average_scores)

# 绘制热图
pdf("figures/3_Merge_ct_merge_HLA_score_Dotplot.pdf", width = 4.5, height = 1.5)
ggplot(average_scores, aes(x = Condition, y = ct_merge)) +
  geom_point(aes(size = percentage_expressing, color = mean_score_normalized)) +
  scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
   #scale_color_gradientn(colors = brewer.pal(9, "Greens")) +
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    panel.border = element_blank(), # 去掉边框
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5), # 分图标题字体设置为5，居中
    text = element_text(size = 0),
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    #axis.text.x = element_text(size = 4), # 设置 X 轴文本
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5),
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个legend之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放
dev.off()

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
average_scores$covid_status <- factor(average_scores$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/3_Merge_ct_merge_HLA_score_Dotplot_with_background.pdf", width = 4.5, height = 1.68)

ggplot(average_scores, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = mean_score_normalized)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new) +
  scale_fill_manual(values = c("red", "blue", "#3A6963")) + # 设置矩形颜色
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(HLA_score, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(HLA_score > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- average_scores %>%
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_score = mean(mean_score, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>% 
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = mean_score / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
final_merged_data$covid_status <- factor(final_merged_data$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/3_Merge_ct_merge_HLA_score_Dotplot_with_background_FoldChange_baseline.pdf", width = 5, height = 1.68)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new, name = "Fold change over baseline",
                       limits = c(0, 2),#limits = c(0.5, 1.5),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  scale_fill_manual(values = c("red", "blue", "#3A6963"))+
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression(Fold change over baseline)")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()



#celltypel2
colorl2 <- c(
  Bcell="#AB3282",  # B细胞类型
  CD4T="#6778AE", # CD4 T细胞类型
  CD8T="#53A85F", # CD8 T细胞类型
  cDC="#9FA3A8",  # CDC细胞类型
  gdT="#3A6963",  #gdT
  HPSC="#585658",#HPSC
  ILC="#F1BC62", #ILC
  MAIT="#5F3D69",#MAIT
  Monocyte="#E95C39",# Mono细胞类型 
  NK="#E1A111", # NK细胞类型
  pDC="#968175", 
 # Platelet="#BD956A",
  Tdn="#585658"# 其他细胞类型
)

merged_data <- merged_data %>%
  group_by(covid_status,Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(covid_status,Condition, ct_merge) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition) 
head(merged_data)

#celltypel2
pdf("figures/4_Line_Condition_ct_merge_free_x.pdf",  30, 16)
 ggplot(merged_data, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(ct_merge))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl2) +
  facet_wrap(covid_status~ ct_merge, scales = "free_x", ncol = 12) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()

merged_data_filter <- merged_data%>% filter(ct_merge %in% c('CD4T','CD8T','NK','Bcell','Monocyte'))
#celltypel2
pdf("figures/4_Line_Condition_ct_merge_filter_free.pdf", 15, 12)
 ggplot(merged_data_filter, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(ct_merge))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl2) +
  facet_wrap(covid_status~ ct_merge, scales = "free", ncol = 5) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()

merged_data_2 <- rbind(mRNA_split_data, VIA_data_merge,cvv_data_merge, challenge_data_merge,cv_19_split_data)
unique(merged_data_2$covid_status)


colorl3=c("Inactivated Vaccine"="#ffd401",
          "mRNA Vaccine"="#00BFFF",
          'mRNA Vaccine(VIA)'="#00BFFF",
          "Sustained infection"="#E57E7EFF", 
          "Transient infection"="#E57E7EFF", 
          "Abortive infection"="#E57E7EFF",
         'COVID_19'="#E57E7EFF")
merged_data_l2 <- merged_data_2 %>%
  group_by(covid_status,Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(covid_status,Condition, ct_merge.l2) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition) 
head(merged_data_l2)

merged_data_filter_l2 <- merged_data_l2%>% filter(ct_merge.l2 %in% c('CD4_Treg'))
#ct_merge.l2
pdf("figures/4_Line_CD4_Treg_filter_free.pdf", 3, 12)
 ggplot(merged_data_filter_l2, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(covid_status))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl3) +
  facet_wrap(covid_status~ ct_merge.l2, scales = "free", ncol = 1) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()
pdf("figures/4_Line_CD4_Treg_filter_free_2.pdf", 21, 2)
 ggplot(merged_data_filter_l2, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(covid_status))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl3) +
  facet_wrap(covid_status~ ct_merge.l2, scales = "free", ncol = 7) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()

merged_data_filter_l2 <- merged_data_l2%>% filter(ct_merge.l2 %in% c('Plasma'))
#ct_merge.l2
pdf("figures/4_Line_Plasma_filter_free.pdf", 3, 12)
 ggplot(merged_data_filter_l2, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(covid_status))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl3) +
  facet_wrap(covid_status~ ct_merge.l2, scales = "free", ncol = 1) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()
pdf("figures/4_Line_Plasma_filter_free_2.pdf", 21, 2)
 ggplot(merged_data_filter_l2, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(covid_status))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl3) +
  facet_wrap(covid_status~ ct_merge.l2, scales = "free", ncol = 7) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()

head(merged_data)

unique(merged_data$ct_merge)

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- merged_data %>%
    filter(covid_status == status) %>%
    select(-ifn_stim1, -HLA_score, -Count_condition, -Count_condition_celltype)
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_prop = mean(Celltype_prop, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>%
    group_by(ct_merge) %>%
    # 计算 p-value
    mutate(p_value = ifelse(
      n_distinct(Condition) > 1 & sum(Condition == baseline_condition) > 1 & sum(Condition != baseline_condition) > 1,
      
      # 使用 filter 筛选基准条件和其他条件的数据进行 t 检验
      t.test(filter(., Condition == baseline_condition)$Celltype_prop, 
             filter(., Condition != baseline_condition)$Celltype_prop)$p.value, NA )) %>%
    ungroup() %>%
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = Celltype_prop / Baseline_prop,
           FDR = p.adjust(p_value, method = "BH")) %>%
    # 移除不需要的列
    select(-Baseline_prop) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

final_merged_data$ct_merge <- factor(final_merged_data$ct_merge,
                               levels = rev(c('CD4T','CD8T','MAIT','gdT','Tdn','NK','Bcell','Monocyte','cDC','pDC','ILC','HPSC')))
unique(final_merged_data$ct_merge)
final_merged_data$covid_status <- factor(final_merged_data$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine','Sustained infection','Transient infection','Abortive infection'))
unique(final_merged_data$covid_status)


write.csv(final_merged_data,"figures/5_final_merged_data_FC.csv")

final_merged_data <- read.csv("figures/5_final_merged_data_FC.csv")
head(final_merged_data)

summary(final_merged_data$FoldChange_baseline)

summary(final_merged_data$FDR)

RdBu_pro <- c('#2166AC','#92C5DE','#F4A582','#D6604D','#B2182B')

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)


# 绘制热图
pdf("figures/5_Merge_ct_merge_Pro_Dotplot_with_background.pdf", width = 4.5, height = 1.68)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = FoldChange_baseline, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_pro) +
  scale_fill_manual(values = c("red", "blue", "#3A6963")) + # 设置矩形颜色
  guides(size = guide_legend(title = "Fold change over baseline"),
         color = guide_colorbar(title = "Fold change over baseline")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius( range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸 #limits = c(0, 0.5),
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

final_merged_data$ct_merge <- factor(final_merged_data$ct_merge,
                               levels = rev(c('CD4T','CD8T','MAIT','gdT','Tdn','NK','Bcell','Monocyte','cDC','pDC','ILC','HPSC')))
unique(final_merged_data$ct_merge)
final_merged_data$covid_status <- factor(final_merged_data$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine','Sustained infection','Transient infection','Abortive infection'))
unique(final_merged_data$covid_status)

# 创建热图数据
# 确保 FoldChange_baseline 是连续的
final_merged_data$FoldChange_baseline <- as.numeric(final_merged_data$FoldChange_baseline) 
Condition_seq <- c('A0','A1','A2','B0','B1','B2','C0','C1','C2','C3','D-1','D3','D7','D10','D14','D28')
final_merged_data$Condition	<- factor(final_merged_data$Condition,level=Condition_seq )
plot_merged_data <- final_merged_data 
unique(plot_merged_data$Condition)

rev(brewer.pal(11, "RdBu"))

library(ggplot2)
library(RColorBrewer)

# 优化后的颜色渐变调色板（减少节点复杂度）
RdBu_pro <- c("#053061", "#2166AC", "#F7F7F7", "#B2182B", "#67001F")

# 绘制热图
pdf("figures/5_Merge_ct_merge_Pro_Heatmap_optimized_Heatmap.pdf", width = 4.5, height = 1.68,
    useDingbats = TRUE, compress = TRUE)

heatmap_plot <- ggplot(plot_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_tile(aes(fill = FoldChange_baseline)) +
  scale_fill_gradientn(
    colors = RdBu_pro, 
    name = "Fold change over baseline",
    limits = c(0, 3.2),
    oob = scales::squish
  ) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 3, hjust = 0.5, face = "bold"),
    axis.line = element_blank(),  # 删除 X 和 Y 轴的线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5),
    axis.text.y = element_text(size = 4),
    axis.ticks.y = element_blank(),  # 删除 y 轴刻度线
    axis.ticks.x = element_blank(),  # 删除 x 轴刻度线
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"),
    legend.key.size = unit(0.05, "inch"),  # 减小 legend 的关键元素大小
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~covid_status, ncol = 5, scales = "free_x")

print(heatmap_plot)
dev.off()

library(ggplot2)
library(RColorBrewer)


# 选择调色板，这里使用 RdBu 调色板的颜色
#RdBu_pro <- rev(brewer.pal(11, "RdBu"))
RdBu_pro <- c('#053061','#2166AC','#F7F7F7','#FDDBC7','#F4A582','#D6604D','#B2182B','#67001F')

# 绘制热图
pdf("figures/5_Merge_ct_merge_Pro_Heatmap_with_A0.pdf", width = 4.5, height = 1.68)

heatmap_plot <- ggplot(plot_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_tile(aes(fill = FoldChange_baseline)) + # 使用 geom_tile 绘制热图
  scale_fill_gradientn(colors = RdBu_pro, name = "Fold change over baseline",
                       limits = c(0, 3.2),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放
dev.off()
# 提取 legend
legend <- cowplot::get_legend(heatmap_plot)

# 保存 legend 到单独的 PDF 文件
pdf("legend_only.pdf", width = 2, height = 1.5)
grid::grid.draw(legend)
dev.off()
print(heatmap_plot)

# 创建热图数据
# 确保 FoldChange_baseline 是连续的
final_merged_data$FoldChange_baseline <- as.numeric(final_merged_data$FoldChange_baseline) 
Condition_seq <- c('A0','A1','A2','B0','B1','B2','C0','C1','C2','C3','D-1','D3','D7','D10','D14','D28')
final_merged_data$Condition	<- factor(final_merged_data$Condition,level=Condition_seq )
plot_merged_data <- final_merged_data %>% filter(!Condition %in% c('A0','D-1'))
unique(plot_merged_data$Condition)

library(ggplot2)
library(RColorBrewer)


# 选择调色板，这里使用 RdBu 调色板的颜色
#RdBu_pro <- rev(brewer.pal(11, "RdBu"))
RdBu_pro <- c('#053061','#2166AC','#F7F7F7','#FDDBC7','#F4A582','#D6604D','#B2182B','#67001F')

# 绘制热图
pdf("figures/5_Merge_ct_merge_Pro_Heatmap.pdf", width = 4.5, height = 1.68)

heatmap_plot <- ggplot(plot_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_tile(aes(fill = FoldChange_baseline)) + # 使用 geom_tile 绘制热图
  scale_fill_gradientn(colors = RdBu_pro, name = "Fold change over baseline") + # 设置连续的颜色渐变
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

print(heatmap_plot)
dev.off()

data_for_pie <- read.csv("data/1_Merged_data_CVV_VIA_Challenge_with_Other.csv")
head(data_for_pie)

unique(data_for_pie$ct_merge)

library(ggplot2)
library(dplyr)
library(tidyr)

# 假设你的数据框名为 final_merged_data
# 计算每个 covid_status 的细胞类型占比
data_for_pie <- data_for_pie %>%
  group_by(covid_status, ct_merge) %>%
  summarise(count = n()) %>%
  group_by(covid_status) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
data_for_pie$ct_merge <- factor(data_for_pie$ct_merge,
                               levels = c('CD4T','CD8T','MAIT','gdT','Tdn','NK','Bcell','Monocyte',
                                          'cDC','pDC','ILC','HPSC','Platelet'))
unique(data_for_pie$ct_merge)
head(data_for_pie)

data_for_pie$covid_status <- factor(data_for_pie$covid_status,
                                level=c('Sustained infection','Transient infection',
                                        'Abortive infection','Inactivated Vaccine','mRNA Vaccine'))

library(ggplot2)

# 绘制饼图
pdf("figures/6_Merge_ct_merge_Pie.pdf", 4, 2)

ggplot(data_for_pie, aes(x = "", y = percentage, fill = ct_merge)) +
  geom_bar(width = 1, stat = "identity", alpha = 0.95, color = "white") +  # 设置透明度为0.75
  coord_polar(theta = "y") +
  facet_wrap(~ covid_status, ncol = 5) +  # 设置排列为5列
  scale_fill_manual(values = colorl2) +
  theme_void() +
  labs(fill = "Cell Type", title = "Cell Type Proportion by COVID Status") +
  theme(
    #legend.direction = "horizontal",  # 图例排列方式为水平
    #legend.box = "horizontal",  # 图例框为水平
    legend.key.size = unit(0.2, "cm"),  # 图例键的大小
    legend.text = element_text(size = 3.5),  # 图例文本大小
    legend.title = element_text(size = 3.5),  # 图例标题大小
    legend.spacing.x = unit(0.1, "cm"),  # 图例项之间的水平间距
    legend.position = 'right',#"bottom",  # 图例在底部
    legend.margin = margin(t = 0.3, b = 0.3, unit = "cm"),  # 图例的边距
    strip.text = element_text(size = 4),  # 每个饼图的标题字体大小
      plot.margin = unit(c(1, 1, 1, 1), "char"),
    plot.title = element_text(size = 5, hjust = 0.5, vjust = 10, face = "bold")  # 总图标题的字体大小，并调整位置
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))  # 设置图例为2行排列

dev.off()

library(ggplot2)

# 绘制饼图
pdf("figures/6_Merge_ct_merge_Pie_big.pdf", 15, 15)

ggplot(data_for_pie, aes(x = "", y = percentage, fill = ct_merge)) +
  geom_bar(width = 1, stat = "identity", alpha = 1, color = "white") +  # 设置透明度为0.75
  coord_polar(theta = "y") +
  facet_wrap(~ covid_status, ncol = 3) +  # 设置排列为5列
  scale_fill_manual(values = colorl2) +
  theme_void() +
  labs(fill = "Cell Type", title = "Cell Type Proportion by COVID Status") +
  theme(
    #legend.direction = "horizontal",  # 图例排列方式为水平
    #legend.box = "horizontal",  # 图例框为水平
    legend.key.size = unit(1, "cm"),  # 图例键的大小
    legend.text = element_text(size = 16),  # 图例文本大小
    legend.title = element_text(size = 18,hjust = 0.1,  face = "bold"),  # 图例标题大小
    legend.spacing.x = unit(1, "cm"),  # 图例项之间的水平间距
    legend.position = 'right',#"bottom",  # 图例在底部
    legend.margin = margin(t = 0.3, b = 0.3,l=0.6, unit = "cm"),  # 图例的边距
    strip.text = element_text(size = 26,  face = "bold"),  # 每个饼图的标题字体大小
      plot.margin = unit(c(1, 1, 1, 1), "char"),
    plot.title = element_text(size = 28, hjust = 0.5, vjust = -126, face = "bold")  # 总图标题的字体大小，并调整位置
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))  # 设置图例为2行排列

dev.off()

merged_data_2 <- rbind(mRNA_split_data, VIA_data_merge,cvv_data_merge, challenge_data_merge,cv_19_split_data)
unique(merged_data_2$covid_status)

merged_data_2 <- merged_data_2 %>%
  filter(! ct_merge %in% c('Platelet')) #三个数据统一有的细胞类型
merged_data_2$ct_merge <- factor(merged_data$ct_merge,
                               levels = rev(c('CD4T','CD8T','MAIT','gdT','Tdn',
                                              'NK','Bcell','Monocyte','cDC','pDC',
                                              'ILC','HPSC')))
unique(merged_data_2$ct_merge)
merged_data_2$covid_status <- factor(merged_data_2$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine','mRNA Vaccine(VIA)',
                                        'Sustained infection','Transient infection',
                                        'Abortive infection','COVID_19'))
unique(merged_data_2$covid_status)

unique(merged_data_2$Condition)

# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data_2 %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(VIA1, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(VIA1 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                          "mRNA Vaccine(VIA)" = "Day0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1",
                          "COVID_19" = "convalescent")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- average_scores %>%
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_score = mean(mean_score, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>% 
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = mean_score / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c('Inactivated Vaccine','mRNA Vaccine','mRNA Vaccine(VIA)',
                                        'Sustained infection','Transient infection',
                                        'Abortive infection','COVID_19'),
                        levels = c('Inactivated Vaccine','mRNA Vaccine','mRNA Vaccine(VIA)',
                                        'Sustained infection','Transient infection',
                                        'Abortive infection','COVID_19')), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue","blue", "#AB3282", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
final_merged_data$covid_status <- factor(final_merged_data$covid_status, 
                                      levels = c('Inactivated Vaccine','mRNA Vaccine',
                                                 'mRNA Vaccine(VIA)',
                                        'Sustained infection','Transient infection',
                                        'Abortive infection','COVID_19'))

# 绘制热图
pdf("figures/7_Merge_All_VIA_Dotplot_with_background_FoldChange_baseline.pdf", width = 7.2, height = 2.2)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new, name = "Fold change over baseline",
                       limits = c(0, 2),#limits = c(0.5, 1.5),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  scale_fill_manual(values = c("red", "blue", "#3A6963"))+
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression(Fold change over baseline)")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 7, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()



# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(VIA1, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(VIA1 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- average_scores %>%
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_score = mean(mean_score, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>% 
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = mean_score / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
final_merged_data$covid_status <- factor(final_merged_data$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/7_Merge_ct_merge_VIA_Dotplot_with_background_FoldChange_baseline.pdf", width = 5.5, height = 1.72)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new, name = "Fold change over baseline",
                       limits = c(0, 2),#limits = c(0.5, 1.5),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  scale_fill_manual(values = c("red", "blue", "#3A6963"))+
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression(Fold change over baseline)")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

merged_data <- read.csv("data/1_Merged_data_CVV_VIA_Challenge.csv")

merged_data <- merged_data %>%
  filter(! ct_merge %in% c('Platelet')) #三个数据统一有的细胞类型
merged_data$ct_merge <- factor(merged_data$ct_merge,
                               levels = rev(c('CD4T','CD8T','MAIT','gdT','Tdn',
                                              'NK','Bcell','Monocyte','cDC','pDC',
                                              'ILC','HPSC')))
unique(merged_data$ct_merge)
merged_data$covid_status <- factor(merged_data$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine','Sustained infection','Transient infection','Abortive infection'))
unique(merged_data$covid_status)
Condition_seq <- c('Day0','Day2','Day11','Day28','A0','A1','A2','B0','B1','B2','C0','C1','C2','D-1','D3','D7','D10','D14','D28')
merged_data$Condition	<- factor(merged_data$Condition,level=Condition_seq)
#head(merged_data)

# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(Apoptosis1, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(Apoptosis1 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- average_scores %>%
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_score = mean(mean_score, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>% 
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = mean_score / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
final_merged_data$covid_status <- factor(final_merged_data$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/8_Merge_ct_merge_Apoptosis_Dotplot_with_background_FoldChange_baseline.pdf", width = 5, height = 1.68)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new, name = "Fold change over baseline",
                       limits = c(0, 2),#limits = c(0.5, 1.5),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  scale_fill_manual(values = c("red", "blue", "#3A6963"))+
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression(Fold change over baseline)")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(ATG71, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(ATG71 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- average_scores %>%
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_score = mean(mean_score, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>% 
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = mean_score / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
final_merged_data$covid_status <- factor(final_merged_data$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/9_Merge_ct_merge_ATG7_Dotplot_with_background_FoldChange_baseline.pdf", width = 5.5, height = 1.72)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new, name = "Fold change over baseline",
                       limits = c(0, 2),#limits = c(0.5, 1.5),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  scale_fill_manual(values = c("red", "blue", "#3A6963"))+
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression(Fold change over baseline)")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

# 对每个 covid_status 内部的数据进行归一化
average_scores <- merged_data %>%
  group_by(covid_status, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(mac_ATG71, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(mac_ATG71 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

library(dplyr)
# 定义基准条件
baseline_conditions <- c("Inactivated Vaccine" = "A0",
                         "mRNA Vaccine" = "A0",
                         "Sustained infection" = "D-1",
                         "Transient infection" = "D-1",
                         "Abortive infection" = "D-1")
# 对所有 covid_status 进行循环计算
final_data <- lapply(names(baseline_conditions), function(status) {
  # 筛选出对应的 covid_status
  status_data <- average_scores %>%
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_score = mean(mean_score, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>% 
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = mean_score / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

# 创建 rect.data 数据框来定义矩形的 ymin 和 ymax
rect.data <- data.frame(
  covid_status = factor(c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"),
                        levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection")), # 定义每个状态的颜色顺序
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue", "#AB3282", "#AB3282", "#AB3282") # 定义每个状态的颜色
)

# 将 average_scores 的 covid_status 设置为因子，确保绘图时顺序正确
final_merged_data$covid_status <- factor(final_merged_data$covid_status, 
                                      levels = c("Inactivated Vaccine", "mRNA Vaccine", "Sustained infection", "Transient infection", "Abortive infection"))

# 绘制热图
pdf("figures/10_Merge_ct_merge_macroautophagy_ATG7_Dotplot_with_background_FoldChange_baseline.pdf", width = 5.5, height = 1.72)

ggplot(final_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE) + # 添加背景矩形
  geom_point(aes(size = percentage_expressing, color = FoldChange_baseline)) +
  #scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  scale_color_gradientn(colors = RdBu_new, name = "Fold change over baseline",
                       limits = c(0, 2),#limits = c(0.5, 1.5),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  scale_fill_manual(values = c("red", "blue", "#3A6963"))+
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression(Fold change over baseline)")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 4, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
    axis.line = element_line(colour = "black", size = 0.5), # 保留 X 轴和 Y 轴的线
    panel.grid = element_blank(), # 去除背景格线
    axis.text.x = element_text(size = 4, angle = 270, hjust = 0.1, vjust = 0.5), # 设置 X 轴文本
    axis.text.y = element_text(size = 4), # 设置 Y 轴文本
    axis.ticks.y = element_line(size = 0.3), # 调整 y 轴刻度线的宽度
    axis.ticks.x = element_line(size = 0.3), # 调整 x 轴刻度线的宽度
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    legend.position = 'right',
    legend.spacing = unit(0.03, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.07, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 1.3)) + # 调整点的范围，缩小点的尺寸
  facet_wrap(~covid_status, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()


