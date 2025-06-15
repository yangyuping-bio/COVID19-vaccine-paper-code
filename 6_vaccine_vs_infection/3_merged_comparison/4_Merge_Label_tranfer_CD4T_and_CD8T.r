library(dplyr)
library(ggplot2)
library(ggplot2)
library(RColorBrewer)

mRNA_CD8T <- read.csv('data/5_merge_mRNA_CD8T_transfer.csv')
mRNA_CD4T <- read.csv('data/5_merge_mRNA_CD4T_transfer.csv')
head(mRNA_CD8T,n=1)
head(mRNA_CD4T,n=1)
Challenge_CD8T <- read.csv('data/5_merge_Challenge_CD8T_transfer.csv')
Challenge_CD4T <- read.csv('data/5_merge_Challenge_CD4T_transfer.csv')
head(Challenge_CD8T,n=1)
head(Challenge_CD4T,n=1)
CVV_CD8T <- read.csv('data/5_merge_CVV_CD8T_transfer.csv')
CVV_CD4T <- read.csv('data/5_merge_CVV_CD4T_transfer.csv')
head(CVV_CD8T,n=1)
head(CVV_CD4T,n=1)

mRNA_CD8T$covid_status <- 'mRNA Vaccine'
mRNA_CD8T$cell_type <- mRNA_CD8T$celltype.l2
mRNA_CD4T$covid_status <- 'mRNA Vaccine'
mRNA_CD4T$cell_type <- mRNA_CD4T$celltype.l2

CVV_CD8T$covid_status <- 'Inactivated Vaccine'
CVV_CD4T$covid_status <- 'Inactivated Vaccine'
CVV_CD8T$cell_type <- CVV_CD8T$ct_level_4
CVV_CD4T$cell_type <- CVV_CD4T$ct_level_4
CVV_CD8T$predicted.id <- CVV_CD8T$ct_level_4
CVV_CD4T$predicted.id <- CVV_CD4T$ct_level_4

Challenge_CD8T$Condition <- Challenge_CD8T$time_point
Challenge_CD4T$Condition <- Challenge_CD4T$time_point

colnames(mRNA_CD8T)
colnames(CVV_CD8T)
colnames(Challenge_CD8T)

merge_seq <- c('predicted.id','Condition','Count_condition',
               'Count_condition_celltype','Celltype_prop','covid_status','cell_type')
mRNA_CD8T_merge <- mRNA_CD8T[,merge_seq ]
cvv_CD8T_merge <- CVV_CD8T[,merge_seq ]
challenge_CD8T_merge <- Challenge_CD8T[,merge_seq ]
merged_CD8T <- rbind(mRNA_CD8T_merge, cvv_CD8T_merge, challenge_CD8T_merge)
merged_CD8T$ct_level1 <- 'CD8T'
head(merged_CD8T,n=2)

merge_seq <- c('predicted.id','Condition','Count_condition',
               'Count_condition_celltype','Celltype_prop','covid_status','cell_type')
mRNA_CD4T_merge <- mRNA_CD4T[,merge_seq ]
cvv_CD4T_merge <- CVV_CD4T[,merge_seq ]
challenge_CD4T_merge <- Challenge_CD4T[,merge_seq ]
merged_CD4T <- rbind(mRNA_CD4T_merge, cvv_CD4T_merge, challenge_CD4T_merge)
merged_CD4T$ct_level1 <- 'CD4T'
head(merged_CD4T,n=2)

merged_data <- rbind(merged_CD8T,merged_CD4T)
head(merged_data,n=2)

merged_data$covid_status <- factor(merged_data$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine','Sustained infection','Transient infection','Abortive infection'))

write.csv(merged_data,'data/1_merged_data_label_Transfer_CD4T_and_CD8T.csv')

merged_data <- read.csv('data/1_merged_data_label_Transfer_CD4T_and_CD8T.csv')

colorl2 <- c(
  CD8_NELL2 ="#23452F",
  CD8_Tem="#99AF09",
  CD8_Temra="#68A180",
  CD8_Tn="#D8E7A9",
 CD4_T_CCR6="#476D87", 
 CD4_Tcm="#C1E6F3", 
 CD4_Tem="#91D0BE", 
 CD4_Temra="#58A4C3", 
 CD4_Th2="#57C3F3", 
 CD4_Tn="#6778AE", 
 CD4_Treg="#2e409a"
)

merged_data$predicted.id <- factor(merged_data$predicted.id,
                                level=c('CD8_Tn','CD8_NELL2','CD8_Tem','CD8_Temra',
           'CD4_Tn','CD4_Tcm','CD4_T_CCR6','CD4_Th2','CD4_Treg','CD4_Tem', 'CD4_Temra'))

#data_for_pie <- read.csv("data/1_Merged_data_CVV_VIA_Challenge_with_Other.csv")
#head(data_for_pie)

library(ggplot2)
library(dplyr)
library(tidyr)

# 计算每个 covid_status 的细胞类型占比
data_for_pie <- merged_data %>%
  group_by(ct_level1,covid_status, predicted.id) %>%
  summarise(count = n()) %>%
  group_by(ct_level1,covid_status) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

head(data_for_pie)

data_for_pie$covid_status <- factor(data_for_pie$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine',
                                       'Sustained infection','Transient infection',
                                        'Abortive infection'))
data_for_pie$ct_level1 <- factor(data_for_pie$ct_level1,
                                level=c('CD8T','CD4T'))
data_for_pie$predicted.id <- factor(data_for_pie$predicted.id,
                                level=c('CD8_Tn','CD8_NELL2','CD8_Tem','CD8_Temra',
           'CD4_Tn','CD4_Tcm','CD4_T_CCR6','CD4_Th2','CD4_Treg','CD4_Tem', 'CD4_Temra'))

library(ggplot2)

# 绘制饼图
pdf("figures_CD4T_CD8T/2_Merge_predicted_Pie_big.pdf", 22, 15)

ggplot(data_for_pie, aes(x = "", y = percentage, fill = predicted.id)) +
  geom_bar(width = 1, stat = "identity", alpha = 1, color = "white") +  # 设置透明度为0.75
  coord_polar(theta = "y") +
  facet_wrap(ct_level1~ covid_status, ncol = 5) +  # 设置排列为5列
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

tail(merged_data)

#检查数据
merged_data <- merged_data %>%
  group_by(ct_level1,covid_status,Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(ct_level1,covid_status,Condition, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition) 
tail(merged_data)

unique(merged_data$covid_status)

summary(merged_data$Count_condition_celltype)

pdf("figures_CD4T_CD8T/1_Line_Condition_predicted.id_filter_free.pdf", 27, 12)
 ggplot(merged_data, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(predicted.id))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl2) +
  facet_wrap(covid_status~ predicted.id, scales = "free") +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()

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
    filter(covid_status == status) 
  # 获取对应的基准条件
  baseline_condition <- baseline_conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_condition) %>%
    group_by(predicted.id) %>%
    summarise(Baseline_score = mean(Celltype_prop, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "predicted.id") %>% 
    # 计算 FoldChange 
    mutate(FoldChange_baseline = Celltype_prop / Baseline_score) %>%
    # 移除不需要的列
    select(-Baseline_score) %>%
    # 保留 covid_status 以便最后合并
    mutate(covid_status = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data,n=1)
#write.csv(final_merged_data,"data/1_final_merged_data_FC.csv")
#final_merged_data <- read.csv("data/1_final_merged_data_FC.csv")

# 创建热图数据
# 确保 FoldChange_baseline 是连续的
final_merged_data$covid_status <- factor(final_merged_data$covid_status,
                                level=c('Inactivated Vaccine','mRNA Vaccine','Sustained infection','Transient infection','Abortive infection'))
final_merged_data$predicted.id <- factor(final_merged_data$predicted.id,
                                         level=rev(unique(final_merged_data$predicted.id)))
final_merged_data$FoldChange_baseline <- as.numeric(final_merged_data$FoldChange_baseline) 
Condition_seq <- c('A0','A1','A2','B0','B1','B2','C0','C1','C2','C3','D-1','D3','D7','D10','D14','D28')
final_merged_data$Condition	<- factor(final_merged_data$Condition,level=Condition_seq )
plot_merged_data <- final_merged_data 
unique(plot_merged_data$Condition)

final_merged_data$predicted.id <- factor(final_merged_data$predicted.id,
                                level=rev(c('CD8_Tn','CD8_NELL2','CD8_Tem','CD8_Temra',
           'CD4_Tn','CD4_Tcm','CD4_T_CCR6','CD4_Th2','CD4_Treg','CD4_Tem', 'CD4_Temra')))

 rev(brewer.pal(11, "RdBu"))

#RdBu_pro <- rev(brewer.pal(11, "RdBu"))
color5 <- c('#053061','#2166AC','#D1E5F0','#FDDBC7','#F4A582','#D6604D','#B2182B')
RdBu_pro <- colorRampPalette(color5)(20)  

# 绘制热图
pdf("figures_CD4T_CD8T/1_Merge_predicted_Heatmap.pdf", width = 4.2, height = 1.5)

heatmap_plot <- ggplot(final_merged_data , aes(x = Condition, y = predicted.id)) +
  geom_tile(aes(fill = FoldChange_baseline)) + # 使用 geom_tile 绘制热图
  scale_fill_gradientn(colors = RdBu_pro, name = "Fold change over baseline",
                       limits = c(0, 2.33),
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + # 设置连续的颜色渐变
  theme(
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    strip.text = element_text(size = 3, hjust = 0.5, face = "bold"), # 分图标题字体设置为4，居中
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






