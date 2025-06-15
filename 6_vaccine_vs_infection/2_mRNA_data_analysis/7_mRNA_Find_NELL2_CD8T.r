library(Seurat)
library("tidyverse")
library("doParallel")
library(BiocParallel)
 library(Seurat)
 library(future)


source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
source("/share/home/qlab/projects/project_cvv/yyp_results_3/yyp_function.R") #加载到当前的R环境中

getwd()

dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/4_new_celltype'

mRNA_obj<- readRDS(paste0(dir_path,"/1_mRNA_final_CCA_integrated.rds"))

unique(mRNA_obj$celltype.l1)
unique(mRNA_obj$group)
unique(mRNA_obj$celltype.l2)
colnames(mRNA_obj@meta.data)

Idents(mRNA_obj) <-  "group"
targetcell  <- subset(mRNA_obj,
             idents = c('Vax','HC','booster'))
Idents(targetcell) <-  "celltype.l1"
targetcell  <- subset(targetcell,
             idents = c('CD8 T'))
unique(targetcell$celltype.l2)
nrow(targetcell@meta.data)

results_path <- paste0(getwd(),'/NELL2_figures')
results_path 

# 加载必要的库
library(Seurat)
library(dplyr)

unique(targetcell$celltype.l2) #targetcell为全部的CD8T

# 加载参考数据集和目标数据集（假设数据已标准化）
# Reference 是参考数据集，包含细胞类型标注，目标数据集是 Query
query <- targetcell #已标注细胞
reference <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/4_Tcell_CD8T/2_Reumap/Trajectory_CD8T.rds")

# 检查参考数据集中是否有细胞类型标签
table(reference$ct_level_4)

common_genes <- intersect(rownames(reference), rownames(query))
length(common_genes)
# 确保两个数据集使用相同的基因集
reference <- subset(reference, features = common_genes)
query <- subset(query, features = common_genes)

# 标准化和 PCA 降维
reference <- NormalizeData(reference) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
query <- NormalizeData(query) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# 找到锚点（anchors）
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  dims = 1:30  # 使用 PCA 的前 30 个主成分
)
# 转移标签
query_res <- TransferData(
  anchorset = anchors,
  refdata = reference$ct_level_4,  # 使用参考数据集中的细胞类型标签
  dims = 1:30
)
# 转移结果存储在 query 的 metadata 中
query@meta.data$predicted.id <- query_res$predicted.id
colnames(query@meta.data)

colnames(query_res)
nrow(query_res)
nrow(query@meta.data)

DefaultAssay(query) <- "RNA"                         
idents_name_list <- c("celltype.l1",
                      "celltype.l2",
                      "ct_merge",
                     "predicted.id")
pdf(file.path(results_path, "3_CCA_predicted_DimPlot.pdf"),10,5)
for (idents_name in idents_name_list) {
  Idents(query) <- idents_name
  p <- DimPlot(object = query, reduction = 'pca',label = T,raster = FALSE) 
  print(p)
}
dev.off()

table(query@meta.data$predicted.id)

table(query@meta.data$predicted.id,query@meta.data$celltype.l2)

# 验证标注结果
table(query@meta.data$predicted.id,query@meta.data$timepoint)  # 检查分布

write.csv(query@meta.data,'data/3_mRNA_CD8T_NELL2_Predicted.csv')

query_meta.data <- read.csv('data/3_mRNA_CD8T_NELL2_Predicted.csv')
nrow(query_meta.data)
colnames(query_meta.data)

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
    'booster_d120' = 'C3'#,           
    # 'acute' = 'acute',
    #'cv19_d11' = 'acute',
    #'convalescent' = 'convalescent'
)

# 使用 dplyr::mutate 和 recode 应用映射
query_meta.data <- query_meta.data %>%
  mutate(Condition = recode(timepoint, !!!time_mapping))
# 查看映射结果
unique(query_meta.data[, c("timepoint", "Condition")])

colnames(query_meta.data)

table(query_meta.data$predicted.id)

table(query_meta.data$predicted.id,query_meta.data$celltype.l2)

colpalette <- c(CD8_NELL2 ="#53A85F",
           CD8_Tem="#99AF09",
           CD8_Temra="#68A180",
           CD8_Tn="#D8E7A9")

pdf("NELL2_figures/4_mRNA_CD8T_Line_Condition_predicted_id.pdf", height = 3, width = 20)
#predicted.id
meta_data_1.df <- query_meta.data %>%
  group_by(Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(Condition, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(predicted.id))) +
  geom_point(size=5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ predicted.id, scales = "free_y", ncol = 4) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 
dev.off()

colnames(meta_data_1.df)[c(23,105:109)]

head(meta_data_1.df[colnames(meta_data_1.df)[c(23,105:109)]],n=1)
CD8T_transfer <- meta_data_1.df[colnames(meta_data_1.df)[c(23,105:109)]]
write.csv(CD8T_transfer,'data/5_merge_mRNA_CD8T_transfer.csv')





colpalette_2 <- c('CD8/CD28 T EM' ="#53A85F",
           'CD8 T EM'="#99AF09",
           'CD8 T Eff'="#68A180",
           'Naive CD8 T'="#D8E7A9")
pdf("NELL2_figures/4_mRNA_CD8T_Line_Condition_celltype_l2.pdf", height = 3, width = 20)
#predicted.id
meta_data_1.df <- query_meta.data %>%
  group_by(Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(Condition, celltype.l2) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(celltype.l2))) +
  geom_point(size=5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette_2) +
  facet_wrap(~ celltype.l2, scales = "free_y", ncol = 4) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 
dev.off()

# 保存 PDF
pdf("NELL2_figures/4_mRNA_CD8T_predicted_pie.pdf", height = 6, width = 6)

# 按 ct_level_4 计算比例
meta_data_1.df <- query_meta.data %>%
  group_by(predicted.id) %>%
  summarise(Count_celltype = n()) %>%
  mutate(Celltype_prop = Count_celltype / sum(Count_celltype))

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(predicted.id))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colpalette) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Overall Cell Type Proportion") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()

colpalette_2 <- c('CD8/CD28 T EM' ="#53A85F",
           'CD8 T EM'="#99AF09",
           'CD8 T Eff'="#68A180",
           'Naive CD8 T'="#D8E7A9")

# 保存 PDF
pdf("NELL2_figures/4_mRNA_CD8T_celltype_l2_pie.pdf", height = 6, width = 6)

# 按 ct_level_4 计算比例
meta_data_1.df <- query_meta.data %>%
  group_by(celltype.l2) %>%
  summarise(Count_celltype = n()) %>%
  mutate(Celltype_prop = Count_celltype / sum(Count_celltype))

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(celltype.l2))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colpalette_2) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Overall Cell Type Proportion") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()





CVV_CD8T <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/4_Tcell_CD8T/2_Reumap/Trajectory_CD8T.rds")

pdf("NELL2_figures/4_CVV_CD8T_Line_Condition_ct_level_4.pdf",  height = 3, width = 20)
#ct_level_4
meta_data_1.df <- CVV_CD8T@meta.data %>%
  group_by(Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(Condition, ct_level_4) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(ct_level_4))) +
  geom_point(size=5) +
  geom_line(size=1.4) +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ ct_level_4, scales = "free_y", ncol = 4) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 
dev.off()

head(meta_data_1.df[c('Condition',
                          'ct_level_4',
                           'Count_condition',
                           'Count_condition_celltype',
                           'Celltype_prop')],n=1)

CD8T_transfer <- meta_data_1.df[c('Condition',
                          'ct_level_4',
                           'Count_condition',
                           'Count_condition_celltype',
                           'Celltype_prop')]
write.csv(CD8T_transfer,'data/5_merge_CVV_CD8T_transfer.csv')





# 保存 PDF
pdf("NELL2_figures/4_CVV_CD8T_Line_Condition_ct_level_4_pie.pdf", height = 6, width = 6)

# 按照 Condition 和 ct_level_4 计算比例
meta_data_1.df <- CVV_CD8T@meta.data %>%
  group_by(Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup() %>% 
  group_by(Condition, ct_level_4) %>%
  mutate(Count_condition_celltype = n(), Celltype_prop = Count_condition_celltype / Count_condition) %>%
  distinct(Condition, ct_level_4, .keep_all = TRUE)

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(ct_level_4))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ Condition, ncol = 3) +
  scale_fill_manual(values = colpalette) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Cell Type Proportion by Condition") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()

# 保存 PDF
pdf("NELL2_figures/4_CVV_CD8T_ct_level_4_pie.pdf", height = 6, width = 6)

# 按 ct_level_4 计算比例
meta_data_1.df <- CVV_CD8T@meta.data %>%
  group_by(ct_level_4) %>%
  summarise(Count_celltype = n()) %>%
  mutate(Celltype_prop = Count_celltype / sum(Count_celltype))

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(ct_level_4))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colpalette) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Overall Cell Type Proportion") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()


