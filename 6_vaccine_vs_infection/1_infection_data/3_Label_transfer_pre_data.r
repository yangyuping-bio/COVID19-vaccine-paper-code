# 加载必要的包
library(Seurat)
 #library(rhdf5)
#library(feather)
#library(Matrix)
library("tidyverse")
library(dplyr)

seurat_obj <- readRDS("data/1_challenge_seurat_obj.rds")

seurat_obj

unique(seurat_obj$cell_type)

# 设置为正确的 Assay
DefaultAssay(seurat_obj) <- "RNA"  # 或者其他包含数据的 Assay
Idents(seurat_obj) <-  "cell_type"
targetcell  <- subset(seurat_obj,
             idents = c('T CD4 Naive', 'T CD4 Helper', 'T Reg','T CD4 CTL'))
unique(targetcell$cell_type)
nrow(targetcell@meta.data)

query <- targetcell
unique(targetcell$cell_type)

data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3"
reference <-readRDS(file.path(data_save_path, "4_Tcell_CD4T/3_ReCCA_with_Treg/CD4T_ReCCA_with_Treg.rds"))
unique(reference$ct_level_4)

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

# 将 query@meta.data$cell_type 从 factor 转换为 character
query@meta.data$cell_type <- as.character(query@meta.data$cell_type)

# 检查数据类型
unique(query@meta.data$cell_type)

table(query@meta.data$cell_type)
table(query@meta.data$predicted.id,query@meta.data$cell_type)
write.csv(query@meta.data,'data/4_CD4T_Lable_transfer_predicted.csv') 

CD4T_res <- read.csv('data/4_CD4T_Lable_transfer_predicted.csv') 
table(CD4T_res$predicted.id,CD4T_res$time_point)
table(CD4T_res$predicted.id,CD4T_res$cell_type)
head(CD4T_res,n=2)

CD4T_res$time_point <- factor(CD4T_res$time_point,
                              levels=c('D-1','D3','D7',
                                       'D10','D14','D28'))

meta_data_1.df <- CD4T_res %>%
  group_by(covid_status,time_point) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(covid_status,time_point, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

colnames(meta_data_1.df)[c(6,7,12,17:20)]

head(meta_data_1.df[colnames(meta_data_1.df)[c(6,7,12,17:20)]],n=1)
CD4T_transfer <- meta_data_1.df[colnames(meta_data_1.df)[c(6,7,12,17:20)]]
write.csv(CD4T_transfer,'data/5_merge_Challenge_CD4T_transfer.csv')



# 设置为正确的 Assay
DefaultAssay(seurat_obj) <- "RNA"  # 或者其他包含数据的 Assay
Idents(seurat_obj) <-  "cell_type"
targetcell  <- subset(seurat_obj,
             idents = c('T CD8 Naive', 'T CD8 Memory', 'T CD8 CTL'))
unique(targetcell$cell_type)
nrow(targetcell@meta.data)

query <- targetcell

unique(targetcell$cell_type)

reference <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/4_Tcell_CD8T/2_Reumap/Trajectory_CD8T.rds")

unique(reference$ct_level_4)

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

# 将 query@meta.data$cell_type 从 factor 转换为 character
query@meta.data$cell_type <- as.character(query@meta.data$cell_type)

# 检查数据类型
unique(query@meta.data$cell_type)

table(query@meta.data$cell_type)
table(query@meta.data$predicted.id,query@meta.data$cell_type)
write.csv(query@meta.data,'data/4_CD8T_Lable_transfer_predicted.csv') 

# 设置为正确的 Assay
DefaultAssay(seurat_obj) <- "RNA"  # 或者其他包含数据的 Assay
Idents(seurat_obj) <-  "cell_type"
targetcell  <- subset(seurat_obj,
             idents = c('Monocyte CD14+', 'Monocyte CD16+'))
unique(targetcell$cell_type)
nrow(targetcell@meta.data)

query <- targetcell

data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3"
reference <- readRDS(file.path(data_save_path, "2_Mono/2_Reumap/Trajectory_Mono.rds"))

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


# 将 query@meta.data$cell_type 从 factor 转换为 character
query@meta.data$cell_type <- as.character(query@meta.data$cell_type)

# 检查数据类型
unique(query@meta.data$cell_type)

table(query@meta.data$predicted.id)

table(query@meta.data$predicted.id,query@meta.data$cell_type)

write.csv(query@meta.data,'data/4_Mono_Lable_transfer_predicted.csv') 





CD8T_res <- read.csv('data/4_CD8T_Lable_transfer_predicted.csv') 
table(CD8T_res$predicted.id,CD8T_res$time_point)
table(CD8T_res$predicted.id,CD8T_res$cell_type)
head(CD8T_res,n=2)

colpalette <- c(CD8_NELL2 ="#53A85F",
           CD8_Tem="#99AF09",
           CD8_Temra="#68A180",
           CD8_Tn="#D8E7A9")
CD8T_res$time_point <- factor(CD8T_res$time_point,
                              levels=c('D-1','D3','D7',
                                       'D10','D14','D28'))

library("tidyverse")
pdf("figures/5_Challenge_CD8T_Line_Time_predicted_id.pdf", height = 3, width = 20)
#predicted.id
meta_data_1.df <- CD8T_res %>%
  group_by(time_point) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(time_point, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = time_point, y = Celltype_prop, group = 1, color=factor(predicted.id))) +
  geom_point(size=5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ predicted.id, scales = "free_y", ncol = 3) +
  labs(y = "Proportion", x = "time_point") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 

dev.off()

library("tidyverse")
pdf("figures/5_Challenge_CD8T_Line_Time_status_predicted_id.pdf", height = 9, width = 17)
#predicted.id
meta_data_1.df <- CD8T_res %>%
  group_by(covid_status,time_point) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(covid_status,time_point, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = time_point, y = Celltype_prop, group = 1, color=factor(predicted.id))) +
  geom_point(size=5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette) +
  facet_wrap(covid_status ~ predicted.id, scales = "free_y", ncol = 3) +
  labs(y = "Proportion", x = "time_point") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 
dev.off()

colnames(meta_data_1.df)[c(6,7,12,17:20)]

head(meta_data_1.df[colnames(meta_data_1.df)[c(6,7,12,17:20)]],n=1)
CD8T_transfer <- meta_data_1.df[colnames(meta_data_1.df)[c(6,7,12,17:20)]]
write.csv(CD8T_transfer,'data/5_merge_Challenge_CD8T_transfer.csv')

# 保存 PDF
pdf("figures/5_Challenge_CD8T_predicted_pie.pdf", height = 6, width = 6)

# 按 ct_level_4 计算比例
meta_data_1.df <- CD8T_res %>%
  group_by(predicted.id) %>%
  summarise(Count_celltype = n()) %>%
  mutate(Celltype_prop = Count_celltype / sum(Count_celltype))

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(predicted.id))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colpalette) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Predicted Proportion") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()

unique(CD8T_res$cell_type)

colpalette_2 <- c(
           'T CD8 Memory'="#99AF09",
           'T CD8 CTL'="#68A180",
           'T CD8 Naive'="#D8E7A9")

# 保存 PDF
pdf("figures/5_Challenge_CD8T_cell_type_pie.pdf", height = 6, width = 6)

meta_data_1.df <- CD8T_res %>%
  group_by(cell_type) %>%
  summarise(Count_celltype = n()) %>%
  mutate(Celltype_prop = Count_celltype / sum(Count_celltype))

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(cell_type))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colpalette_2) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Cell Type Proportion") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()

Mono_res <- read.csv('data/4_Mono_Lable_transfer_predicted.csv') 
table(Mono_res$predicted.id,Mono_res$time_point)
table(Mono_res$predicted.id,Mono_res$cell_type)
head(Mono_res,n=2)

#Mono
colpalette <- c(Mono_CD14 ="#E95C59",
  'Mono_CD14_ATG7+' ="#F1BB72", 
  Mono_CD14_CD16 ="#F3B1A0",
  Mono_CD16 ="#E5D2DD", 
  'Mono_CD16_ATG7+' ="#476D87")

unique(Mono_res$time_point)
Mono_res$time_point <- factor(Mono_res$time_point,
                              levels=c('D-1','D3','D7',
                                       'D10','D14','D28'))

library("tidyverse")
pdf("figures/5_Challenge_Mono_Line_Time_predicted_id.pdf", height = 3, width = 25)
#predicted.id
meta_data_1.df <- Mono_res %>%
  group_by(time_point) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(time_point, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = time_point, y = Celltype_prop, group = 1, color=factor(predicted.id))) +
  geom_point(size=5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ predicted.id, scales = "free_y", ncol = 5) +
  labs(y = "Proportion", x = "time_point") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 
dev.off()

unique(Mono_res$covid_status)

table(Mono_res$covid_status,Mono_res$time_point)

library("tidyverse")
pdf("figures/5_Challenge_Mono_Line_Time_status_predicted_id.pdf", height = 9, width = 25)
#predicted.id
meta_data_1.df <- Mono_res %>%
  group_by(covid_status,time_point) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(covid_status,time_point, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = time_point, y = Celltype_prop, group = 1, color=factor(predicted.id))) +
  geom_point(size=5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette) +
  facet_wrap(covid_status ~ predicted.id, scales = "free_y", ncol = 5) +
  labs(y = "Proportion", x = "time_point") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 
dev.off()

colnames(meta_data_1.df)

head(meta_data_1.df[colnames(meta_data_1.df)[c(6,7,12,17:20)]],n=1)
Mono_transfer <- meta_data_1.df[colnames(meta_data_1.df)[c(6,7,12,17:20)]]
write.csv(Mono_transfer,'data/5_merge_Challenge_Mono_transfer.csv')

library("tidyverse")
pdf("figures/5_Challenge_Mono_Line_Time_predicted_id.pdf", height = 3, width = 25)
#predicted.id
meta_data_1.df <- Mono_res %>%
  group_by(time_point) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(time_point, predicted.id) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition)

 ggplot(meta_data_1.df, aes(x = time_point, y = Celltype_prop, group = 1, color=factor(predicted.id))) +
  geom_point(size=5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ predicted.id, scales = "free_y", ncol = 5) +
  labs(y = "Proportion", x = "time_point") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #legend.position = "none"
       ) 
dev.off()

# 保存 PDF
pdf("figures/5_Challenge_Mono_predicted_pie.pdf", height = 6, width = 6)

# 按 ct_level_4 计算比例
meta_data_1.df <- Mono_res %>%
  group_by(predicted.id) %>%
  summarise(Count_celltype = n()) %>%
  mutate(Celltype_prop = Count_celltype / sum(Count_celltype))

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(predicted.id))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colpalette) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Predicted Proportion") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()

unique(Mono_res$cell_type)

colpalette_2 <- c('Monocyte CD14+' ="#E95C59",
 'Monocyte CD16+' ="#E5D2DD")

# 保存 PDF
pdf("figures/5_Challenge_Mono_cell_type_pie.pdf", height = 6, width = 6)

meta_data_1.df <- Mono_res %>%
  group_by(cell_type) %>%
  summarise(Count_celltype = n()) %>%
  mutate(Celltype_prop = Count_celltype / sum(Count_celltype))

# 绘制饼图
ggplot(meta_data_1.df, aes(x = "", y = Celltype_prop, fill = factor(cell_type))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = colpalette_2) +
  labs(fill = "Cell Type", y = NULL, x = NULL, title = "Cell Type Proportion") +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))

# 关闭 PDF
dev.off()




