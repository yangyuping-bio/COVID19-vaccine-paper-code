library(Seurat)
library(ComplexHeatmap)
library(cols4all)
library("GSEABase")
library("tidyverse")
library(org.Hs.eg.db)
library(tidyverse)
library(org.Hs.eg.db) ##加载人类
 require(clusterProfiler)
 require(clusterProfiler)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

source("/share/home/qlab/projects/project_cvv/yyp_results_3/yyp_function.R") #加载到当前的R环境中

getwd()

mRNA_obj <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/data/GSE247917_cv19_combined_obj.rds")

nrow(mRNA_obj@meta.data)

str(mRNA_obj@meta.data)

colnames(mRNA_obj@meta.data)

head(mRNA_obj@meta.data)

unique(mRNA_obj@meta.data$timepoint)
unique(mRNA_obj@meta.data$timepoint)

table(mRNA_obj@meta.data$celltype.l1)

table(mRNA_obj@meta.data$celltype.l2)

table(mRNA_obj@meta.data$celltype.l1,mRNA_obj@meta.data$timepoint)

summary(mRNA_obj@meta.data$celltype.l2)

getwd()

mRNA_obj@meta.data$timepoint <- factor(
  mRNA_obj@meta.data$timepoint,
  levels = c('HC','Vax_d7','Vax_d10','Vax_d14','Vax_d20','Vax_d21',
             'Vax_d28','Vax_d29','Vax_d35','Vax_d36',
             'booster_d0','booster_d7',
             'booster_d28','booster_d120',             
             'acute','cv19_d11','convalescent'))

library("tidyverse")
library("ggpubr")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(RColorBrewer)
display.brewer.all()

color_l3 <- c(
  "#E5D2DD", "#F3B1A0", "#E59CC4", "#AB3282", # B细胞类型
  "#CCE0F5", "#C1E6F3", "#91D0BE", "#58A4C3", "#57C3F3", "#6778AE", "#476D87",# CD4 T细胞类型
  "#D6E7A9", "#C5DEBA", "#68A180", "#53A85F", # CD8 T细胞类型
  "#E0D4CA", "#9FA3A8", # CDC细胞类型
  "#3A6963", "#23452F", #gdT
  "#585658",#HPSC
  "#F1CC92","#F1BC62", #ILC
  "#5F3D69",#MAIT
  "#E95C39", "#A01319", "#B53E2B", "#712820", "#A05401",# Mono细胞类型 
  "#E4C955", "#E1A111", "#AA9A59", # NK细胞类型
  "#968175","#8C549C", "#BD956A","#585658"# 其他细胞类型
)

colpalette <- color_l3  
pdf("figures/1_Line_timepoint_celltype.l1.pdf", height = 10, width = 16)
#celltype.l1
meta_data.df <-data.frame(mRNA_obj@meta.data)
meta_data.df <- meta_data.df %>%
  group_by(timepoint) %>%
  mutate(Count_timepoint = n()) %>%
  ungroup%>% 
  group_by(timepoint, celltype.l1) %>%
  mutate(Count_timepoint_celltype = n(),Celltype_prop = Count_timepoint_celltype/Count_timepoint) 

 ggplot(meta_data.df, aes(x = timepoint, y = Celltype_prop, group = 1, color=factor(celltype.l1))) +
  geom_point(size=3.5) +
  geom_line(size=1.5) +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ celltype.l1, scales = "free_y", ncol = 3) +
  labs(y = "Proportion", x = "timepoint") +
  theme_minimal(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.position = "none"
       ) 
dev.off()

meta_data.df <-data.frame(mRNA_obj@meta.data)
library(ggplot2)
library(dplyr)
colpalette <- color_l3

# 更新 meta_data.df 数据
meta_data.df <- meta_data.df %>%
  group_by(timepoint) %>%
  mutate(Count_timepoint = n()) %>%
  ungroup() %>% 
  group_by(timepoint, celltype.l2) %>%
  mutate(Count_timepoint_celltype = n(), Celltype_prop = Count_timepoint_celltype / Count_timepoint) 


meta_data.df$timepoint <- factor(
  meta_data.df$timepoint,
  levels = c('HC','Vax_d7','Vax_d10','Vax_d14','Vax_d20','Vax_d21',
             'Vax_d28','Vax_d29','Vax_d35','Vax_d36',
             'booster_d0','booster_d7',
             'booster_d28','booster_d120',             
             'acute','cv19_d11','convalescent'))

# 绘制图形
pdf("figures/1_Line_timepoint_celltype_l2.pdf", height = 18, width = 16)
#figures/
ggplot(meta_data.df, aes(x = timepoint, y = Celltype_prop, group = 1, color = factor(celltype.l2))) +
  geom_point(size = 3.5) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ celltype.l2, scales = "free_y", ncol = 3) +
  labs(y = "Proportion", x = "timepoint") +
  theme_minimal(base_size = 26) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

dev.off()









Idents(mRNA_obj) <- 'celltype.l1'
CD4T <- subset(mRNA_obj,idents=c('CD4 T'))
nrow(CD4T@meta.data)

CD4_Temra_GO <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/2_VIA_paper_data/2_Module_score_VIA/data/1_CD4_Temra_vs_CD4_Tem_GO_results.csv")
rownames(CD4_Temra_GO) <- CD4_Temra_GO$X 
 CD4_Temra_GO$X <- NULL
colnames(CD4_Temra_GO)
head(CD4_Temra_GO,n=1)

CD4_Treg_GO <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/2_VIA_paper_data/2_Module_score_VIA/data/2_CD4_Treg_vs_Others_GO_results.csv")
rownames(CD4_Treg_GO) <- CD4_Treg_GO$X 
CD4_Treg_GO$X <- NULL
colnames(CD4_Treg_GO)
head(CD4_Treg_GO,n=1)

head(CD4_Treg_GO$Description,n=10)

CD4T <- analyze_and_visualize_GO_genes (CD4_Temra_GO, CD4T,  ctrl_genes = NULL,
                                        Description_name="leukocyte mediated immunity",
                                        Celltype=  'celltype.l2',
                                         Condition = 'timepoint',
                                        column_names_rot =45,
                                         column_names_centered = FALSE,
                                   module_score_name = 'Temra_act', module_score_real_name = 'Temra_act1',
                                   output_heatmap_file = "figures/1_CD4_Temra_heatmap_GO_genes_Module_2.pdf",
                                   col_fun_values = c(0, 0.5, 1, 1.5, 2, 2.5), 
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 8, heatmap_height = 4)

CD4T <- analyze_and_visualize_GO_genes (CD4_Temra_GO, CD4T,  ctrl_genes = NULL,
                                        Description_name="leukocyte mediated immunity",
                                        Celltype=  'celltype.l2',
                                         Condition = 'day',
                                        column_names_rot =45,
                                         column_names_centered = FALSE,
                                   module_score_name = 'Temra_act', module_score_real_name = 'Temra_act1',
                                   output_heatmap_file = "figures/1_CD4_Temra_heatmap_GO_genes_Module_day.pdf",
                                   col_fun_values = c(0, 0.5, 1, 1.5, 2, 2.5), 
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 8, heatmap_height = 4)

#展示GO某个的基因Moudule_scor时期变化
CD4T <- analyze_and_visualize_GO_genes (CD4_Treg_GO, CD4T, 
                                        Description_name="negative regulation of T cell activation",
                                        ctrl_genes = NULL, 
                                        Celltype=  'celltype.l2',
                                         Condition = 'timepoint',
                                        column_names_rot =45,
                                         column_names_centered = FALSE,
                                   module_score_name = 'Treg_act', module_score_real_name = 'Treg_act1',
                                   output_heatmap_file = "figures/2_CD4_Treg_heatmap_GO_genes_Module.pdf",
                                    col_fun_values = c(0, 0.5, 0.8, 1, 1.25, 1.5),  
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 10, heatmap_height = 4)

#展示GO某个的基因Moudule_scor时期变化
CD4T <- analyze_and_visualize_GO_genes (CD4_Treg_GO, CD4T, 
                                        Description_name="negative regulation of T cell activation",
                                        ctrl_genes = NULL, 
                                        Celltype=  'celltype.l2',
                                         Condition = 'day',
                                        column_names_rot =45,
                                         column_names_centered = FALSE,
                                   module_score_name = 'Treg_act', module_score_real_name = 'Treg_act1',
                                   output_heatmap_file = "figures/2_CD4_Treg_heatmap_GO_genes_Module_day.pdf",
                                    col_fun_values = c(0, 0.5, 0.8, 1, 1.25, 1.5),  
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 10, heatmap_height = 4)



Idents(mRNA_obj) <- 'celltype.l1'
CD8T <- subset(mRNA_obj,idents=c('CD8 T'))
nrow(CD8T@meta.data)

CD8_NELL2_GO <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/2_VIA_paper_data/2_Module_score_VIA/data/2_CD8_NELL2_vs_CD8_Tn_GO_results.csv")
rownames(CD8_NELL2_GO) <- CD8_NELL2_GO$X 
CD8_NELL2_GO$X <- NULL
colnames(CD8_NELL2_GO)
#head(CD8_NELL2_GO,n=1)

CD8_Temra_GO <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/2_VIA_paper_data/2_Module_score_VIA/data/1_CD8_Temra_vs_CD8_Tem_GO_results.csv")
rownames(CD8_Temra_GO) <- CD8_Temra_GO$X 
 CD8_Temra_GO$X <- NULL
colnames(CD8_Temra_GO)
#head(CD8_Temra_GO,n=1)

head(CD8_Temra_GO$Description,n=10)

#展示GO某个的基因Moudule_scor时期变化
CD8T <- analyze_and_visualize_GO_genes (CD8_Temra_GO, CD8T, 
                                        Description_name="cell killing",
                                        ctrl_genes = NULL, 
                                        Celltype=  'celltype.l2',
                                         Condition = 'timepoint',
                                        column_names_rot =45,
                                         column_names_centered = FALSE,
                                   module_score_name = 'Temra_act', module_score_real_name = 'Temra_act1',
                                   output_heatmap_file = "figures/3_CD8_Temra_heatmap_GO_genes_Module.pdf",
                                    col_fun_values = c(0, 0.8, 1.2, 1.8, 2, 2.3),  
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 10, heatmap_height = 4)

#展示GO某个的基因Moudule_scor时期变化
CD8T <- analyze_and_visualize_GO_genes (CD8_NELL2_GO, CD8T, 
                                        Description_name="T cell differentiation",
                                        ctrl_genes = NULL, 
                                        Celltype=  'celltype.l2',
                                         Condition = 'timepoint',
                                        column_names_rot =45,
                                         column_names_centered = FALSE,
                                   module_score_name = 'NELL2_act', module_score_real_name = 'NELL2_act1',
                                   output_heatmap_file = "figures/3_CD8_NELL2_heatmap_GO_genes_Module.pdf",
                                    col_fun_values = c(0, 2, 3.5, 4, 4.5, 4.8),  
                                   col_fun_colors = c('#2166AC', '#92C5DE', '#fffde7', '#F4A582', '#FD8D3C', '#B2182B'),
                                   heatmap_width = 10, heatmap_height = 4)

colnames(mRNA_obj@meta.data)

write.csv(mRNA_obj@meta.data,'mRNA_obj_meta_data.csv')

unique(mRNA_obj$celltype.l1)

nclos_clonotype <-  c("lightgrey", "darkgrey", "#E0C5BEFF",
                      "#F39B7FB2", "#E57E7EFF", "#DE5C00FF", "#B22C2CFF", "#573333FF")

# 定义函数来处理细胞对象的克隆型数据
process_clonotype_data <- function(cell_type_obj,ct_seq) {
  cell_type_obj$celltype.l2 <- factor(
    cell_type_obj$celltype.l2,
    levels = ct_seq  # 设置 celltype.l2 因子顺序
  )
  cell_type_obj$donor_id <- factor(
    cell_type_obj$donor_id,
    levels = unique(cell_type_obj$donor_id)# 设置 donor_id 因子顺序
  ) 
  cell_type_obj$experiment <- factor(
    cell_type_obj$experiment,
    levels = unique(cell_type_obj$experiment)# 设置 experiment 因子顺序
  )
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

Idents(mRNA_obj) <- 'celltype.l1'
CD4T <- subset(mRNA_obj,idents=c('CD4 T'))
nrow(CD4T@meta.data)
unique(CD4T$celltype.l2)


dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/8_ATG7_NELL2'
CD4T <- read.csv(paste0(dir_path,'/data/3_mRNA_CD4T_Predicted.csv'))
nrow(query_meta.data)
table(CD4T$celltype.l2)
table(CD4T$predicted.id)

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
CD4T <- CD4T %>%
  mutate(Condition = recode(timepoint, !!!time_mapping))
# 查看映射结果
unique(CD4T[, c("timepoint", "Condition")])

#CD4T <- CD4T %>%filter(Condition != "C3")

# CD4T 调用：
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
#ct_seq <- c("Naive CD4 T", "CD4 T CM", "Treg", "CD4 T activated")
ct_seq <- c('CD4_Tn','CD4_Tcm','CD4_T_CCR6','CD4_Th2','CD4_Treg','CD4_Tem','CD4_Temra')
clonotype_num <- as.data.frame(table(CD4T$TCRab_TRB.CDR3.1))
CD4T$Clonotype_num <- clonotype_num$Freq[match(CD4T$TCRab_TRB.CDR3.1, clonotype_num$Var1)]
CD4T <- process_clonotype_data(CD4T,ct_seq)



color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
ct_seq <- c('CD4_Tn','CD4_Tcm','CD4_T_CCR6','CD4_Th2','CD4_Treg','CD4_Tem','CD4_Temra')

CD4T$predicted.id <- factor(CD4T$predicted.id,
                            levels=ct_seq)
data <- as.data.frame(table(CD4T$predicted.id, CD4T$Clonotype_num_ld))

colnames(data) <- c("predicted.id", "Clonotype_num_ld", "Count")

pdf("figures/6_CD4T_Bar_Clonotype_Levels_all_2.pdf",3.6,2.5)
ggplot(data, aes(x = predicted.id, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(
       x = "predicted.id", 
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
ggplot(data, aes(x = predicted.id, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_mapping) +
  labs(
       x = "predicted.id", 
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
CD4T_f <- CD4T%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD4T_f$Clonotype_num_ld <- factor(CD4T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5",  "n > 5")))
data <- as.data.frame(table(CD4T_f$predicted.id, CD4T_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("predicted.id", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/6_CD4T_Bar_Clonotype_Levels_2_predicted.id.pdf",3.6,2.5)
ggplot(data, aes(x = predicted.id, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number", 
       x = "predicted.id", 
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

# Idents(mRNA_obj) <- 'celltype.l1'
# CD8T <- subset(mRNA_obj,idents=c('CD8 T'))
# nrow(CD8T@meta.data)
# unique(CD8T$celltype.l2)

dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/8_ATG7_NELL2'
CD8T <- read.csv(paste0(dir_path,'/data/3_mRNA_CD8T_NELL2_Predicted.csv'))
nrow(query_meta.data)
table(CD8T$celltype.l2)
table(CD8T$predicted.id)

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
CD8T <- CD8T %>%
  mutate(Condition = recode(timepoint, !!!time_mapping))
# 查看映射结果
#unique(CD8T[, c("timepoint", "Condition")])
#CD8T <- CD8T %>%filter(Condition != "C3")

# CD8T 调用：
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
#ct_seq <- c("CD8 T Eff", "CD8/CD28 T EM", "Naive CD8 T", "CD8 T EM")
ct_seq <- c("CD8_Tn","CD8_NELL2", "CD8_Tem", "CD8_Temra")
clonotype_num <- as.data.frame(table(CD8T$TCRab_TRB.CDR3.1))
CD8T$Clonotype_num <- clonotype_num$Freq[match(CD8T$TCRab_TRB.CDR3.1, clonotype_num$Var1)]
CD8T <- process_clonotype_data(CD8T,ct_seq)


#CD8T
color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
CD8T$predicted.id <- factor(CD8T$predicted.id,level=ct_seq)
data <- as.data.frame(table(CD8T$predicted.id, CD8T$Clonotype_num_ld))
colnames(data) <- c("predicted.id", "Clonotype_num_ld", "Count")

pdf("figures/6_CD8T_Bar_Clonotype_Levels_all_predicted.id.pdf",3.6,2.5)
ggplot(data, aes(x = predicted.id, y = Count, fill = Clonotype_num_ld)) +
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
ggplot(data, aes(x = predicted.id, y = Count, fill = Clonotype_num_ld)) +
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
CD8T_f <- CD8T%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD8T_f$Clonotype_num_ld <- factor(CD8T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5",  "n > 5")))
data <- as.data.frame(table(CD8T_f$predicted.id, CD8T_f$Clonotype_num_ld))

# 将列重命名为更友好的名称
colnames(data) <- c("predicted.id", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/6_CD8T_Bar_Clonotype_Levels_predicted.id.pdf",3.6,2.5)
ggplot(data, aes(x = predicted.id, y = Count, fill = Clonotype_num_ld)) +
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





# 定义函数来处理细胞对象的克隆型数据
process_clonotype_data <- function(cell_type_obj,ct_seq) {
  cell_type_obj$celltype.l2 <- factor(
    cell_type_obj$celltype.l2,
    levels = ct_seq  # 设置 celltype.l2 因子顺序
  )
  cell_type_obj$donor_id <- factor(
    cell_type_obj$donor_id,
    levels = unique(cell_type_obj$donor_id)# 设置 donor_id 因子顺序
  ) 
  cell_type_obj$experiment <- factor(
    cell_type_obj$experiment,
    levels = unique(cell_type_obj$experiment)# 设置 experiment 因子顺序
  )
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

nclos_clonotype <-  c("lightgrey", "darkgrey", "#E0C5BEFF",
                      "#F39B7FB2", "#E57E7EFF", "#DE5C00FF", "#B22C2CFF", "#573333FF")


mRNA_obj <- read.csv('/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/1_data_info_and_module_scores/mRNA_obj_meta_data.csv')

unique(mRNA_obj$group)

colnames(mRNA_obj)

mRNA_obj <- mRNA_obj%>% filter(group %in% c('HC','Vax','booster'))

CD4T <- mRNA_obj%>% filter(celltype.l1 == 'CD4 T')
CD8T <- mRNA_obj%>% filter(celltype.l1 == 'CD8 T')

nrow(CD4T)
nrow(CD8T)

# CD4T 调用：
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
ct_seq <- c("Naive CD4 T", "CD4 T CM", "Treg", "CD4 T activated")
clonotype_num <- as.data.frame(table(CD4T$TCRab_TRB.CDR3.1))
CD4T$Clonotype_num <- clonotype_num$Freq[match(CD4T$TCRab_TRB.CDR3.1, clonotype_num$Var1)]
CD4T <- process_clonotype_data(CD4T,ct_seq)


# CD8T 调用：
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
ct_seq <- c("Naive CD8 T","CD8/CD28 T EM","CD8 T EM","CD8 T Eff")
clonotype_num <- as.data.frame(table(CD8T$TCRab_TRB.CDR3.1))
CD8T$Clonotype_num <- clonotype_num$Freq[match(CD8T$TCRab_TRB.CDR3.1, clonotype_num$Var1)]
CD8T <- process_clonotype_data(CD8T,ct_seq)

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(CD4T$celltype.l2, CD4T$Clonotype_num_ld))
colnames(data) <- c("celltype.l2", "Clonotype_num_ld", "Count")

pdf("figures/6_mRNA_CD4T_Bar_Clonotype_Levels_all.pdf",3.6,2.5)
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
CD4T_f <- CD4T@meta.data%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD4T_f$Clonotype_num_ld <- factor(CD4T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5",  "n > 5")))
data <- as.data.frame(table(CD4T_f$celltype.l2, CD4T_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("celltype.l2", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/6_mRNA_CD4T_Bar_Clonotype_Levels.pdf",3.6,2.5)
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

#CD8T
color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(CD8T$celltype.l2, CD8T$Clonotype_num_ld))
colnames(data) <- c("celltype.l2", "Clonotype_num_ld", "Count")

pdf("figures/6_mRNA_CD8T_Bar_Clonotype_Levels_all.pdf",3,2)
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
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1,color='black')# 设置 X 轴文本
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
CD8T_f <- CD8T%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD8T_f$Clonotype_num_ld <- factor(CD8T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5",  "n > 5")))
data <- as.data.frame(table(CD8T_f$celltype.l2, CD8T_f$Clonotype_num_ld))

# 将列重命名为更友好的名称
colnames(data) <- c("celltype.l2", "Clonotype_num_ld", "Count")
# 绘制堆叠条形图
pdf("figures/6_mRNA_CD8T_Bar_Clonotype_Levels.pdf",3,2.2)
ggplot(data, aes(x = celltype.l2, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(#title = "Clonotype Number", 
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
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1,color='black'), # 设置 X 轴文本
      legend.spacing = unit(0.05, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.15, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) 
dev.off()

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
library(ggplot2)

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
CD4T <- CD4T %>%
  mutate(Condition = recode(timepoint, !!!time_mapping))
summary(CD4T$timepoint)
summary(CD4T$Condition)
# 查看映射结果
unique(CD4T[, c("timepoint", "Condition")])


# 使用 dplyr::mutate 和 recode 应用映射
CD8T <- CD8T %>%
  mutate(Condition = recode(timepoint, !!!time_mapping))
summary(CD8T$timepoint)
summary(CD8T$Condition)
# 查看映射结果
unique(CD8T[, c("timepoint", "Condition")])

CD4T$Condition <- factor(CD4T$Condition,
                         levels=c('A0','A1','A2',
                                  'B0','B1','B2',
                                  'C0','C1','C2','C3'))
#,'acute','convalescent'
CD8T$Condition <- factor(CD8T$Condition,
                         levels=c('A0','A1','A2',
                                  'B0','B1','B2',
                                  'C0','C1','C2','C3'))

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(CD4T$Condition, CD4T$Clonotype_num_ld))
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/7_CD4T_Bar_Clonotype_Levels_Condition_all.pdf",3.7,1.6)
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
    axis.text.x = element_text(size = 8, color='black'),# 设置 X 轴文本
      legend.spacing = unit(0.05, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.15, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
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
    axis.text.x = element_text(size = 8,  color='black')# 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)
CD4T_f <- CD4T%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD4T_f$Clonotype_num_ld <- factor(CD4T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(CD4T_f$Condition, CD4T_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/7_CD4T_Bar_Clonotype_Levels_Condition.pdf",3.7,1.6)
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(#title = "Clonotype Number", 
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
    axis.text.x = element_text(size = 8, color='black'), # 设置 X 轴文本
       legend.spacing = unit(0.05, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.15, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) 
dev.off()

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(CD8T$Condition, CD8T$Clonotype_num_ld))
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/7_CD8T_Bar_Clonotype_Levels_Condition_all.pdf",3.7,1.6)
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
    axis.text.x = element_text(size = 8,  color='black'),# 设置 X 轴文本
      legend.spacing = unit(0.05, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.15, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
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
    axis.text.x = element_text(size = 8,  color='black')# 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)
CD8T_f <- CD8T%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD8T_f$Clonotype_num_ld <- factor(CD8T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(CD8T_f$Condition, CD8T_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/7_CD8T_Bar_Clonotype_Levels_Condition.pdf",3.7,1.6)
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(#title = "Clonotype Number", 
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
    axis.text.x = element_text(size = 8, color='black'), # 设置 X 轴文本 # angle = 45,hjust = 1,
      legend.spacing = unit(0.05, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.15, "inch"),
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) 
dev.off()

table(CD8T$predicted.id)
table(CD8T$predicted.id,CD8T$Clonotype_num_ld)

CD8_Temra_Tem <- CD8T%>% filter(predicted.id %in% c('CD8_Tem','CD8_Temra') )
nrow(CD8_Temra_Tem)

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(CD8_Temra_Tem$Condition, CD8_Temra_Tem$Clonotype_num_ld))
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/1_CD8_Temra_Tem_Bar_Clonotype_Levels_Condition_all.pdf",3.6,2.2)
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number in CD8", 
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
    axis.text.x = element_text(size = 8, color='black')# 设置 X 轴文本
  ) 
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number in CD8", 
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
CD4T_f <-CD8_Temra_Tem%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD4T_f$Clonotype_num_ld <- factor(CD4T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(CD4T_f$Condition, CD4T_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/1_CD8_Temra_Tem_Bar_Clonotype_Levels_Condition.pdf",3.6,2.2)
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number in CD8", 
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
    axis.text.x = element_text(size = 8, color='black'), # 设置 X 轴文本
  ) 
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity") + #position = "fill"
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number in CD8", 
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
    axis.text.x = element_text(size = 8, color='black'), # 设置 X 轴文本
  ) 
dev.off()

table(CD8_Temra_Tem$Clonotype_num_ld)

#估算下TCR richness在每个节点的变化，richness = unique（TCR）/all TCR
library(dplyr)
# 提取相关元数据
richness_df <- as.data.frame(CD8_Temra_Tem[ ,c("TCRab_TRB.clonotype.1","Condition","donor_id")]) %>%
  na.omit()
colnames(richness_df) <- c("tcr_clonotype_id_new","Condition","Donor")
head(richness_df)
nrow(richness_df)

#估算下TCR richness在每个节点的变化，richness = unique（TCR）/all TCR
library(dplyr)
# 提取相关元数据
richness_df <- as.data.frame(CD8_Temra_Tem[ ,c("TCRab_TRB.clonotype.1","Condition","donor_id")]) %>%
  na.omit()
colnames(richness_df) <- c("tcr_clonotype_id_new","Condition","Donor")
# 计算 unique TCR 数量
richness_unique <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())

# 计算总的 TCR 数量
richness_total <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])

# 合并并计算 richness
richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(richness$Var2)
richness$richness <- richness$unique_TCR / richness$total_TCR
# 整理列名
names(richness) <- c("Condition", "Donor", "unique_TCR", "total_TCR", "richness")
# 按 Condition 分组计算均值/中位数
richness$Condition <- factor(richness$Condition, level = unique(richness$Condition))
richness <- as.data.frame(richness) %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),
         median_richness = median(richness, na.rm = TRUE))
# 查看结果
head(richness)
# 导出为 CSV
write.csv(richness, "figures/2_richness_CD8_Tem_Temra.csv", row.names = FALSE)




library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

pdf("figures/2_CD8T_Temra_Tem_line_richness_C3.pdf", width = 3, height = 2.2)

# 去除 Condition 为 "C3" 的样本
filtered_richness <- richness #%>%
  #filter(Condition != "C3")

# 计算平均值和标准误
summary_df <- filtered_richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

# 绘图
ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#53A85F", size = 0.5) +
  geom_point(color = "#53A85F", size = 1.8) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.15, linetype = "dashed", alpha = 0.6, color = "#53A85F") +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  theme_bw(base_size = 10) +
  labs(title = "CD8+ Temra and Tem Richness",
       y = "Richness", x = "Condition")

dev.off()



# dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/8_ATG7_NELL2'
# CD4T <- read.csv(paste0(dir_path,'/data/3_mRNA_CD4T_Predicted.csv'))
# nrow(query_meta.data)
# table(CD4T$celltype.l2)
# table(CD4T$predicted.id)



table(CD4T$predicted.id)
table(CD4T$predicted.id,CD4T$Clonotype_num_ld)

CD4_Temra_Tem <- CD4T%>% filter(predicted.id %in% c('CD4_Tem','CD4_Temra') )
nrow(CD4_Temra_Tem)

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(CD4_Temra_Tem$Condition, CD4_Temra_Tem$Clonotype_num_ld))
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/8_CD4_Temra_Tem_Bar_Clonotype_Levels_Condition_all_Predicted.pdf",3.6,2.2)
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number in CD4", 
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
    axis.text.x = element_text(size = 8, color='black')# 设置 X 轴文本
  ) 
ggplot(data, aes(x = Condition, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number in CD4", 
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
CD4T_f <- CD4_Temra_Tem%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD4T_f$Clonotype_num_ld <- factor(CD4T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(CD4T_f$Condition, CD4T_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/8_CD4_Temra_Tem_Bar_Clonotype_Levels_Condition_Predicted.pdf",3.6,2.5)
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
    axis.text.x = element_text(size = 8, color='black'), # 设置 X 轴文本
  ) 
dev.off()

table(CD4_Temra_Tem$Clonotype_num_ld)

#估算下TCR richness在每个节点的变化，richness = unique（TCR）/all TCR
library(dplyr)
# 提取相关元数据
richness_df <- as.data.frame(CD4_Temra_Tem[ ,c("TCRab_TRB.clonotype.1","Condition","donor_id")]) %>%
  na.omit()
colnames(richness_df) <- c("tcr_clonotype_id_new","Condition","Donor")
head(richness_df)
nrow(richness_df)

#估算下TCR richness在每个节点的变化，richness = unique（TCR）/all TCR
library(dplyr)
# 提取相关元数据
richness_df <- as.data.frame(CD4_Temra_Tem[ ,c("TCRab_TRB.clonotype.1","Condition","donor_id")]) %>%
  na.omit()
colnames(richness_df) <- c("tcr_clonotype_id_new","Condition","Donor")
# 计算 unique TCR 数量
richness_unique <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  filter(Freq > 0) %>%
  summarize(unique_TCR = n())

# 计算总的 TCR 数量
richness_total <- as.data.frame(table(richness_df$tcr_clonotype_id_new, richness_df$Condition, richness_df$Donor)) %>%
  group_by(Var2, Var3) %>%
  mutate(total_TCR = sum(Freq))
richness_total <- unique(richness_total[, c(2, 3, 5)])

# 合并并计算 richness
richness <- left_join(richness_unique, richness_total)
richness$unique_TCR[is.na(richness$unique_TCR)] <- 0
richness$Var2 <- factor(richness$Var2)
richness$richness <- richness$unique_TCR / richness$total_TCR
# 整理列名
names(richness) <- c("Condition", "Donor", "unique_TCR", "total_TCR", "richness")
# 按 Condition 分组计算均值/中位数
richness$Condition <- factor(richness$Condition, level = unique(richness$Condition))
richness <- as.data.frame(richness) %>%
  group_by(Condition) %>%
  mutate(mean_richness = mean(richness, na.rm = TRUE),
         median_richness = median(richness, na.rm = TRUE))
# 查看结果
head(richness)
# 导出为 CSV
write.csv(richness, "figures/9_richness_CD4_Tem_Temra_Predicted.csv", row.names = FALSE)




library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)

pdf("figures/9_CD4T_Temra_Tem_line_richness_Predicted_C3.pdf", width = 3, height = 2.2)

# 去除 Condition 为 "C3" 的样本
filtered_richness <- richness #%>%
  #filter(Condition != "C3")

# 计算平均值和标准误
summary_df <- filtered_richness %>%
  group_by(Condition) %>%
  summarize(mean_richness = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            n = n(),
            sem = sd / sqrt(n))

# 绘图
ggplot(summary_df, aes(x = Condition, y = mean_richness, group = 1)) +
  geom_line(color = "#E95C39", size = 0.5) +
  geom_point(color = "#E95C39",size = 1.8) +
  geom_errorbar(aes(ymin = mean_richness - sem, ymax = mean_richness + sem),
                width = 0.15, linetype = "dashed", alpha = 0.6, color =  "#E95C39") +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  theme_bw(base_size = 10) +
  labs(title = "CD4+ Temra and Tem Richness",
       y = "Richness", x = "Condition")

dev.off()



unique(CD4T$celltype.l2)

colnames(CD4T@meta.data)

table(CD4T$celltype.l2,CD4T$Clonotype_num_ld)

Idents(CD4T) <- 'celltype.l2'
nrow(CD4T@meta.data)
CD4_Treg <- subset(CD4T,idents=c('Treg'))
nrow(CD4_Treg@meta.data)

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(CD4_Treg$Condition, CD4_Treg$Clonotype_num_ld))
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/8_CD4_Treg_Bar_Clonotype_Levels_Condition_all.pdf",3.6,2.5)
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
    axis.text.x = element_text(size = 8, color='black')# 设置 X 轴文本
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
    axis.text.x = element_text(size = 8,  color='black')# 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)
CD4T_f <- CD4T@meta.data%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
CD4T_f$Clonotype_num_ld <- factor(CD4T_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(CD4T_f$Condition, CD4T_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/8_CD4_Treg_Bar_Clonotype_Levels_Condition.pdf",3.6,2.5)
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
    axis.text.x = element_text(size = 8, color='black'), # 设置 X 轴文本
  ) 
dev.off()

table(mRNA_obj@meta.data$celltype.l1,mRNA_obj@meta.data$celltype.l2)

unique(mRNA_obj@meta.data$celltype.l1)
unique(mRNA_obj@meta.data$celltype.l2)

# 加载 dplyr 包
library(dplyr)

# 定义映射关系，包含 ILC 类别
cell_type_mapping <- c(
  'Classical Monocytes' = 'Monocyte',
  'CD8 T Eff' = 'CD8T',
  'MAIT' = 'MAIT',
  'Resting B' = 'Bcell',
  'Memory B' = 'Bcell',
  'Naive CD4 T' = 'CD4T',
  'Non-classical Monocytes' = 'Monocyte',
  'NK' = 'NK',
  'CD4 T CM' = 'CD4T',
  'CD8/CD28 T EM' = 'CD8T',
  'pDCs' = 'pDC',
  'Naive CD8 T' = 'CD8T',
  'CD8 T EM' = 'CD8T',
  'DCs' = 'cDC',
  'HSPC' = 'HPSC',
  'NK MThi' = 'NK',
  'NEAT1hi Monocytes' = 'Monocyte',
  'NKG7+ Monocytes' = 'Monocyte',
  'gdT' = 'gdT',
  'Treg' = 'CD4T',
  'CD4 T activated' = 'CD4T',
  'T Proliferating' = 'Proliferating',
  'B Exhausted' = 'Bcell',
  'Plasmablasts' = 'Bcell',
  'Activated B' = 'Bcell'
    #ILC和Tdn 缺少
)

# 使用 dplyr::mutate 和 recode 应用映射
mRNA_obj@meta.data <- mRNA_obj@meta.data %>%
  mutate(celltype_l1 = recode(celltype.l2, !!!cell_type_mapping))

# 查看映射结果
unique(mRNA_obj@meta.data[, c("celltype.l2", "celltype_l1")])

summary(mRNA_obj@meta.data$celltype_l1)

summary(mRNA_obj@meta.data$celltype.l1)












