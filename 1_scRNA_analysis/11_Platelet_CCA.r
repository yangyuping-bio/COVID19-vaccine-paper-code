# 0.环境准备 --------------------------------------------------------------------
library("tidyverse")
library("ggpubr")
library("hrbrthemes")
library("gridExtra")
options(tibble.width = Inf, tibble.print_max = 10000)
library("Matrix")
library("doParallel")
library(BiocParallel)
library(Seurat)
library(future)
library(scales)
library(SingleCellExperiment)
library(DropletUtils)
library(DoubletFinder)
library(reticulate)
library("GSEABase")

# 设置路径
# 设置路径
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/1_Platelet/1_Reumap"
data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/4_Cell_type_level_1_2"
# setwd()指定路径，setwd("~/project/")

#运行环境设置
source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
# set.seed(2)
# system("hostname") #计算机主机名
# use_python("~/miniconda3/envs/SC/bin/python", required = FALSE)
# use_condaenv("SC",conda="~/miniconda3/bin/conda")
# scr = import("scrublet")

targetcell <- readRDS(file.path(data_save_path, "4_integrated_level_1_2.rds"))

unique(targetcell$ct_level_1)

# ---- subset cell
Idents(targetcell) <- "ct_level_1"
targetcell <- subset(targetcell, idents = c("Platelet"))

targetcell@meta.data %>%
  dplyr::group_by(ct_level_1) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)
targetcell@meta.data %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num) 

#sample.integrated <- targetcell

# ---- PCA and clustering
DefaultAssay(sample.integrated) <- "integrated"
#sample.integrated <- RunPCA(sample.integrated, npcs = 20)
sample.integrated <- RunUMAP(sample.integrated,  reduction = "pca", dims = 1:20)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:20, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(0.1, 2, 0.1))

# 定义绘制UMAP图的函数
plot_umap <- function(meta_data_col, output_file,width, height,allcolors,lable_size) {
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
  umap <- obj@reductions$umap@cell.embeddings %>%
    as.data.frame() %>%
    cbind(cell_type = meta_data_col)
  # 计算cell_type的中位数
  cell_type_med <- umap %>%
    group_by(cell_type) %>%
    summarise(umap_1 = median(umap_1),umap_2 = median(umap_2))
  pdf(file.path(output_file), width, height)  # 设置输出PDF文件
  # 绘制基础的UMAP图
  p <- ggplot(umap, aes(x = umap_1, y = umap_2, color = cell_type)) +
    geom_point(size = 0.1, alpha = 1) +
    scale_color_manual(values = allcolors) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white'),
      plot.background = element_rect(fill = "white")
    ) +
    theme(
      legend.title = element_blank(),
      legend.key = element_rect(fill = 'white'),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1, 'cm')
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(umap$umap_1), y = min(umap$umap_2),
                     xend = min(umap$umap_1) + 3, yend = min(umap$umap_2)),
                 colour = "black", size = 1, arrow = arrow(length = unit(0.3, "cm"))) +
    geom_segment(aes(x = min(umap$umap_1), y = min(umap$umap_2),
                     xend = min(umap$umap_1), yend = min(umap$umap_2) + 3),
                 colour = "black", size = 1, arrow = arrow(length = unit(0.3, "cm"))) +
    annotate("text", x = min(umap$umap_1) + 1.5, y = min(umap$umap_2) - 1, label = "umap_1",
             color = "black", size = 3, fontface = "bold") +
    annotate("text", x = min(umap$umap_1) - 1, y = min(umap$umap_2) + 1.5, label = "umap_2",
             color = "black", size = 3, fontface = "bold", angle = 90) 

  print(p)
  p2<- p+
    geom_label_repel(aes(label = cell_type), size = lable_size, fontface = "bold", data = cell_type_med,
                     point.padding = unit(0.5, "lines")) +
    theme(legend.position = "none")  # 去掉legend
  print(p2)
  dev.off()
}

plot_umap <- function(meta_data_col, output_file, width, height, allcolors, lable_size) {
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  umap <- obj@reductions$umap@cell.embeddings %>%
    as.data.frame() %>%
    cbind(cell_type = meta_data_col)
  
  # 计算cell_type的中位数
  cell_type_med <- umap %>%
    group_by(cell_type) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
  
  pdf(file.path(output_file), width, height)  # 设置输出PDF文件
  
  # 绘制基础的UMAP图
  p <- ggplot(umap, aes(x = umap_1, y = umap_2, color = cell_type)) +
    geom_point(size = 0.1, alpha = 1) +
    scale_color_manual(values = allcolors) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white'),
      plot.background = element_rect(fill = "white")
    ) +
    theme(
      legend.title = element_blank(),
      legend.key = element_rect(fill = 'white'),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1, 'cm')
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  print(p)
  
  p2 <- p +
    geom_label_repel(
      aes(label = cell_type), 
      size = lable_size, 
      fontface = "bold", 
      data = cell_type_med,
      point.padding = unit(0.5, "lines")
    ) +
    theme(legend.position = "none")  # 去掉legend
  
  print(p2)
  
  dev.off()
}

sample.integrated <- readRDS("Platelet_Reumap.rds")

unique(sample.integrated$ct_level_2)

ncols <- c("#BD956A")
obj <- sample.integrated
plot_umap(obj@meta.data$ct_level_2, "1_Platelet_DimPlot.pdf",1.5, 1.5,ncols,4)

table(obj$Condition)

getwd()

head(meta_data.df )

meta_data.df <- read.csv('platelet_Pro_df.csv')
plot_data <- meta_data.df[,c('Condition','Celltype_prop','ct_level_4')]
head(plot_data)

library(tidyverse)
library(ggpubr)

pdf('1_platelet_bar_plot_simplified.pdf', 1.6, 1)
ggplot(plot_data, aes(x = Condition, y = Celltype_prop)) +
  geom_bar(stat = "identity",color="#BD956A", width = 0.75) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 6) +
  theme(
      panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
        axis.text = element_text(size = 4),
        legend.position = "none")

dev.off()

meta_data.df <- read.csv('platelet_Pro_df.csv')
pdf('1_platelet_bar_plot.pdf',1.5,1.2)

ncols <- c("#C1E6F3")
ggplot(meta_data.df, aes(x = Condition, y = Celltype_prop, group = factor(ct_level_4), fill = factor(ct_level_4))) +
  geom_bar(stat = "identity", position = position_dodge(),  width = 0.8) +
  geom_line(aes(group = factor(ct_level_4)), size = 0.3, color = "#F1BC62", position = position_dodge(width = 0.7)) +
  geom_point(size = 0.5, position = position_dodge(width = 0.5),shape=9,color = "#F1BC62") +
  scale_fill_manual(values = ncols) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 4.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        legend.position = "none"
       )
dev.off()

meta_data.df <- read.csv('platelet_Pro_df.csv')
ncols <- c("#BD956A")
pdf('1_platelet_line_plot.pdf',1.5,0.8)
 ggplot(meta_data.df, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(ct_level_4))) +
  geom_point(size=0.6,shape=9) +
  geom_line(size=0.4) +
  scale_color_manual(values = ncols) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "none"
       ) 
dev.off()

unique(sample.integrated$integrated_snn_res.1)

unique(sample.integrated$integrated_snn_res.0.5)

Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "Reumap_DimPlot_Donor.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Donor',
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "Reumap_DimPlot_Batch.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Batch',
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "Reumap_DimPlot_Sample.pdf"),25,25)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Sample',
             ncol = 6,label = T,raster = FALSE)
print(p)
dev.off()

allcolors <- c(
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
  "#968175", "#51A3CCFF", "#f18800","#85B22CFF","#B22C2CFF",
            "#E57E7EFF", "#BFB2FFFF","#e20612", "#ffd401", 
             "#E5B17EFF", "#942d8d","#573333FF", "#2e409a",
    "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",      
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",        
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",     
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887"
)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.0.1","integrated_snn_res.0.2",
                      "integrated_snn_res.0.3","integrated_snn_res.0.4",
                      "integrated_snn_res.0.5","integrated_snn_res.0.6",
                      "integrated_snn_res.0.7","integrated_snn_res.0.8",
                      "integrated_snn_res.0.9","integrated_snn_res.1",
                      "integrated_snn_res.1.5","integrated_snn_res.2",
                      "Donor","Condition","Sample","Batch","HTO_classification.global",
                      "scrublet_callrate","DF.classifications")
pdf(file.path(results_path, "DimPlot_all.pdf"),6,6)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.0.5","ct_level_2","ct_level_1")
pdf(file.path(results_path, "DimPlot_res.0.5.pdf"),6,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 4,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_B_cell.pdf"), width = 20, height=20)
FeaturePlot(object = sample.integrated ,
            features = c( "CD79A","CD79B","TCL1A","CD37","CD38",
                         "CD27","IGHG1","IGHM","IGHD",#Th22
                          "IGHA1", "IGHA2", 
                         "JCHAIN", "MS4A1","CD19",# plasma
                          "TCL1A","MZB1","HLA-A",
                         "SOX5", 
                         "CD3D" #检查是否混入T细胞
                        ),
            cols = c("grey", "blue"), pt.size = 0.1, 
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_others.pdf"), width = 20, height =20)
FeaturePlot(object = sample.integrated,
            features = c("CD3D","CD4","CD8A","SLC4A10",
                        "CCR7","GZMA","GZMK","GZMB",
                         "JCHAIN","CD1C", "PF4","FOXP3",
                         "TRDC","TRGV9", "TRGV10","TRDV1",
                         "CD79A","TCL1A","IGHA1","MZB1",
                         "CD27","IGHG1","IGHM","IGHD",
                         "CD14","FCGR3A","FCGR3B","LILRA4",
                         "CD68","CD34","CFP","STAB1","MKI67"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",
            raster=FALSE,label= TRUE,label.size = 3)
dev.off()

idents_name_list <- c("integrated_snn_res.0.5","Donor","Condition","Sample","Batch")
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score")
pdf(file.path(results_path, "Reumap_QC_vinplot.pdf"))
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  print(idents_name)
  for (plot_name in plot_list) {
    p <- VlnPlot(sample.integrated, features = c(plot_name), pt.size = 0, raster=FALSE) 
    print(p)
  }
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.0.5"
plan("multicore", workers = 4)
sample.integrated.markers <- FindAllMarkers(sample.integrated, only.pos = TRUE , min.pct = 0.15, logfc.threshold = 0.25)  #,min.pct = 0.1, logfc.threshold = 0.25
write.csv(sample.integrated.markers, file.path(results_path, "Markers_res.0.5.csv"))

table(sample.integrated$Condition)

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "Condition"
plan("multicore", workers = 4)
sample.integrated.markers <- FindAllMarkers(sample.integrated, only.pos = TRUE , min.pct = 0.15, logfc.threshold = 0.25)  #,min.pct = 0.1, logfc.threshold = 0.25
write.csv(sample.integrated.markers, file.path(results_path, "Condition_Markers.csv"))

saveRDS(sample.integrated, "Platelet_Reumap.rds")


