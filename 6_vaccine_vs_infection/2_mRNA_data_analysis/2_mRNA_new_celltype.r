# 加载必要的库
library(Seurat)
library(ComplexHeatmap)
library(cols4all)
library(GSEABase)
library(tidyverse)           # 包含 ggplot2, dplyr, tidyr, 等
library(org.Hs.eg.db)        # 加载人类基因注释数据库
library(clusterProfiler)
library(circlize)
library(ggpubr)              # 提供高级绘图功能
library(RColorBrewer)        # 配色工具
library(grid)
library(scales)
library(Matrix)
library(doParallel)
library(BiocParallel)
library(SingleCellExperiment)
library(reticulate)

source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
source("/share/home/qlab/projects/project_cvv/yyp_results_3/yyp_function.R") #加载到当前的R环境中

getwd()

mRNA_obj <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/data/GSE247917_cv19_combined_obj.rds")

mRNA_obj

sample.combined <- mRNA_obj

msigdb <- getBroadSets("/share/reference/sigdb/msigdb_v7.4.xml")
index <- sapply(msigdb, function(gs) {
  bcCategory(collectionType(gs)) == "c2"
})
geneset.c2 <- msigdb[index]
geneset.interferon <- geneset.c2[["BROWNE_INTERFERON_RESPONSIVE_GENES"]]
interferon.genes <- geneset.interferon@geneIds
interferon.genes <- interferon.genes[interferon.genes %in% rownames(sample.combined)]
cellcycle.genes <- geneset.c2[["KEGG_CELL_CYCLE"]]@geneIds
    
genes.hypoxia <- c(
  "VEGFA", "SLC2A1", "PGAM1", "ENO1",
  "LDHA", "TPI1", "P4HA1", "MRPS17",
  "CDKN3", "ADM", "NDRG1", "TUBB6",
  "ALDOA", "MIF", "ACOT7"
)    
stress.genes <- read_csv("/share/home/qlab/projects/project_cvv/data/stress.genes.csv") %>%
  dplyr::filter(gene %in% rownames(sample.combined)) %>%
  dplyr::pull(gene)
    
# HG genes
hgGenes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
# IG genes
igGenes <- c(
  "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "JCHAIN",
  "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"
)

IG.gene1 <- grep(pattern = "^IGLV", x = rownames(x = sample.combined), value = TRUE)
IG.gene2 <- grep(pattern = "^IGKV", x = rownames(x = sample.combined), value = TRUE)
IG.gene3 <- grep(pattern = "^IGHV", x = rownames(x = sample.combined), value = TRUE)
IG.genes <- union(IG.gene1, c(IG.gene2, IG.gene3))
# TRV genes
TRDV.genes <- grep(pattern = "^TRDV", x = rownames(x = sample.combined), value = TRUE)
TRAV.genes <- grep(pattern = "^TRAV", x = rownames(x = sample.combined), value = TRUE)
TRBV.genes <- grep(pattern = "^TRBV", x = rownames(x = sample.combined), value = TRUE)
TRGV.genes <- grep(pattern = "^TRGV", x = rownames(x = sample.combined), value = TRUE)
TRV.genes <- union(TRDV.genes, c(TRAV.genes, TRBV.genes, TRGV.genes))
# RPS & RPL
RPS.genes <- grep(pattern = "^RPS", x = rownames(x = sample.combined), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(x = sample.combined), value = TRUE)
# MT genes
MT.genes <- grep(pattern = "^MT-", x = rownames(x = sample.combined), value = TRUE)
# HIST genes
HIST.genes <- grep(pattern = "^HIST", x = rownames(x = sample.combined), value = TRUE)

DefaultAssay(sample.combined) <- "RNA"
sample.list <- SplitObject(sample.combined, split.by = "sampleName")
sample.list <- purrr::map(sample.list, function(x) {
  x <- NormalizeData(x, verbose = TRUE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features.select <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)

features.filtered <- setdiff(features.select, c(
  MT.genes, RPL.genes, RPS.genes, hgGenes, HIST.genes,IG.genes,
    TRV.genes, interferon.genes,
    stress.genes,cellcycle.genes,
     genes.hypoxia
))

sample.combined

colnames(sample.combined@meta.data)

# ---- scale and PCA for each sample
sample.list <- purrr::map(sample.list, function(x) {
  x <- ScaleData(x,
    vars.to.regress = c("nFeature_RNA"),
    features = features.filtered
  )
  x <- RunPCA(x, features = features.filtered, npcs = 30, verbose = FALSE)
})
#save(sample.list, file = "sample_list.rda")

sample.list[c(1,3,14,34)]

# ---- integrate samples
# load("sample_list.rda")
plan("multicore", workers = 4)
plan()
options(future.globals.maxSize = 100 * 1024^3)
gc()
plan()
sample.anchors <- FindIntegrationAnchors(
  sample.list,
  reduction = "cca", #reduction ="cca",
  # k.anchor = 5,
  # k.filter = 200,
  reference = c(1,3,14,34),
  dims = 1:30,
  anchor.features = features.filtered,
  scale = FALSE
)

plan("multicore", workers = 1)
plan()
options(future.globals.maxSize = 100 * 1024^3)
gc()
plan()

sample.integrated <- IntegrateData(
  anchorset = sample.anchors,
  dims = 1:30,
   k.weight = 50,
  verbose = TRUE
)

# ---- scale expression
sample.integrated <- FindVariableFeatures(sample.integrated, selection.method = "vst", nfeatures = 3000, verbose = F)
var.features <- VariableFeatures(sample.integrated)
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- ScaleData(
  sample.integrated,
  features = var.features
)

var.features.filtered <- var.features[!var.features %in% c(
     MT.genes, RPL.genes, RPS.genes, hgGenes, HIST.genes,IG.genes,
  interferon.genes,
    stress.genes,cellcycle.genes,
     genes.hypoxia
)]

# ---- PCA and clustering
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- RunPCA(sample.integrated, npcs = 30, features = var.features.filtered)

library(reticulate)
# 导入 Python 的 `pkg_resources` 模块
pkg_resources <- import("pkg_resources")
# 再次导入 `umap` 并手动注入 `pkg_resources` 属性
umap <- import("umap")
umap$pkg_resources <- pkg_resources

sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn",
                             reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(0.1, 3, 0.1))
saveRDS(sample.integrated, file = "1_mRNA_CCA_integrated.rds")

sample.integrated <- readRDS("1_mRNA_CCA_integrated.rds")

colnames(sample.integrated@meta.data)

# 删除低质量细胞数据
Idents(sample.integrated) <- 'celltype.l2'
sample.integrated <- subset(sample.integrated, 
                            idents = setdiff(levels(sample.integrated), 'NKG7+ Monocytes'))

### Reumap
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- RunPCA(sample.integrated, npcs = 30)
sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn",
                             reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(3, 10, 1))

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.4","integrated_snn_res.5",
                     "integrated_snn_res.6","integrated_snn_res.7",
                      "integrated_snn_res.8","integrated_snn_res.9",
                     "integrated_snn_res.10")
pdf(file.path(results_path, "1_CCA_DimPlot_10_first.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("celltype.l1",
                     "celltype.l2")
pdf(file.path(results_path, "1_CCA_DimPlot_celltype_first.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

unique(sample.integrated$celltype.l2)

# 删除低质量细胞数据
Idents(sample.integrated) <- "integrated_snn_res.4"
sample.integrated <- subset(sample.integrated, 
                            idents = setdiff(levels(sample.integrated), 
                                            c(51,63,66,53,61)))

unique(sample.integrated$integrated_snn_res.4)

### Reumap
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- RunPCA(sample.integrated, npcs = 30)
sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn",
                             reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(1, 10, 1))

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.1","integrated_snn_res.2",
                     "integrated_snn_res.3",
                      "integrated_snn_res.4","integrated_snn_res.5",
                     "integrated_snn_res.6","integrated_snn_res.7",
                      "integrated_snn_res.8","integrated_snn_res.9",
                     "integrated_snn_res.10")
pdf(file.path(results_path, "1_CCA_DimPlot_10_second.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

# 删除低质量细胞数据
Idents(sample.integrated) <- "integrated_snn_res.6"
sample.integrated <- subset(sample.integrated, 
                            idents = setdiff(levels(sample.integrated), 
                                            c(82,79)))

### Reumap
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- RunPCA(sample.integrated, npcs = 30)
sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn",
                             reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(1, 10, 1))

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.1","integrated_snn_res.2",
                     "integrated_snn_res.3",
                      "integrated_snn_res.4","integrated_snn_res.5",
                     "integrated_snn_res.6","integrated_snn_res.7",
                      "integrated_snn_res.8","integrated_snn_res.9",
                     "integrated_snn_res.10")
pdf(file.path(results_path, "1_CCA_DimPlot_10_third.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()
DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("celltype.l1",
                     "celltype.l2")
pdf(file.path(results_path, "1_CCA_DimPlot_celltype_third.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

# 删除低质量细胞数据
Idents(sample.integrated) <- "integrated_snn_res.10"
sample.integrated <- subset(sample.integrated, 
                            idents = setdiff(levels(sample.integrated), 
                                            c(117,124,128,130)))

### Reumap
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- RunPCA(sample.integrated, npcs = 30)
sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn",
                             reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(1, 10, 1))

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.1","integrated_snn_res.2",
                     "integrated_snn_res.3",
                      "integrated_snn_res.4","integrated_snn_res.5",
                     "integrated_snn_res.6","integrated_snn_res.7",
                      "integrated_snn_res.8","integrated_snn_res.9",
                     "integrated_snn_res.10")
pdf(file.path(results_path, "1_CCA_DimPlot_10_final.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()
DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("celltype.l1",
                     "celltype.l2")
pdf(file.path(results_path, "1_CCA_DimPlot_celltype_final.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

getwd()
results_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/6_mRNA_data/4_new_celltype/figures'

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("celltype.l1",
                      "celltype.l2",
                     "integrated_snn_res.2")
pdf(file.path(results_path, "1_CCA_Celltype_DimPlot.pdf"),13,6)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

idents_name_list <- c("integrated_snn_res.2","")
plot_list <- c("nFeature_RNA")

pdf(file.path(results_path, "1_CCA_QC_vinplot.pdf"))
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  print(idents_name)
  for (plot_name in plot_list) {
    p <- VlnPlot(sample.integrated, features = c(plot_name), pt.size = 0, raster=FALSE) 
    print(p)
  }
}
dev.off()

Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "2_CCA_FeaturePlot.pdf"), width = 80, height =100)
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

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "3_FeaturePlot_T_cell.pdf"), width = 40, height=60)
FeaturePlot(object = sample.integrated ,
            features = c("CD3D","CD3E","CD4","CD8A", # T_cell
                         "CCR7","CD27", # naive_T
                          "TRDV2","TRGV9", "TRGV10","TRDV1",# gdT
                          "CD25","FOXP3","CTLA4", # Treg
                         "ICOS", # CD4_Tfh
                         "IL17A",# Th17
                         "IFNG","CXCR3","CCR5", #Th1
                          "GATA3", #Th2
                        "TNFRSF4","NKG7","GNLY","SELL", # 细胞毒性,判断发育状态:Tcm,Tem,Temra
                         "GZMA","GZMK","GZMB","GZMH","GZMM",# 细胞毒性,判断发育状态:Tcm,Tem,Temra
                         "MK167", #pro T，增殖型T
                         "PRF1", #活化T
                         "NELL2",#CD8_NELL2
                         "BTLA", "NEDD4", #Anergic T
                         "FCGR3A","FGFBP2","CX3CR1",#Temra
                          "SLC4A10" # MAIT
                        ),
            cols = c("grey", "blue"), pt.size = 0.1, 
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "3_FeaturePlot_B_cell.pdf"), width = 40, height=40)
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

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "3_FeaturePlot_Mono_cell.pdf"), width = 40, height=40)
FeaturePlot(object = sample.integrated ,
            features = c("CD3D","CD4","CD8A", #T_cell
                         "CD14","FCGR3A","ATG7", #Mono
                         "CST3", "LYZ","CD68", "CD163",  #Mono
                         "SELL","S100A8","S100A12","CLEC12A","CXCL14",
                         "CX3CR1","KLF3","FLI1",#"IFTM1","FAM65B",
                         "CD300E","CEBPD","ZFP36L2",
                         "FCN1","APOBEC3A","THBS1",
                         "MARCO","CD33",
                         "PRF1", "GNLY", "KLRC4","KLRK1"
                        ),
            cols = c("grey", "blue"), pt.size = 0.1, 
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "3_FeaturePlot_Mono_others.pdf"), width = 40, height=40)
FeaturePlot(object = sample.integrated ,
            features = c("CD3D","CD4","CD8A", #T_cell
                         "CD14","FCGR3A","ATG7", #Mono
                         "CST3", "LYZ","CD68", "CD163",  #Mono
                         "SELL","S100A8","S100A12","CLEC12A","CXCL14",
                         "CX3CR1","KLF3","FLI1",#"IFTM1","FAM65B",
                         "CD300E","CEBPD","ZFP36L2",
                         "FCN1","APOBEC3A","THBS1",
                         "MARCO","CD33",
                         "PRF1", "GNLY", "KLRC4","KLRK1"
                        ),
            cols = c("grey", "blue"), pt.size = 0.1, 
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "3_FeaturePlot_Others_cell.pdf"), width = 40, height=30)
FeaturePlot(object = sample.integrated ,
            features = c("PF4","PPBP", # platelet
                         "CLEC9A", "XCR1", #cDC1
                         "CD1C","CLEC10A", #cDC2
                         "CCL19","CCR7", #cDC3
                         "CD34" #HPSC
                        ),
            cols = c("grey", "blue"), pt.size = 0.1, 
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
plan("multicore", workers = 2)
sample.integrated.markers <- FindAllMarkers(sample.integrated, only.pos = TRUE , min.pct = 0.1, logfc.threshold = 0.05)  #,min.pct = 0.1, logfc.threshold = 0.25
write.csv(sample.integrated.markers, file.path(results_path, "3_Markers_res.2.csv"))



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
  #'NKG7+ Monocytes' = 'Monocyte',
  'gdT' = 'gdT',
  'Treg' = 'CD4T',
  'CD4 T activated' = 'CD4T',
  'T Proliferating' = 'Proliferating',#需要再细分
  'B Exhausted' = 'Bcell',
  'Plasmablasts' = 'Bcell',
  'Activated B' = 'Bcell'
    #ILC和Tdn 缺少
)

# 使用 dplyr::mutate 和 recode 应用映射
sample.integrated@meta.data <- sample.integrated@meta.data %>%
  mutate(celltype_l1 = recode(celltype.l2, !!!cell_type_mapping))
summary(sample.integrated@meta.data$celltype_l1)
summary(sample.integrated@meta.data$celltype.l1)
# 查看映射结果
unique(sample.integrated@meta.data[, c("celltype.l2", "celltype_l1")])


unique(sample.integrated$celltype_l1)

# 删除celltype_l1为NA的细胞数据
Idents(sample.integrated) <- 'celltype_l1'
sample.integrated <- subset(sample.integrated, 
                            idents = c('Monocyte','MAIT','Bcell',
                                       'CD4T','CD8T',
                                       'pDC','HPSC','NK','gdT',
                                       'cDC','Proliferating'))

unique(sample.integrated$celltype.l2)

#----细胞注释---------
sample.integrated@meta.data$ct_merge  <- sample.integrated$celltype_l1
sample.integrated$ct_merge[sample.integrated$celltype.l1 == 'NK Proliferating' ] <- "NK"
sample.integrated$ct_merge[sample.integrated$celltype_l1 == 'Proliferating' ] <- "ILC"
sample.integrated$ct_merge[sample.integrated$integrated_snn_res.1 == 19 ] <- "NK"
sample.integrated$ct_merge[sample.integrated$integrated_snn_res.10 == 69 ] <- "ILC"

table(sample.integrated$ct_merge)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("ct_merge","celltype.l1","celltype.l2")
pdf(file.path(results_path, "2_CCA_ct_merge_DimPlot.pdf"),12,6)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

saveRDS(sample.integrated, file = "1_mRNA_final_CCA_integrated.rds")

## 在第七部更新ct_merge,更新图

sample.integrated$ct_merge <- factor(sample.integrated$ct_merge,
                                    levels=rev(c('CD4T','CD8T','MAIT','gdT','Tdn',
                                              'NK','Bcell','Monocyte','cDC','pDC',
                                              'ILC','HPSC')))

sample.integrated@meta.data$timepoint <- factor(
  sample.integrated@meta.data$timepoint,
  levels = c('HC','Vax_d7','Vax_d10','Vax_d14','Vax_d20','Vax_d21',
             'Vax_d28','Vax_d29','Vax_d35','Vax_d36',
             'booster_d0','booster_d7',
             'booster_d28','booster_d120',             
             'acute','cv19_d11','convalescent'))

table(sample.integrated@meta.data$timepoint)

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

# add moudule score
DefaultAssay(sample.integrated) <- "RNA"
genes.ifn <- c("IFI35","IFI44L","IFI6","IFIT3",
                        "IRF7","ISG15","MX1","MX2",
                         "OAS1","OAS2")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.ifn),
  ctrl = length(genes.ifn),
  name = "ifn_stim",
  assay = "RNA"
)

# add moudule score
DefaultAssay(sample.integrated) <- "RNA"
genes.HLA<- c("HLA-DQA2","HLA-DQA2")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.HLA),
  ctrl = length(genes.HLA),
  name = "HLA_DQA2",
  assay = "RNA"
)

Idents(sample.integrated) <- 'group'
cv_19_part <- subset(sample.integrated,idents=c('cv_19'))
mRNA_part <- subset(sample.integrated,idents=c('HC','Vax','booster'))

sample.integrated$state <- 'mRNA'
sample.integrated$state[sample.integrated$group == 'cv_19'] <- 'cv_19'
sample.integrated$state[sample.integrated$group == 'Vax'] <- 'mRNA'
sample.integrated$state[sample.integrated$group == 'HC'] <- 'mRNA'
sample.integrated$state[sample.integrated$group == 'booster'] <- 'mRNA'
table(sample.integrated$state)

getwd()

write.csv(sample.integrated@meta.data,  "data/1_mRNA_meta_data_df.csv")

# 加载 scales 包
library(scales)
# 对每个 state 内部的数据进行归一化
average_scores <- data.frame(sample.integrated@meta.data) %>%
  group_by(state, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(ifn_stim1, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(ifn_stim1 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

average_scores <- average_scores %>%
  group_by(state) %>% 
  mutate(mean_score_normalized = rescale(mean_score, to = c(-1, 1))) # 在每个 state1 内归一化
head(average_scores)

# 绘制热图
pdf("figures/4_Merge_ct_merge_IFN_module_Dotplot.pdf", width = 4.5, height = 1.8)
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
  facet_wrap(~state, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放
dev.off()
RdBu_new <- c('#2166AC','#4393C3','#92C5DE','#F4A582','#D6604D','#B2182B')
rect.data <- data.frame(
  state = factor(c("mRNA", "cv_19"), levels = c("mRNA", "cv_19")),
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue")  # 将颜色列表缩短为与 state 数量一致
)

# 将 average_scores 的 state 设置为因子，确保绘图时顺序正确
average_scores$state <- factor(average_scores$state, 
                                      levels = c("mRNA", "cv_19"))

# 绘制热图
pdf("figures/4_Merge_ct_merge_IFN_Dotplot_with_background.pdf", width = 4.5, height = 1.68)

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
  facet_wrap(~state, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()


# 加载 scales 包
library(scales)
# 对每个 state 内部的数据进行归一化
average_scores <- data.frame(sample.integrated@meta.data) %>%
  group_by(state, ct_merge, Condition) %>%
  summarize(
    mean_score = mean(HLA_DQA21, na.rm = TRUE), # 计算每组的平均得分
    percentage_expressing = mean(HLA_DQA21 > 0, na.rm = TRUE) * 100 # 计算表达该特征的细胞百分比
  )

average_scores <- average_scores %>%
  group_by(state) %>% 
  mutate(mean_score_normalized = rescale(mean_score, to = c(-1, 1))) # 在每个 state1 内归一化
head(average_scores)


# 绘制热图
pdf("figures/4_Merge_ct_merge_HLA_DQA21_module_Dotplot.pdf", width = 4.5, height = 1.8)
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
  facet_wrap(~state, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放
dev.off()

RdBu_new <- c('#2166AC','#4393C3','#92C5DE','#F4A582','#D6604D','#B2182B')
rect.data <- data.frame(
  state = factor(c("mRNA", "cv_19"), levels = c("mRNA", "cv_19")),
  ymin = -Inf,
  ymax = Inf,
  colors = c("red", "blue")  # 将颜色列表缩短为与 state 数量一致
)

# 将 average_scores 的 state 设置为因子，确保绘图时顺序正确
average_scores$state <- factor(average_scores$state, 
                                      levels = c("mRNA", "cv_19"))

# 绘制热图
pdf("figures/4_Merge_ct_merge_HLA_DQA21_Dotplot_with_background.pdf", width = 4.5, height = 1.68)

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
  facet_wrap(~state, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

dev.off()

merged_data <- data.frame(sample.integrated@meta.data) %>%
  group_by(state,Condition) %>%
  mutate(Count_Condition = n()) %>%
  ungroup%>% 
  group_by(state,Condition, ct_merge) %>%
  mutate(Count_Condition_celltype = n(),Celltype_prop = Count_Condition_celltype/Count_Condition) 
head(merged_data)

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

#celltypel2
pdf("figures/5_Line_Condition_ct_merge_free.pdf",  14, 14)
 ggplot(merged_data, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(ct_merge))) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_color_manual(values = colorl2) +
  facet_wrap(state~ ct_merge, scales = "free", ncol = 4) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
       ) 
dev.off()

library(dplyr)
# 定义基准条件
baseline_Conditions <- c("mRNA" = "A0",
                         "cv_19" = "acute")
# 对所有 state 进行循环计算
final_data <- lapply(names(baseline_Conditions), function(status) {
  # 筛选出对应的 state
  status_data <- merged_data %>%
    filter(state == status) %>%
    select(-ifn_stim1, -HLA_DQA21, -Count_Condition, -Count_Condition_celltype)
  # 获取对应的基准条件
  baseline_Condition <- baseline_Conditions[status]
  # 计算基准值（Baseline）
  baseline <- status_data %>%
    filter(Condition == baseline_Condition) %>%
    group_by(ct_merge) %>%
    summarise(Baseline_prop = mean(Celltype_prop, na.rm = TRUE))  
  # 将基准值合并回原数据
  status_data <- status_data %>%
    left_join(baseline, by = "ct_merge") %>%
    group_by(ct_merge) %>%
    # 计算 p-value
    mutate(p_value = ifelse(
      n_distinct(Condition) > 1 & sum(Condition == baseline_Condition) > 1 & sum(Condition != baseline_Condition) > 1,
      
      # 使用 filter 筛选基准条件和其他条件的数据进行 t 检验
      t.test(filter(., Condition == baseline_Condition)$Celltype_prop, 
             filter(., Condition != baseline_Condition)$Celltype_prop)$p.value, NA )) %>%
    ungroup() %>%
    # 计算 FoldChange 和 FDR
    mutate(FoldChange_baseline = Celltype_prop / Baseline_prop,
           FDR = p.adjust(p_value, method = "BH")) %>%
    # 移除不需要的列
    select(-Baseline_prop) %>%
    # 保留 state 以便最后合并
    mutate(state = status)
  return(status_data)
})

# 将所有结果合并为一个数据集
final_merged_data <- bind_rows(final_data)

# 查看合并后的数据
head(final_merged_data)

final_merged_data$ct_merge <- factor(final_merged_data$ct_merge,
                               levels = rev(c('CD4T','CD8T','MAIT','gdT','Tdn',
                                              'NK','Bcell','Monocyte',
                                              'cDC','pDC','ILC','HPSC')))
unique(final_merged_data$ct_merge)
final_merged_data$state <- factor(final_merged_data$state,
                                level=c('mRNA','cv_19'))
unique(final_merged_data$state)


write.csv(final_merged_data,"data/5_final_merged_data_FC.csv")

unique(final_merged_data$Condition)

# 创建热图数据
# 确保 FoldChange_baseline 是连续的
final_merged_data$FoldChange_baseline <- as.numeric(final_merged_data$FoldChange_baseline) 
plot_merged_data <- final_merged_data 
unique(plot_merged_data$Condition)

head(plot_merged_data)

summary(plot_merged_data$FoldChange_baseline)

library(ggplot2)
library(RColorBrewer)


# 选择调色板，这里使用 RdBu 调色板的颜色
#RdBu_pro <- rev(brewer.pal(11, "RdBu"))
RdBu_pro <- c('#053061','#2166AC','#F7F7F7','#FDDBC7','#F4A582','#D6604D','#B2182B','#67001F')

# 绘制热图
pdf("figures/5_Merge_ct_merge_Pro_Heatmap_with_A0_width.pdf", width = 3.8, height = 1.72)

heatmap_plot <- ggplot(plot_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_tile(aes(fill = FoldChange_baseline)) + # 使用 geom_tile 绘制热图
  scale_fill_gradientn(colors = RdBu_pro, 
                       name = "Fold change over baseline",
                       limits = c(0, 2.5), 
                    oob = scales::squish  # 将超出范围的数据压缩到边界值
                      ) + 
   # 设置连续的颜色渐变 #
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
  facet_wrap(~state, ncol = 5, scales = "free_x") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

print(heatmap_plot)
dev.off()

library(ggplot2)
library(RColorBrewer)

# 选择调色板，这里使用 RdBu 调色板的颜色
#RdBu_pro <- rev(brewer.pal(11, "RdBu"))
RdBu_pro <- c('#053061','#2166AC','#F7F7F7','#FDDBC7','#F4A582','#D6604D','#B2182B','#67001F')

# 绘制热图
pdf("figures/5_Merge_ct_merge_Pro_Heatmap_with_A0_width_2.pdf", width = 2.85, height = 1.72)

heatmap_plot <- ggplot(plot_merged_data, aes(x = Condition, y = ct_merge)) +
  geom_tile(aes(fill = FoldChange_baseline)) + # 使用 geom_tile 绘制热图
  scale_fill_gradientn(colors = RdBu_pro,
                       name = "Fold change over baseline",
                       limits = c(0, 2.5),
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
  facet_wrap(~state, ncol = 5, scales = "free") # scales 设置为 "free"，让 x 和 y 都可以自由缩放

print(heatmap_plot)
dev.off()



data_for_pie <- merged_data
head(data_for_pie)

library(ggplot2)
library(dplyr)
library(tidyr)

# 假设你的数据框名为 final_merged_data
# 计算每个 state 的细胞类型占比
data_for_pie <- data_for_pie %>%
  group_by(state, ct_merge) %>%
  summarise(count = n()) %>%
  group_by(state) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
data_for_pie$ct_merge <- factor(data_for_pie$ct_merge,
                               levels = c('CD4T','CD8T','MAIT','gdT','Tdn','NK','Bcell','Monocyte',
                                          'cDC','pDC','ILC','HPSC','Platelet'))
unique(data_for_pie$ct_merge)
head(data_for_pie)

data_for_pie$state <- factor(data_for_pie$state,
                                level=c('mRNA',
                                        'cv_19'))

library(ggplot2)

# 绘制饼图
pdf("6_Merge_ct_merge_Pie_big.pdf", 15, 15)

ggplot(data_for_pie, aes(x = "", y = percentage, fill = ct_merge)) +
  geom_bar(width = 1, stat = "identity", alpha = 1, color = "white") +  # 设置透明度为0.75
  coord_polar(theta = "y") +
  facet_wrap(~ state, ncol = 3) +  # 设置排列为5列
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
    plot.title = element_text(size = 28, hjust = 0.5, vjust = -86, face = "bold")  # 总图标题的字体大小，并调整位置
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))  # 设置图例为2行排列

dev.off()









sample.integrated<- readRDS("1_mRNA_final_CCA_integrated.rds")

colnames(sample.integrated@meta.data)

table(sample.integrated$celltype.l1)

table(sample.integrated$celltype.l2)

table(sample.integrated$celltype.l1,sample.integrated$celltype.l2)

# 提取 T Proliferating 细胞
Idents(sample.integrated) <- 'celltype.l2'
t_proliferating <- subset(sample.integrated, idents = "T Proliferating")
nrow(t_proliferating@meta.data)

t_proliferating

DefaultAssay(t_proliferating) <- "RNA"
# 标准化和 PCA
t_proliferating <- NormalizeData(t_proliferating)
t_proliferating <- FindVariableFeatures(t_proliferating)
t_proliferating <- ScaleData(t_proliferating)
t_proliferating <- RunPCA(t_proliferating)

# 聚类和 UMAP
t_proliferating <- FindNeighbors(t_proliferating, dims = 1:10)
t_proliferating <- FindClusters(t_proliferating, resolution = seq(0.1, 1, 0.1)) # 调整分辨率以获得合适的细分
t_proliferating <- RunUMAP(t_proliferating, dims = 1:10)

results_path <- getwd()

DefaultAssay(t_proliferating) <- "RNA"                         
idents_name_list <- c('RNA_snn_res.0.1','RNA_snn_res.0.2',
                      'RNA_snn_res.0.3','RNA_snn_res.0.4',
                      'RNA_snn_res.0.5','RNA_snn_res.0.6',
                      'RNA_snn_res.0.7','RNA_snn_res.0.8',
                      'RNA_snn_res.0.9','RNA_snn_res.1','celltype_l1',
                      'ct_merge','celltype.l1','celltype.l2')
pdf(file.path(results_path, "/figures/7_CCA_DimPlot.pdf"),4,2.5)
for (idents_name in idents_name_list) {
  Idents(t_proliferating) <- idents_name
  p <- DimPlot(object = t_proliferating, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()


table(t_proliferating$RNA_snn_res.0.3)

Idents(t_proliferating) <- "RNA_snn_res.0.3"
pdf(file.path("figures/7_FeaturePlot_T_proliferating_T_genes.pdf"), width = 12, height=15)
FeaturePlot(object = t_proliferating ,
            features = c("CD3D","CD3E","CD4","CD8A", # T_cell
                         "CCR7","CD27", # naive_T
                          "TRDV2","TRGV9", "TRGV10","TRDV1",# gdT
                          "CD25","FOXP3","CTLA4", # Treg
                         "ICOS", # CD4_Tfh
                         "IL17A",# Th17
                         "IFNG","CXCR3","CCR5", #Th1
                          "GATA3", #Th2
                        "TNFRSF4","NKG7","GNLY","SELL", # 细胞毒性,判断发育状态:Tcm,Tem,Temra
                         "GZMA","GZMK","GZMB","GZMH","GZMM",# 细胞毒性,判断发育状态:Tcm,Tem,Temra
                         "MK167", #pro T，增殖型T
                         "PRF1", #活化T
                         "NELL2",#CD8_NELL2
                         "BTLA", "NEDD4", #Anergic T
                         "FCGR3A","FGFBP2","CX3CR1",#Temra
                          "SLC4A10" # MAIT
                        ),
            cols = c("grey", "blue"), pt.size = 0.1, 
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

DefaultAssay(t_proliferating) <- "RNA"
Idents(t_proliferating) <- "RNA_snn_res.0.3"
plan("multicore", workers = 5)
t_proliferating.markers <- FindAllMarkers(t_proliferating, only.pos = TRUE , min.pct = 0.15, logfc.threshold = 0.25)  #,min.pct = 0.1, logfc.threshold = 0.25
write.csv(t_proliferating.markers, file.path("data/7_T_proliferating_RNA_snn_res_0_3.csv"))

unique(sample.integrated$ct_merge)

# 根据 Marker 基因重新命名
t_proliferating$SubClusters <- t_proliferating$RNA_snn_res.0.3
Idents(t_proliferating) <- "RNA_snn_res.0.3"
t_proliferating <- RenameIdents(t_proliferating, 
                                "0" = "Tdn", 
                                "1" = "Tdn", 
                                "2" = "CD8T",
                               "3" = "NK")  # 根据实际结果调整命名
# 提取细分群信息
t_proliferating@meta.data$SubClusters  <- Idents(t_proliferating)
t_proliferating@meta.data$SubClusters  <- as.character(t_proliferating@meta.data$SubClusters)

table(t_proliferating$SubClusters)

# 将细分群信息映射回原始数据
sample.integrated$SubClusters <- NA
sample.integrated$SubClusters[Cells(t_proliferating)] <- t_proliferating$SubClusters

# 验证分群是否添加成功
table(sample.integrated$SubClusters, sample.integrated$celltype.l2)

# 更新ct_merge
sample.integrated$ct_merge[sample.integrated$SubClusters == 'NK'] <- 'NK'
sample.integrated$ct_merge[sample.integrated$SubClusters == 'CD8T'] <- 'CD8T'
sample.integrated$ct_merge[sample.integrated$SubClusters == 'Tdn'] <- 'Tdn'
table(sample.integrated$ct_merge)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("ct_merge","celltype.l1","celltype.l2")
pdf(file.path("figures/7_CCA_ct_merge_DimPlot.pdf"),12,6)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

saveRDS(sample.integrated, file = "1_mRNA_final_CCA_integrated.rds")

sample.integrated

library(ggplot2)
library(ggrepel)
library(dplyr)
# 定义绘制UMAP图的函数
plot_umap <- function(meta_data_col, output_file,width, height,allcolors,lable_size) {
  umap <- sample.integrated@reductions$umap@cell.embeddings %>%
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

colpalette <- c(
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
 # Platelet="#BD956A",
  Tdn="#495608"# 其他细胞类型
)

plot_umap(sample.integrated@meta.data$ct_merge, 
          "figures/7_Umap_ct_merge.pdf",12,8,colpalette,12)



sample.integrated<- readRDS("1_mRNA_final_CCA_integrated.rds")

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

1

# df_genes <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/2_Mono/7_diff_gene_in_Mono/data/macroautophagy_genes.csv")
# macroautophagy_genes <- df_genes$X
# length(macroautophagy_genes)

# # add moudule score
# DefaultAssay(sample.integrated) <- "RNA"
# sample.integrated <- AddModuleScore(
#   object = sample.integrated,
#   features = list(macroautophagy_genes),
#   ctrl = length(macroautophagy_genes),
#   name = "mac_ATG7",
#   assay = "RNA"
# )

macroautophagy_genes <-c("ATG3", "ATG10", "ATG12", "ATG5", "ATG4",
                         "MAP1LC3B", "ATG16L1", "BECN1", "ULK1")
# c("ATG7","ATG3", "ATG10", "ATG12", "ATG5", "ATG4", "MAP1LC3B",
#                           "ATG16L1", "BECN1", "ULK1", "ATG9A", "ATG9B", "WIPI1",
#                           "WIPI2", "GABARAP", "FIP200", "PIK3C3", "SQSTM1",  
#                           "NBR1", "ATG14", "TECPR1", "BNIP3")
length(macroautophagy_genes)

# add moudule score
DefaultAssay(sample.integrated) <- "RNA"
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(macroautophagy_genes),
  ctrl = length(macroautophagy_genes),
  name = "mac_ATG7",
  assay = "RNA"
)

# add moudule score
DefaultAssay(sample.integrated) <- "RNA"
genes.ifn <- c("IFI35","IFI44L","IFI6","IFIT3",
                        "IRF7","ISG15","MX1","MX2",
                         "OAS1","OAS2")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.ifn),
  ctrl = length(genes.ifn),
  name = "ifn_stim",
  assay = "RNA"
)

# add moudule score
DefaultAssay(sample.integrated) <- "RNA"
genes.HLA<- c("HLA-DQA2","HLA-DQA2")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.HLA),
  ctrl = length(genes.HLA),
  name = "HLA_DQA2",
  assay = "RNA"
)


# add moudule score
DefaultAssay(sample.integrated) <- "RNA"
genes.ATG7<- c("ATG7","ATG7")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.ATG7),
  ctrl = length(genes.ATG7),
  name = "ATG7",
  assay = "RNA"
)


VIA_genes <- readRDS("antigen_module_genes.rds")
length(unique(VIA_genes[[1]]))

unique(VIA_genes[[1]])

VIA.genes <- c("ACTG1", "PCLAF", "ACTB", "GAPDH", "CD38", "GZMH", "PFN1", "CORO1A", "PYCARD", "PSMB8", "TMSB10", "ARPC3", 
  "CLIC1", "CD52", "TMSB4X", "APOBEC3H", "ADA", "SEM1", "MYL6", "GZMA", "DECR1", "ARPC2", "NOP10", "APOBEC3C", 
  "LGALS1", "GZMB", "CFL1", "CLTA", "TWF2", "ISG15", "MXD4", "SUB1", "PSMA4", "HLA-DRA", "PPP1R18", "GSTO1", 
  "APOBEC3G", "ARPC5", "UCP2", "PSMA5", "WDR1", "ARHGDIB", "YWHAE", "CCL5", "SLC25A5", "PSMB4", "ITGB1", 
  "DYNLL1", "ARPC4", "ATP5MC2", "FGFBP2", "GTF3C6", "RGS10", "ATP6V1D", "S1PR4", "EZH2", "ATP5F1A", "GMFG", 
  "OSTF1", "PRF1", "ABRACL", "GLRX", "ACTR3", "SIT1", "TXNDC17", "H3F3A", "CCDC28B", "ATP5MG", "DNAJC15", 
  "IDH2", "SEC11A", "FABP5", "CAPZB", "DBI", "IFI16", "DCTN2", "COTL1", "CNN2", "CHMP2A", "ATP5MC3", "ARPC1B", 
  "NDUFS2", "MT1E", "MRPL28", "OAS2", "SERPINB1", "CARHSP1", "TIMD4", "UBE2L6", "MYL6B", "DCTN3", "IRF4", 
  "ANXA5", "ATP5F1B", "GIMAP4", "PSMB2", "NDUFB3", "NCOA4", "YARS", "CARD16", "COX5A", "FKBP1A", "TCEAL8", 
  "CLDND1", "NDUFB5", "ATP5F1C", "MRPL42", "RPA3", "MT1F", "NUDT5", "POLR2G", "LDHB", "GBP1", "ATP5MF", 
  "COX6B1", "COX6C", "SUCLG1", "TAP1", "CLIC3", "TPM4", "GBP2", "NDUFC2", "CHCHD1", "PPP1CA", "TKT", "SH3BGRL3", 
  "EIF4E2", "LAIR2", "MRPL51", "TRAPPC1", "CSK", "BCL2L11", "LIMD2", "JPT1", "GTF2H5", "NDUFA12", "TROAP", 
  "RTRAF", "ZYX", "RAC2", "CHCHD5", "MPG", "RPS6KA1", "SRSF9", "CEBPD", "GIMAP1", "COPZ1", "TALDO1", "CALM3", 
  "EIF2S2", "HLA-DMA", "POMP", "PTRHD1", "COMMD4", "GGCT", "HAVCR2", "LAP3", "FERMT3", "COPS9", "QARS", 
  "AP2S1", "ZNHIT1", "ICA1", "HADHB", "BLOC1S1", "RABL3", "IFI27L2", "RACK1", "CD27", "FKBP3", "MRPL10", 
  "TXNDC9", "FIBP", "PTTG1", "PSMD8", "CXCR3", "PSME1", "LCP1", "ACAA2", "CHI3L2", "NCKAP1L", "PPP4C", "IGBP1", 
  "PSMB3", "LAMTOR2", "ARHGAP30", "TMEM256", "COPB2", "AC010618.1", "AP1S1", "PFDN4", "NT5C3A", "HMGA1", 
  "MT-CO1", "PSMA2", "NSMCE1", "PARK7")

VIA.genes <- unique(VIA_genes[[1]])

sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(VIA.genes),
  ctrl = length(VIA.genes),
  name = "VIA",
  assay = "RNA"
)


# add moudule score
DefaultAssay(sample.integrated) <- "RNA"
genes.Apoptosis <- unique(c("AGO3", "CBX4", "UBC", "DAPK1", "TNFRSF1B",
                      "MDM4", "PIK3R1", "FOS", "SIVA1", "APPL1", 
                      "RAF1", "JAK2", "HIST1H2AC", "HIST1H2BC", 
                      "HIPK1", "HIST1H2BK", "CDKN2D", "JUN", "CBX6", 
                      "PHC3", "MOAP1", "BAK1", "GSK3B", "RPS27A"))
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.Apoptosis),
  ctrl = length(genes.Apoptosis),
  name = "Apoptosis",
  assay = "RNA"
)



Idents(sample.integrated) <- 'group'
cv_19_part <- subset(sample.integrated,idents=c('cv_19'))
mRNA_part <- subset(sample.integrated,idents=c('HC','Vax','booster'))

sample.integrated$state <- 'mRNA'
sample.integrated$state[sample.integrated$group == 'cv_19'] <- 'cv_19'
sample.integrated$state[sample.integrated$group == 'Vax'] <- 'mRNA'
sample.integrated$state[sample.integrated$group == 'HC'] <- 'mRNA'
sample.integrated$state[sample.integrated$group == 'booster'] <- 'mRNA'
table(sample.integrated$state)

colnames(sample.integrated@meta.data)

write.csv(sample.integrated@meta.data,  "data/1_mRNA_meta_data_df.csv")







source("/share/home/qlab/projects/project_cvv/yyp_results_3/yyp_function.R") #加载到当前的R环境中

unique(sample.integrated@meta.data$ct_merge)

Idents(sample.integrated) <-  "ct_merge"
pDC <- subset(sample.integrated,
             idents = c('pDC'))
unique(pDC$ct_merge)
nrow(pDC@meta.data)

pDC$Condition <- factor(
  pDC$Condition,
  levels = c("A0" ,"A1","A2" ,"B0" ,"B1","B2" ,"C0" ,"C1","C2","C3","acute","convalescent" ))
unique(pDC$Condition)
table(pDC$Condition)

saveRDS(pDC,"data/pDC.rds")

colnames(pDC@meta.data)

idents_name_list <- c("Condition")
plot_list <- c("Apoptosis1","ifn_stim1","HLA_DQA21","VIA1","nCount_RNA","nFeature_RNA"
              )

pdf( "figures/9_pDC_QC_vinplot.pdf",6,3.5)
for (idents_name in idents_name_list) {
  Idents(pDC) <- idents_name
  print(idents_name)
  for (plot_name in plot_list) {
    p <- VlnPlot(pDC, features = c(plot_name), pt.size = 0, raster=FALSE) 
    print(p)
  }
}
dev.off()

genes.Apoptosis <- unique(c("AGO3", "CBX4", "UBC", "DAPK1", "TNFRSF1B",
                      "MDM4", "PIK3R1", "FOS", "SIVA1", "APPL1", 
                      "RAF1", "JAK2", "HIST1H2AC", "HIST1H2BC", 
                      "HIPK1", "HIST1H2BK", "CDKN2D", "JUN", "CBX6", 
                      "PHC3", "MOAP1", "BAK1", "GSK3B", "RPS27A"))
# Define gene list
gene_list <- list(
   genes.Apoptosis
)

DefaultAssay(pDC) <- "RNA"
# Generate dot plot
p <- Dotplot_fun(
  seuratObj = pDC,
  genes = gene_list,
    coord_flip = T,
  group.by = "Condition",
    dot.scale = 2
)

# Customize plot
p <- p +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),
    text = element_text(size = 0),
    panel.grid = element_line(colour = "grey", linetype = "dashed", size = 0.1),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = 'right',
    legend.spacing = unit(0.01, "lines"), # 调整两个legend之间的间距
    legend.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.1, unit = "cm"), # 设置legend的外边距为0
    legend.key.size = unit(0.12, "inch"),
    plot.margin = unit(c(2, 0, 0, 0), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 2.5))

p

ggsave("figures/9_pDC_heatmap_select_Apoptosis_genes.pdf",width = 3,height = 3.6)





