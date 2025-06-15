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
library(clustree)

# 设置路径
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/2_QC_raw_data"
data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results/data"
# setwd()指定路径，setwd("~/project/")

#运行环境设置
source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
set.seed(2)
system("hostname") #计算机主机名
use_python("~/miniconda3/envs/SC/bin/python", required = FALSE)
use_condaenv("SC",conda="~/miniconda3/bin/conda")
scr = import("scrublet")


sample.combined <- readRDS("/share/home/qlab/projects/project_cvv/data/0919_sample_combined_yrs.rds")
sample.combined@meta.data$orig.ident%>% unique()


#计算合并数据的质控数据
sample.combined  <- NormalizeData(object = sample.combined, normalization.method = "LogNormalize",
                                  scale.factor = 10000)
sample.combined[["percent.mt"]] <- PercentageFeatureSet(sample.combined, pattern = "^MT-")

mean(sample.combined@meta.data$percent.mt)

sample.combined <- CellCycleScoring(sample.combined, s.features = cc.genes$s.genes,
                                    g2m.features = cc.genes$g2m.genes, set.ident = TRUE)


#参考基因数据
msigdb <- getBroadSets("/share/home/qlab/projects/project_cvv/yyp_results/data/msigdb_v6.1.xml")
index <- sapply(msigdb, function(gs)
  bcCategory(collectionType(gs))=="c2")
geneset.c2 = msigdb[index]

#干扰素基因
geneset.interferon <- geneset.c2[["BROWNE_INTERFERON_RESPONSIVE_GENES"]] #干扰素响应相关的基因
interferon.genes <- geneset.interferon@geneIds
interferon.genes  <-  interferon.genes[ interferon.genes %in% rownames(sample.combined)]
sample.combined <- AddModuleScore(
  object = sample.combined,
  features = list(interferon.genes),
  ctrl = length(interferon.genes),
  name = 'IFN.Score',
  assay = "RNA")

#32 stress genes
stress.genes <- read_csv("/share/home/qlab/projects/project_cvv/data/stress.genes.csv") %>%
  dplyr::filter(gene %in% rownames(sample.combined)) %>%
  dplyr::pull(gene)

#添加基因表达均值分数,Zscore=（X-μ）/σ
sample.combined <- AddModuleScore( object = sample.combined,
                                   features = list(stress.genes),
                                   ctrl = length(stress.genes),
                                   name = 'Stress.Score',
                                   assay = "RNA")

# hypoxia genes
genes.hypoxia<- c("VEGFA", "SLC2A1", "PGAM1", "ENO1",
                  "LDHA", "TPI1", "P4HA1", "MRPS17",
                  "CDKN3", "ADM", "NDRG1", "TUBB6",
                  "ALDOA", "MIF", "ACOT7")

sample.combined <- AddModuleScore(object = sample.combined,
                                    features = list(genes.hypoxia),
                                    ctrl = length(genes.hypoxia),
                                    name = 'Hypoxia.Score',
                                    #search = FALSE,
                                    assay = "RNA")

plan(multisession, workers=6)
options(future.globals.maxSize = 100 * 1000 * 1024^2)
sample.combined <- ScaleData(sample.combined,
                             vars.to.regress = c("Batch","nFeature_RNA","nCount_RNA",
                                                 "percent.mt","G2M.Score","S.Score"))

# 寻找变量特征
plan(multisession, workers=6)
options(future.globals.maxSize = 100 * 1000 * 1024^2)
sample.combined <- FindVariableFeatures(sample.combined, selection.method = "vst", nfeatures = 2500)
sample.combined <- RunPCA(sample.combined, npcs = 36)
sample.combined <- RunUMAP(sample.combined, dims = 1:30, umap.method = "umap-learn")
sample.combined <- FindNeighbors(sample.combined, dims = 1:30, force.recalc = T)
sample.combined <- FindClusters(sample.combined, resolution = c(0.1, 0.5, 1, 2, 4 ,8,10) )

saveRDS(sample.combined, file.path(results_path, "2_QC_raw_umap.rds"))

1

pdf(file.path(results_path, "PCA_raw.pdf"))
DimHeatmap(object = sample.combined, dims = 1:9, cells = 1000, balanced = TRUE,raster = FALSE)
DimHeatmap(object = sample.combined, dims = 10:18, cells = 1000, balanced = TRUE,raster = FALSE)
DimHeatmap(object = sample.combined, dims = 19:27, cells = 1000, balanced = TRUE,raster = FALSE)
DimHeatmap(object = sample.combined, dims = 28:35, cells = 1000, balanced = TRUE,raster = FALSE)
ElbowPlot(object = sample.combined, ndims = 35, reduction = "pca")
dev.off()

pdf(file.path(results_path, "Clustree_raw.pdf"),width=12,height=10)
clustree(sample.combined)
dev.off()

idents_name_list <- c("RNA_snn_res.0.1","RNA_snn_res.0.5","RNA_snn_res.1","RNA_snn_res.2")
pdf(file.path(results_path, "Umap_raw_DimPlot_0.1_2.pdf"),12,8)
for (idents_name in idents_name_list) {
  Idents(sample.combined) <- idents_name
  p <- DimPlot(object = sample.combined, reduction = 'umap',label = T,raster = FALSE)
  print(p)
}
dev.off()

idents_name_list <- c("RNA_snn_res.4","RNA_snn_res.8","RNA_snn_res.10")
pdf(file.path(results_path, "Umap_raw_DimPlot_4_10.pdf"),20,8)
for (idents_name in idents_name_list) {
  Idents(sample.combined) <- idents_name
  p <- DimPlot(object = sample.combined, reduction = 'umap',label = T,raster = FALSE)
  print(p)
}
dev.off()

idents_name_list <- c("Donor","Condition","Sample","Batch","HTO_classification.global",
                      "scrublet_callrate", "DF.classifications")
pdf(file.path(results_path, "Umap_raw_DimPlot.pdf"),12,8)
for (idents_name in idents_name_list) {
  Idents(sample.combined) <- idents_name
  p <- DimPlot(object = sample.combined, reduction = 'umap',label = F,raster = FALSE)
  print(p)
}
dev.off()

idents_name_list <- c("RNA_snn_res.2","Donor","Condition","Sample","Batch")
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score","Stress.Score1",
               "IFN.Score1","scrublet_score","DF.pANN")
pdf(file.path(results_path, "QC_raw_vinplot.pdf"),20,10)
for (idents_name in idents_name_list) {
  Idents(sample.combined) <- idents_name
  print(idents_name)
  for (plot_name in plot_list) {
    p <- VlnPlot(sample.combined, features = c(plot_name), pt.size = 0, raster=FALSE)
    print(p)
  }
}
dev.off()

idents_name_list <- c("Donor","Condition","Sample","Batch","HTO_classification.global")
pdf(file.path(results_path, "QC_raw_scatter_all.pdf"))
for (idents_name in idents_name_list) {
  Idents(sample.combined) <- idents_name
  #pdf(file.path(results_path, paste0("1_QC_raw_scatter_", idents_name, ".pdf")))
  p <- FeatureScatter(object = sample.combined, feature1 = 'nFeature_RNA', feature2 = 'nCount_RNA',
                      pt.size = 0.1,raster=FALSE) + geom_smooth()
  print(p)
  p <- FeatureScatter(object = sample.combined, feature1 = 'nFeature_RNA', feature2 = 'percent.mt',
                      pt.size = 0.1,raster=FALSE) + geom_smooth() +xlim(0,1000)
  print(p)
  p <- FeatureScatter(object = sample.combined, feature1 = 'nCount_RNA', feature2 = 'percent.mt',
                      pt.size = 0.1,raster=FALSE) +geom_smooth() +xlim(0,2000)
  print(p)
  #dev.off()
}
dev.off()


plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score","Stress.Score1",
               "IFN.Score1","scrublet_score","DF.pANN")
pdf(file.path(results_path, "Umap_raw_FeaturePlot.pdf"))
for (plot_name in plot_list) {
  p <- FeaturePlot(object = sample.combined,features = plot_name,cols = c("grey", "blue"), pt.size = 0.1,
                   reduction = "umap",raster=FALSE)
  print(p)
}
dev.off()

Idents(sample.combined) <- "RNA_snn_res.2"
pdf(file.path(results_path, "RNA_FeaturePlot_Genes.pdf"), width = 20, height =20)
FeaturePlot(object = sample.combined,
            features = c("CD3D","CD4","CD8A","SLC4A10"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
FeaturePlot(object = sample.combined ,
            features = c("CCR7","GZMA","GZMK","GZMB"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
FeaturePlot(object = sample.combined,
            features = c("JCHAIN","CD1C", "PF4","FOXP3"),
            cols = c("grey", "blue"), pt.size = 0.1,
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
FeaturePlot(object = sample.combined,
            features = c("TRDC","TRGV9", "TRGV10","TRDV1"),
            cols = c("grey", "blue"), pt.size = 0.1,
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
FeaturePlot(object = sample.combined,
            features = c("CD79A","TCL1A","IGHA1","MZB1"),
            cols = c("grey", "blue"), pt.size = 0.1,
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
FeaturePlot(object = sample.combined,
            features = c("CD27","IGHG1","IGHM","IGHD"),
            cols = c("grey", "blue"), pt.size = 0.1,
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
FeaturePlot(object = sample.combined,
            features = c("CD14","FCGR3A","FCGR3B","LILRA4"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
FeaturePlot(object = sample.combined,
            features = c("CD68","CD34","CFP","STAB1","MKI67"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

sample.combined <- FindClusters(sample.combined, resolution = c(1.2) )

idents_name_list <- c("RNA_snn_res.1.2")
pdf(file.path(results_path, "Umap_raw_DimPlot_1.2.pdf"),12,8)
for (idents_name in idents_name_list) {
  Idents(sample.combined) <- idents_name
  p <- DimPlot(object = sample.combined, reduction = 'umap',label = T,raster = FALSE)
  print(p)
}
dev.off()

idents_name_list <- c("RNA_snn_res.1.2")
pdf(file.path(results_path, "QC_raw_vinplot_1.2.pdf"),20,10)
for (idents_name in idents_name_list) {
  Idents(sample.combined) <- idents_name
  print(idents_name)
  for (plot_name in plot_list) {
    p <- VlnPlot(sample.combined, features = c(plot_name), pt.size = 0, raster=FALSE)
    print(p)
  }
}
dev.off()



Idents(sample.combined) <- "RNA_snn_res.2"
pdf(file.path(results_path, "RNA_FeaturePlot_Genes_HBB.pdf"), width = 20, height =20)
FeaturePlot(object = sample.combined,
            features = c("HBA1", "HBA2", "HBB", "HBD", "HBE1",
                         "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",
            raster=FALSE,label= TRUE,label.size = 3)

dev.off()


