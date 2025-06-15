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
#library(clustree)

# 设置路径
# 设置路径
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/2_Mono/1_CCA"
data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/4_Cell_type_level_1_2"
# setwd()指定路径，setwd("~/project/")

#运行环境设置
source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
set.seed(2)
system("hostname") #计算机主机名
use_python("~/miniconda3/envs/SC/bin/python", required = FALSE)
use_condaenv("SC",conda="~/miniconda3/bin/conda")
scr = import("scrublet")

targetcell <- readRDS(file.path(data_save_path, "4_integrated_level_1_2.rds"))

unique(targetcell$ct_level_1)

# ---- subset cell
Idents(targetcell) <- "ct_level_1"
targetcell <- subset(targetcell, idents = c("Monocyte"))

targetcell@meta.data %>%
  dplyr::group_by(ct_level_2) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)
targetcell@meta.data %>%
  dplyr::group_by(Batch) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num) 

#msigdb <- getBroadSets("/share/home/qlab/projects/project_cvv/yyp_results/data/msigdb_v6.1.xml")
msigdb <- getBroadSets("/share/reference/sigdb/msigdb_v7.4.xml")
index <- sapply(msigdb, function(gs) {
  bcCategory(collectionType(gs)) == "c2"
})
geneset.c2 <- msigdb[index]
geneset.interferon <- geneset.c2[["BROWNE_INTERFERON_RESPONSIVE_GENES"]]
interferon.genes <- geneset.interferon@geneIds
interferon.genes <- interferon.genes[interferon.genes %in% rownames(targetcell)]
cellcycle.genes <- geneset.c2[["KEGG_CELL_CYCLE"]]@geneIds
    
genes.hypoxia <- c(
  "VEGFA", "SLC2A1", "PGAM1", "ENO1",
  "LDHA", "TPI1", "P4HA1", "MRPS17",
  "CDKN3", "ADM", "NDRG1", "TUBB6",
  "ALDOA", "MIF", "ACOT7"
)    
stress.genes <- read_csv("/share/home/qlab/projects/project_cvv/data/stress.genes.csv") %>%
  dplyr::filter(gene %in% rownames(targetcell)) %>%
  dplyr::pull(gene)
    
# HG genes
hgGenes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
# IG genes
igGenes <- c(
  "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "JCHAIN",
  "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"
)

IG.gene1 <- grep(pattern = "^IGLV", x = rownames(x = targetcell), value = TRUE)
IG.gene2 <- grep(pattern = "^IGKV", x = rownames(x = targetcell), value = TRUE)
IG.gene3 <- grep(pattern = "^IGHV", x = rownames(x = targetcell), value = TRUE)
IG.genes <- union(IG.gene1, c(IG.gene2, IG.gene3))
# TRV genes
TRDV.genes <- grep(pattern = "^TRDV", x = rownames(x = targetcell), value = TRUE)
TRAV.genes <- grep(pattern = "^TRAV", x = rownames(x = targetcell), value = TRUE)
TRBV.genes <- grep(pattern = "^TRBV", x = rownames(x = targetcell), value = TRUE)
TRGV.genes <- grep(pattern = "^TRGV", x = rownames(x = targetcell), value = TRUE)
TRV.genes <- union(TRDV.genes, c(TRAV.genes, TRBV.genes, TRGV.genes))
# RPS & RPL
RPS.genes <- grep(pattern = "^RPS", x = rownames(x = targetcell), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(x = targetcell), value = TRUE)
# MT genes
MT.genes <- grep(pattern = "^MT-", x = rownames(x = targetcell), value = TRUE)
# HIST genes
HIST.genes <- grep(pattern = "^HIST", x = rownames(x = targetcell), value = TRUE)

DefaultAssay(targetcell) <- "RNA"
sample.list <- SplitObject(targetcell, split.by = "Batch")
sample.list <- purrr::map(sample.list, function(x) {
  x <- NormalizeData(x, verbose = TRUE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features.select <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)

# ---- discard the genes not needed
genes.cc <- extract_cellcycle(targetcell, features.select, cores = 10, assay = "RNA", cutoff = 0.05)
genes.stress <- extract_stress(targetcell, features.select, cores = 10, assay = "RNA", cutoff = 0.1)
genes.ifn <- extract_ifn(targetcell, features.select, cores = 10, assay = "RNA", cutoff = 0.1)
features.filtered <- setdiff(features.select, c(
  MT.genes, RPL.genes, RPS.genes, hgGenes, HIST.genes,IG.genes, TRV.genes, stress.genes,cellcycle.genes,
  genes.cc[1:500], genes.stress[1:10], genes.ifn[1:40], genes.hypoxia
))

# ---- scale and PCA for each sample
sample.list <- purrr::map(sample.list, function(x) {
  x <- ScaleData(x,
    vars.to.regress = c("nFeature_RNA", "percent.mt", "G2M.Score", "S.Score"),
    features = features.filtered
  )
  x <- RunPCA(x, features = features.filtered, npcs = 30, verbose = FALSE)
})

# ---- integrate samples
# load("sample_list.rda")
plan("multicore", workers = 4)
plan()
options(future.globals.maxSize = 60 * 1024^3)
gc()
plan()
sample.anchors <- FindIntegrationAnchors(
  sample.list,
  reduction = "cca", #reduction ="rpca",
  k.anchor = 5,
  k.filter = 200, #default:k.weight = 200
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
  k.weight = 100, #default:k.weight = 100
  verbose = TRUE
)

# ---- scale expression
sample.integrated <- FindVariableFeatures(sample.integrated, selection.method = "vst", nfeatures = 3500, verbose = F)
var.features <- VariableFeatures(sample.integrated)
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- ScaleData(
  sample.integrated,
  features = var.features
)

# ---- variable genes
genes.stress <- extract_stress(sample.integrated, var.features, cores = 10, assay = "integrated", cutoff = 0.1)
genes.cc <- extract_cellcycle(sample.integrated, var.features, cores = 10, assay = "integrated", cutoff = 0.01)
genes.ifn <- extract_ifn(sample.integrated, var.features, cores = 10, assay = "integrated", cutoff = 0.01)
var.features.filtered <- var.features[!var.features %in% c(
  genes.stress[1:10], genes.cc[1:200],genes.hypoxia, genes.ifn[1:40],
  MT.genes, RPL.genes, RPS.genes, hgGenes, HIST.genes,IG.genes, TRV.genes, stress.genes,cellcycle.genes
)]

# ---- PCA and clustering
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- RunPCA(sample.integrated, npcs = 30, features = var.features.filtered)
sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn", reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(0.1, 2, 0.1))

Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "CCA_DimPlot_Donor.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Donor',
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "CCA_DimPlot_Batch.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Batch',
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "CCA_DimPlot_Sample.pdf"),25,25)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Sample',
             ncol = 6,label = T,raster = FALSE)
print(p)
dev.off()

allcolors <- c(
    "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
    "#968175", "#51A3CCFF", "#f18800","#85B22CFF","#B22C2CFF",
    "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
            "#E57E7EFF", "#BFB2FFFF","#e20612", "#ffd401", 
             "#E5B17EFF", "#942d8d","#573333FF", "#2e409a",
    "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",      
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",        
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",     
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887"
)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.1","integrated_snn_res.0.5","ct_level_2","ct_level_1")
pdf(file.path(results_path, "DimPlot_res.1.pdf"),6,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.0.1","integrated_snn_res.0.2",
                      "integrated_snn_res.0.3","integrated_snn_res.0.4",
                      "integrated_snn_res.0.5","integrated_snn_res.0.6",
                      "integrated_snn_res.0.7","integrated_snn_res.0.8",
                      "integrated_snn_res.0.9","integrated_snn_res.1",
                      "integrated_snn_res.1.5","integrated_snn_res.2",
                     "integrated_snn_res.2.5","integrated_snn_res.3",
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
Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_Mono_cell.pdf"), width = 40, height=40)
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

Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_others.pdf"), width = 40, height =50)
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
pdf(file.path(results_path, "CCA_QC_vinplot.pdf"))
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

sample.integrated <- readRDS("CCA_Monocyte_Batch.rds")

#----细胞注释---------
Idents(sample.integrated)  <- "ct_level_2"
sample.integrated@meta.data$ct_level_3  <- Idents(sample.integrated)
sample.integrated@meta.data$ct_level_3  <- as.character(sample.integrated@meta.data$ct_level_3)

sample.integrated$ct_level_3[sample.integrated$integrated_snn_res.0.5 %in% c(10)] <- "Mono_ATG7+"
sample.integrated$ct_level_3[sample.integrated$integrated_snn_res.0.5 %in% c(2,14)] <- "Mono_CD16"
sample.integrated$ct_level_3[sample.integrated$integrated_snn_res.0.5 %in% c(1,3,5,0,11,9,6)] <- "Mono_CD14"
sample.integrated$ct_level_3[sample.integrated$integrated_snn_res.3 == 34] <- "Mono_CD14_CD16"
sample.integrated$ct_level_3[sample.integrated$integrated_snn_res.0.5 %in% c(8,13,7,4,16,12)] <- "Doublet_Mono"
sample.integrated$ct_level_3[sample.integrated$integrated_snn_res.0.5 %in% c(15)] <- "Doublet_Mono_Platelet"

unique(sample.integrated@meta.data$ct_level_3)

allcolors <- c(
    "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
    "#968175", "#51A3CCFF", "#f18800","#85B22CFF","#B22C2CFF",
    "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
            "#E57E7EFF", "#BFB2FFFF","#e20612", "#ffd401", 
             "#E5B17EFF", "#942d8d","#573333FF", "#2e409a",
    "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",      
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",        
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",     
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887"
)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("ct_level_3","ct_level_2","ct_level_1")
pdf(file.path(results_path, "ct_level_3_DimPlot.pdf"),6,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

saveRDS(sample.integrated, "CCA_Monocyte_Batch.rds")

unique(sample.integrated@meta.data$ct_level_3)

sample.integrated@meta.data %>%
  dplyr::group_by(ct_level_3) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)

# ---- subset cell
Idents(sample.integrated) <- "ct_level_3"
targetcell_D<- subset(sample.integrated, idents = c("Doublet_Mono_Platelet",
                                                    "Doublet_Mono"))
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/0_merge_data"
saveRDS(targetcell_D, file = paste0(results_path,"/2_CCA_Monocyte_batch_Doublet.rds"))

# Reumap ------------------------------------------------------------------

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
#library(clustree)

# 设置路径
# 设置路径
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/2_Mono/2_Reumap"
data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/2_Mono/1_CCA"
# setwd()指定路径，setwd("~/project/")

#运行环境设置
source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
# set.seed(2)
# system("hostname") #计算机主机名
# use_python("~/miniconda3/envs/SC/bin/python", required = FALSE)
# use_condaenv("SC",conda="~/miniconda3/bin/conda")
# scr = import("scrublet")

Monocyte <- readRDS(file.path(data_save_path, "CCA_Monocyte_Batch.rds"))

Monocyte@meta.data %>%
  dplyr::group_by(ct_level_3) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)

targetcell<- Monocyte
unique(targetcell$ct_level_1)
unique(targetcell$ct_level_2)
unique(targetcell$ct_level_3)

# ---- subset cell
Idents(targetcell) <- "ct_level_3"
targetcell<- subset(targetcell, idents = c("Mono_CD14","Mono_CD16","Mono_ATG7+","Mono_CD14_CD16"))

unique(targetcell$ct_level_1)
unique(targetcell$ct_level_2)
unique(targetcell$ct_level_3)

targetcell@meta.data %>%
  dplyr::group_by(ct_level_3) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)

sample.integrated <- RunUMAP(targetcell, umap.method = "umap-learn", reduction = "pca", dims = 1:20)

#msigdb <- getBroadSets("/share/home/qlab/projects/project_cvv/yyp_results/data/msigdb_v6.1.xml")
msigdb <- getBroadSets("/share/reference/sigdb/msigdb_v7.4.xml")
index <- sapply(msigdb, function(gs) {
  bcCategory(collectionType(gs)) == "c2"
})
geneset.c2 <- msigdb[index]
geneset.interferon <- geneset.c2[["BROWNE_INTERFERON_RESPONSIVE_GENES"]]
interferon.genes <- geneset.interferon@geneIds
interferon.genes <- interferon.genes[interferon.genes %in% rownames(targetcell)]
cellcycle.genes <- geneset.c2[["KEGG_CELL_CYCLE"]]@geneIds

genes.hypoxia <- c(
  "VEGFA", "SLC2A1", "PGAM1", "ENO1",
  "LDHA", "TPI1", "P4HA1", "MRPS17",
  "CDKN3", "ADM", "NDRG1", "TUBB6",
  "ALDOA", "MIF", "ACOT7"
)    
stress.genes <- read_csv("/share/home/qlab/projects/project_cvv/data/stress.genes.csv") %>%
  dplyr::filter(gene %in% rownames(targetcell)) %>%
  dplyr::pull(gene)

# HG genes
hgGenes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
# IG genes
igGenes <- c(
  "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "JCHAIN",
  "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"
)

IG.gene1 <- grep(pattern = "^IGLV", x = rownames(x = targetcell), value = TRUE)
IG.gene2 <- grep(pattern = "^IGKV", x = rownames(x = targetcell), value = TRUE)
IG.gene3 <- grep(pattern = "^IGHV", x = rownames(x = targetcell), value = TRUE)
IG.genes <- union(IG.gene1, c(IG.gene2, IG.gene3))
# TRV genes
TRDV.genes <- grep(pattern = "^TRDV", x = rownames(x = targetcell), value = TRUE)
TRAV.genes <- grep(pattern = "^TRAV", x = rownames(x = targetcell), value = TRUE)
TRBV.genes <- grep(pattern = "^TRBV", x = rownames(x = targetcell), value = TRUE)
TRGV.genes <- grep(pattern = "^TRGV", x = rownames(x = targetcell), value = TRUE)
TRV.genes <- union(TRDV.genes, c(TRAV.genes, TRBV.genes, TRGV.genes))
# RPS & RPL
RPS.genes <- grep(pattern = "^RPS", x = rownames(x = targetcell), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(x = targetcell), value = TRUE)
# MT genes
MT.genes <- grep(pattern = "^MT-", x = rownames(x = targetcell), value = TRUE)
# HIST genes
HIST.genes <- grep(pattern = "^HIST", x = rownames(x = targetcell), value = TRUE)

# ---- scale expression
plan(multisession, workers=5)
options(future.globals.maxSize = 100 * 1000 * 1024^2)
sample.integrated <- ScaleData(targetcell,
                               vars.to.regress = c("Batch","nFeature_RNA","nCount_RNA",
                                                   "percent.mt","G2M.Score","S.Score"))

sample.integrated <- FindVariableFeatures(sample.integrated, selection.method = "vst", nfeatures = 3500)
var.features <- VariableFeatures(sample.integrated)
var.features.filtered <- var.features[!var.features %in% c(
  MT.genes, RPL.genes, RPS.genes, hgGenes, HIST.genes,IG.genes, TRV.genes, stress.genes,cellcycle.genes)]

# ---- PCA and clustering
sample.integrated <- RunPCA(sample.integrated, npcs = 30, features= var.features.filtered)
sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn", reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(0.1, 2, 0.1))

Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "DimPlot_Donor.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Donor',
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "DimPlot_Batch.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Batch',
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "DimPlot_Sample.pdf"),25,25)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Sample',
             ncol = 6,label = T,raster = FALSE)
print(p)
dev.off()

allcolors <- c(
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
  "#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF", "#BFB2FFFF",
  "#e20612", "#ffd401","00b0eb", "#f18800",
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#968175","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",      
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
                      "integrated_snn_res.2.5","integrated_snn_res.3",
                      "Donor","Condition","Sample","Batch","HTO_classification.global",
                      "scrublet_callrate","DF.classifications")
pdf(file.path(results_path, "DimPlot_all.pdf"),5,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.0.5","integrated_snn_res.1","integrated_snn_res.2")
pdf(file.path(results_path, "DimPlot_res.0.5_1_2.pdf"),4,4)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("ct_level_3","ct_level_1","ct_level_2")
pdf(file.path(results_path, "DimPlot_ct_level_3.pdf"),5,4)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_T_cell.pdf"), width = 20, height=30)
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

My study explores the differences in immune responses during a three-dose COVID vaccination process.
Our data originates from 12 healthy individuals who have received between two and three vaccine doses. 
We divided the samples into nine conditions (from A0 to C2:A0,A1,A2,B0,B1,B2,C0,C1,C2) according to the sampling time. 
Using several methods, we obtained scRNA data, TCR and BCR data, antibody titer data, and metabolomics data for our samples.


Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_others.pdf"), width = 20, height =32)
FeaturePlot(object = sample.integrated,
            features = c("CD3D","CD4","CD8A","SLC4A10",
                         "CCR7","GZMA","GZMK","GZMB",
                         "JCHAIN","CD1C", "PF4","PPBP","FOXP3",
                         "TRDC","TRGV9", "TRGV10","TRDV1",
                         "CD79A","TCL1A","IGHA1","MZB1",
                         "CD27","IGHG1","IGHM","IGHD",
                         "CD14","FCGR3A","FCGR3B","LILRA4",
                         "CD68","CD34","CFP","STAB1","MKI67"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",
            raster=FALSE,label= TRUE,label.size = 3)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_Mono_cell.pdf"), width = 40, height=40)
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

idents_name_list <- c("integrated_snn_res.1","Donor","Condition","Sample","Batch")
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score")
pdf(file.path(results_path, "QC_vinplot.pdf"))
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

sample.integrated <- readRDS("Trajectory_Mono.rds")

#----细胞注释---------
Idents(sample.integrated)  <- "ct_level_3"
sample.integrated@meta.data$ct_level_4  <- Idents(sample.integrated)
sample.integrated@meta.data$ct_level_4  <- as.character(sample.integrated@meta.data$ct_level_4)

sample.integrated$ct_level_4[sample.integrated$integrated_snn_res.0.5 == 14] <- "Mono_CD16_ATG7+"
sample.integrated$ct_level_4[sample.integrated$ct_level_4 == "Mono_ATG7+"] <- "Mono_CD14_ATG7+"
# sample.integrated$ct_level_4[sample.integrated$integrated_snn_res.0.1  == 3] <- "Mono_ATG7+"
# sample.integrated$ct_level_4[sample.integrated$integrated_snn_res.0.1  == 5] <- "Mono_IFI44L"
# sample.integrated$ct_level_4[sample.integrated$integrated_snn_res.0.1 == 1] <- "Mono_CD16"
# sample.integrated$ct_level_4[sample.integrated$integrated_snn_res.0.1 == 6] <- "D_Mono_platelet"
# sample.integrated$ct_level_4[sample.integrated$integrated_snn_res.0.1 == 4] <- "D_T_Mono"
# sample.integrated$ct_level_4[sample.integrated$integrated_snn_res.1 == 7] <- "Mono_CD14_CD16"

sample.integrated@meta.data %>%
  dplyr::group_by(ct_level_4) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)

allcolors <- c(
  "#53A85F", "#F1BB72", "#F3B1A0","#E5D2DD", "#D6E7A3", "#57C3F3", "#476D87"
)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("ct_level_4","ct_level_3","ct_level_2","ct_level_1")
pdf(file.path(results_path, "ct_level_4_DimPlot.pdf"),8,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.1,label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

#saveRDS(targetcell_use, file = "Reumap_Mono.rds")

sample.integrated@meta.data %>%
  dplyr::group_by(ct_level_4) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)

library(SingleCellExperiment)
library(slingshot, quietly = FALSE)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(tradeSeq)
library(magrittr)
library(Seurat)
library("tidyverse")
library(slingshot)

targetcell_use <- sample.integrated

Idents(targetcell_use) <- "ct_level_4"
sds <- slingshot(
  Embeddings(targetcell_use, "umap"),
  clusterLabels = targetcell_use$ct_level_4,
  start.clus = "Mono_CD14",
  end.clus = c("Mono_CD16","Mono_CD14_ATG7+","Mono_CD16_ATG7+"),
  use.median = T,#FALSE,
  #omega = T,#FALSE,
  #omega_scale = 1.5,
  stretch = 0
)
lineages <- slingLineages(sds)
lineages

curves <- slingCurves(sds, as.df = TRUE)
curved <- rbind(curves[curves$Lineage == 1, ], curves[curves$Lineage == 2, ], curves[curves$Lineage == 3, ])
curved$Lineage <- paste0("Lineage", curved$Lineage)
table(curved$Lineage)

targetcell_use$Pseudotime1 <- slingPseudotime(sds)[, 1]
targetcell_use$Pseudotime1 <- (targetcell_use$Pseudotime1 - min(targetcell_use$Pseudotime1, na.rm = T)) / (max(targetcell_use$Pseudotime1, na.rm = T) - min(targetcell_use$Pseudotime1, na.rm = T)) * 100
targetcell_use$Pseudotime2 <- slingPseudotime(sds)[, 2]
targetcell_use$Pseudotime2 <- (targetcell_use$Pseudotime2 - min(targetcell_use$Pseudotime2, na.rm = T)) / (max(targetcell_use$Pseudotime2, na.rm = T) - min(targetcell_use$Pseudotime2, na.rm = T)) * 100
targetcell_use$Pseudotime3 <- slingPseudotime(sds)[, 3]
targetcell_use$Pseudotime3 <- (targetcell_use$Pseudotime3 - min(targetcell_use$Pseudotime3, na.rm = T)) / (max(targetcell_use$Pseudotime3, na.rm = T) - min(targetcell_use$Pseudotime3, na.rm = T)) * 100
# # # targetcell_use$Pseudotime4 <- slingPseudotime(sds)[, 4]
# # targetcell_use$Pseudotime4 <- (targetcell_use$Pseudotime4 - min(targetcell_use$Pseudotime4, na.rm = T)) / (max(targetcell_use$Pseudotime4, na.rm = T) - min(targetcell_use$Pseudotime4, na.rm = T)) * 100
# targetcell_use$Pseudotime5 <- slingPseudotime(sds)[, 5]
# targetcell_use$Pseudotime5 <- (targetcell_use$Pseudotime5 - min(targetcell_use$Pseudotime5, na.rm = T)) / (max(targetcell_use$Pseudotime5, na.rm = T) - min(targetcell_use$Pseudotime5, na.rm = T)) * 100

targetcell_use$Pseudotime <- rowMeans(targetcell_use@meta.data[, c("Pseudotime1", "Pseudotime2", "Pseudotime3")], na.rm = T)
head(targetcell_use@meta.data[, c("Pseudotime1","Pseudotime2","Pseudotime3",  "Pseudotime")])

data <- FetchData(targetcell_use,
                  vars = c("UMAP_1", "UMAP_2", "Pseudotime", "Pseudotime1","Pseudotime2","Pseudotime3")
)

head(data)

ncols <- c(
  "#53A85F", "#F1BB72", "#F3B1A0","#E5D2DD", "#D6E7A3", "#57C3F3", "#476D87"
)
pdf(file.path(results_path, "slingshot_Mono_line3.pdf"), width = 8.5, height = 6)
DimPlot(object = targetcell_use, reduction = "umap", label = F, group.by = "ct_level_4", pt.size = 0.5) +
  geom_path(aes_string("UMAP_1", "UMAP_2", linetype = "Lineage"), curved, size = 1, arrow = arrow(length = unit(0.15, "inches"))) +
  scale_fill_manual(values = ncols) +
  scale_color_manual(values = ncols) #+
#NoLegend()
FeaturePlot(targetcell_use, feature = "Pseudotime") + scale_color_distiller(palette = "YlOrRd", na.value = "grey90") #+ NoLegend()
dev.off()

saveRDS(targetcell_use, file = "Trajectory_Mono.rds")

library("tidyverse")
library("ggpubr")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(RColorBrewer)
display.brewer.all()

targetcell_use <- readRDS("Trajectory_Mono.rds")

meta_data.df <- targetcell_use@meta.data

#colpalette <- colorRampPalette(brewer.pal(8, "Set3"))(40)  
colpalette <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
  "#968175"
)
pdf("1_Line_Mono_ct_level_4.pdf", height = 10, width = 8)
#ct_level_4
meta_data.df <- meta_data.df %>%
  group_by(Condition) %>%
  mutate(Count_condition = n()) %>%
  ungroup%>% 
  group_by(Condition, ct_level_4) %>%
  mutate(Count_condition_celltype = n(),Celltype_prop = Count_condition_celltype/Count_condition) 

ggplot(meta_data.df, aes(x = Condition, y = Celltype_prop, group = 1, color=factor(ct_level_4))) +
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values = colpalette) +
  facet_wrap(~ ct_level_4, scales = "free_y", ncol = 1) +
  labs(y = "Proportion", x = "Condition") +
  theme_minimal(base_size = 22) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  ) 
dev.off()



library(patchwork)
library(tidyverse)

n_condition =  c("#57C3F3", "#44b1f0", "#00b0eb","#E5B17EFF","#F4A582", "#f7aa5d", "#cc340c","#B22C2CFF","#931635")

n_condition =  c("#90C3F1","#00A0FA","#5072B2","#F0E442", "#F4A582", "#f18800","#E63863", "#bc340c","#831635")

pdf(file.path(results_path, "2_Mono_celltype_density_2.pdf"), width = 8, height = 4)
# lineage 1
lineage_plot <- ggplot(targetcell_use@meta.data[!is.na(targetcell_use$Pseudotime1), ]) +
  geom_density(aes(x = Pseudotime1, color = Condition)) +
  # scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(colour = "Condition") +
  scale_color_manual(values = n_condition) +
  xlab("Lineage1:Mono_CD14 -> Mono_CD14_CD16 -> Mono_CD16")
celltype_plot <- ggplot(targetcell_use@meta.data[!is.na(targetcell_use$Pseudotime1), ]) +
  geom_jitter(aes(x = Pseudotime1, y = 0.00025, color = ct_level_4), height = 0.0015, size = 0.01) +
  scale_color_manual(values = ncols) +
  theme_bw() +
  theme(
    legend.position = "None",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )
(celltype_plot / lineage_plot) + plot_layout(heights = c(0.15, 2))


# lineage 2
lineage_plot <- ggplot(targetcell_use@meta.data[!is.na(targetcell_use$Pseudotime2), ]) +
  geom_density(aes(x = Pseudotime2, color = Condition)) +
  # scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(colour = "Condition") +
  scale_color_manual(values = n_condition) +
  xlab("Lineage2:Mono_CD14 -> Mono_CD14_CD16 -> Mono_CD16_ATG7+")
celltype_plot <- ggplot(targetcell_use@meta.data[!is.na(targetcell_use$Pseudotime2), ]) +
  geom_jitter(aes(x = Pseudotime2, y = 0.00025, color = ct_level_4), height = 0.0015, size = 0.01) +
  scale_color_manual(values = ncols) +
  theme_bw() +
  theme(
    legend.position = "None",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )
(celltype_plot / lineage_plot) + plot_layout(heights = c(0.15, 2))


# lineage 3
lineage_plot <- ggplot(targetcell_use@meta.data[!is.na(targetcell_use$Pseudotime3), ]) +
  geom_density(aes(x = Pseudotime3, color = Condition)) +
  # scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(colour = "Condition") +
  scale_color_manual(values = n_condition) +
  xlab("Lineage3:Mono_CD14 -> Mono_CD14_ATG7+")
celltype_plot <- ggplot(targetcell_use@meta.data[!is.na(targetcell_use$Pseudotime3), ]) +
  geom_jitter(aes(x = Pseudotime3, y = 0.00025, color = ct_level_4), height = 0.0015, size = 0.01) +
  scale_color_manual(values = ncols) +
  theme_bw() +
  theme(
    legend.position = "None",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )
(celltype_plot / lineage_plot) + plot_layout(heights = c(0.15, 2))
dev.off()

gene <-  c("CD14","FCGR3A","ATG7", #Mono
           "CST3", "LYZ","CD68", "CD163",  #Mono
           "SELL","S100A8","S100A12","CLEC12A","CXCL14",
           "CX3CR1","KLF3","FLI1",#"IFTM1","FAM65B",
           "CD300E","CEBPD","ZFP36L2",
           "FCN1","APOBEC3A","THBS1",
           "MARCO","CD33",
           "PRF1", "GNLY", "KLRC4","KLRK1",
           "CD3D","CD4","CD8A"#T_cell
)

pdf(file.path(results_path, "2_Mono_feature.pdf"), width = 25, height = 35)
Idents(targetcell_use) <- "ct_level_4"
DefaultAssay(targetcell_use) <- "RNA"
FeaturePlot(targetcell_use, features = gene, label = T, ncol = 4)
DefaultAssay(targetcell_use) <- "integrated"
dev.off()

matrix <- GetAssayData(targetcell_use, assay = "RNA")
expr_matrix <- as.data.frame(t(rbind(
  matrix[gene, , drop = FALSE],
  Pseudotime1 = targetcell_use$Pseudotime1,
  Pseudotime2 = targetcell_use$Pseudotime2,
  Pseudotime3 = targetcell_use$Pseudotime3
)))
table(rownames(expr_matrix) %in% colnames(targetcell_use))
expr_matrix$Donor <- targetcell_use$Donor
expr_matrix$Condition <- targetcell_use$Condition
expr_matrix$ct_level_4 <- targetcell_use$ct_level_4
expr_matrix$S.Score <- targetcell_use$S.Score
expr_matrix$G2M.Score <- targetcell_use$G2M.Score
expr_matrix$IFN.Score1 <- targetcell_use$IFN.Score1
expr_matrix$Stress.Score1 <- targetcell_use$Stress.Score1
head(expr_matrix)

df <- pivot_longer(expr_matrix, gene, names_to = "feature", values_to = "expr")
head(df)

# ---- plot
pdf(file.path(results_path, "2_Mono_Cell_Cycle_trend_gam.pdf"), width = 8, height = 5)
# ---- plot S.Score
ggplot(df) +
  geom_smooth(aes(x = Pseudotime1, y = S.Score, color = "#e20612"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime2, y = S.Score, color = "#ffd401"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime3, y = S.Score, color = "#00b0eb"), method = "gam", se = F) + #gam
  xlab("PseudoTime") +
  ylab("S.Score") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_identity(
    name = "Lineages",
    labels = c("Lineage1: Mono_CD16","Lineage2: Mono_CD16_ATG7+","Lineage3: Mono_CD14_ATG7+"),
    guide = "legend"
  ) 

ggplot(df) +
  geom_smooth(aes(x = Pseudotime1, y = S.Score, color = "#e20612"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime2, y = S.Score, color = "#ffd401"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime3, y = S.Score, color = "#00b0eb"), method = "gam", se = F) + #gam
  xlab("PseudoTime") +
  ylab("S.Score") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_identity(
    name = "Lineages",
    labels = c("Lineage1: Mono_CD16","Lineage2: Mono_CD16_ATG7+","Lineage3: Mono_CD14_ATG7+"),
    guide = "legend"
  ) +
  facet_wrap(~Condition, nrow = 3, scales = "free")

# ---- plot G2M.Score
ggplot(df) +
  geom_smooth(aes(x = Pseudotime1, y = G2M.Score, color = "#e20612"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime2, y = G2M.Score, color = "#ffd401"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime3, y = G2M.Score, color = "#00b0eb"), method = "gam", se = F) + #gam
  xlab("PseudoTime") +
  ylab("G2M.Score") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_identity(
    name = "Lineages",
    labels = c("Lineage1: Mono_CD16","Lineage2: Mono_CD16_ATG7+","Lineage3: Mono_CD14_ATG7+"),
    guide = "legend"
  ) 

ggplot(df) +
  geom_smooth(aes(x = Pseudotime1, y = G2M.Score, color = "#e20612"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime2, y = G2M.Score, color = "#ffd401"), method = "gam", se = F) + #gam
  geom_smooth(aes(x = Pseudotime3, y = G2M.Score, color = "#00b0eb"), method = "gam", se = F) + #gam
  xlab("PseudoTime") +
  ylab("G2M.Score") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_identity(
    name = "Lineages",
    labels = c("Lineage1: Mono_CD16","Lineage2: Mono_CD16_ATG7+","Lineage3: Mono_CD14_ATG7+"),
    guide = "legend"
  ) +
  facet_wrap(~Condition, nrow = 3, scales = "free")

dev.off()

# ---- plot
pdf(file.path(results_path, "2_Mono_exp_trend_loess.pdf"), width = 16, height = 8)
ggplot(df) +
  geom_smooth(aes(x = Pseudotime1, y = expr, color = "#e20612"), method = "loess", se = F) + #gam
  geom_smooth(aes(x = Pseudotime2, y = expr, color = "#ffd401"), method = "loess", se = F) + #gam
  geom_smooth(aes(x = Pseudotime3, y = expr, color = "#00b0eb"), method = "loess", se = F) + #gam
  xlab("PseudoTime") +
  ylab("Expression") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_identity(
    name = "Lineages",
    labels = c("Lineage1: Mono_CD16","Lineage2: Mono_CD16_ATG7+","Lineage3: Mono_CD14_ATG7+"),
    guide = "legend"
  ) +
  facet_wrap(~feature, nrow = 4, scales = "free")
dev.off()

pdf(file.path(results_path, "2_Mono_exp_trend_condition_loess_3.pdf"), width = 8, height = 14)
# lineage 1
ggplot(df[!is.na(df$Pseudotime1), ]) +
  #geom_smooth(aes(x = Pseudotime1, y = expr, color = Condition), method = "gam", se = F) +
  geom_smooth(aes(x = Pseudotime1, y = expr, color = Condition), method = "loess", se = F) +
  xlab("Lineage1:  Mono_CD16") +
  ylab("Expression") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = n_condition) +
  facet_wrap(vars(feature), nrow = 6, scales = "free")
# lineage 2
ggplot(df[!is.na(df$Pseudotime2), ]) +
  geom_smooth(aes(x = Pseudotime2, y = expr, color = Condition), method = "loess", se = F) + #"gam"
  xlab("Lineage1:  Mono_CD16_ATG7+") +
  ylab("Expression") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = n_condition) +
  facet_wrap(vars(feature), nrow = 6, scales = "free")
# lineage 3
ggplot(df[!is.na(df$Pseudotime3), ]) +
  geom_smooth(aes(x = Pseudotime3, y = expr, color = Condition), method = "loess", se = F) + #"gam"
  xlab("Lineage3: Mono_CD14_ATG7+") +
  ylab("Expression") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = n_condition) +
  facet_wrap(vars(feature), nrow = 6, scales = "free")
dev.off()








