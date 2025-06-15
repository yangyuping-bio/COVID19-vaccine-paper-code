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
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/1_Tcell_Left"
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
targetcell <- subset(targetcell, idents = c("CD4T","CD8T","ILC",
                                           "ILC2","gdT_V9","gdT","Tdn",
                                           "MAIT","NK")) 

targetcell@meta.data %>%
  dplyr::group_by(ct_level_1) %>%
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
plan("multicore", workers = 3)
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
options(future.globals.maxSize = 200 * 1024^3)
gc()
plan()

sample.integrated <- IntegrateData(
  anchorset = sample.anchors,
  dims = 1:30,
  k.weight = 100, #default:k.weight = 100
  verbose = TRUE
)

# ---- scale expression
sample.integrated <- FindVariableFeatures(sample.integrated, selection.method = "vst", nfeatures = 4000, verbose = F)
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

# add moudule score
genes.via <- c("ACTG1","PCLAF","ACTB","GAPDH","CD38","GZMH","PFN1","CORO1A","PYCARD","PSMB8","TMSB10","ARPC3",
          "CLIC1","CD52","TMSB4X","APOBEC3H","ADA","SEM1","MYL6","GZMA","DECR1","ARPC2","NOP10","APOBEC3C","LGALS1","GZMB",
          "CFL1","CLTA","TWF2","ISG15","MXD4","SUB1","PSMA4","HLA-DRA","PPP1R18","GSTO1","APOBEC3G","ARPC5","UCP2",
          "PSMA5","WDR1","ARHGDIB","YWHAE","CCL5","SLC25A5","PSMB4","ITGB1","DYNLL1","ARPC4","ATP5MC2")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.via),
  ctrl = length(genes.via),
  name = "VIA.Score.Tcell",
  assay = "RNA"
)
genes.vib <- c("TYMS","PCNA","CLSPN","MCM4","NASP","RRM2","GINS2","FEN1","MCM2","SLBP",
         "HELLS","USP1","DTL",
         "RRM1","CDC6","MCM6","GMNN","HMGB2","TOP2A","CKS1B","MKI67","TMPO","CDK1","CCNB2")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.vib),
  ctrl = length(genes.vib),
  name = "VIB.Score.Tcell",
  assay = "RNA"
)

# ---- PCA and clustering
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- RunPCA(sample.integrated, npcs = 30, features = var.features.filtered)
sample.integrated <- RunUMAP(sample.integrated, umap.method = "umap-learn", reduction = "pca", dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1, resolution = seq(0.1, 3, 0.1))

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

Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "CCA_DimPlot_Donor.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Donor',cols=allcolors,
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "CCA_DimPlot_Batch.pdf"),15,15)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Batch',cols=allcolors,
             ncol = 3,label = T,raster = FALSE)
print(p)
dev.off()
pdf(file.path(results_path, "CCA_DimPlot_Sample.pdf"),25,25)
p <- DimPlot(object = sample.integrated, reduction = 'umap',split.by = 'Sample',cols=allcolors,
             ncol = 6,label = T,raster = FALSE)
print(p)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.0.7","ct_level_2","ct_level_1","integrated_snn_res.1")
pdf(file.path(results_path, "DimPlot_res.0.7.pdf"),7.5,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 4,raster = FALSE) 
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
pdf(file.path(results_path, "DimPlot_all.pdf"),7,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               label = T,label.size = 2,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_T_cell.pdf"), width = 30, height=50)
FeaturePlot(object = sample.integrated ,
            features = c("CD3D","CD3E","CD4","CD8A", # T_cell
                         "CCR7","CD27", # naive_T
                          "TRDV2","TRGV9", "TRGV10","TRDV1",# gdT
                          "CD25","FOXP3","CTLA4", # Treg
                         "ICOS", # CD4_Tfh
                         "IL17A",# Th17
                         "IFNG","CXCR3","CCR5", #Th1
                         "CCR6", #检查CD4_CCR6
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

Idents(sample.integrated) <- "integrated_snn_res.0.5"
pdf(file.path(results_path, "FeaturePlot_others.pdf"), width = 30, height =40)
FeaturePlot(object = sample.integrated,
            features = c("CD3D","CD4","CD8A","SLC4A10",
                        "CCR7","GZMA","GZMK","GZMB",
                         "JCHAIN","CD1C", "PF4","PPBP","FOXP3",
                         "TRDC","TRGV9", "TRGV10","TRDV1",
                         "CD79A","TCL1A","IGHA1","MZB1",
                         "CD27","IGHG1","IGHM","IGHD","SOX5",
                         "CD14","FCGR3A","FCGR3B","LILRA4",
                         "CD68","CD34","CFP","STAB1","MKI67"),
            cols = c("grey", "blue"), pt.size = 0.1, reduction = "umap",
            raster=FALSE,label= TRUE,label.size = 3)
dev.off()

str(sample.integrated@meta.data)

idents_name_list <- c("integrated_snn_res.0.5","ct_level_2","Donor","Condition")
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","VIA.Score.Tcell1","VIB.Score.Tcell1","S.Score","G2M.Score")
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
plan("multicore", workers = 12)
sample.integrated.markers <- FindAllMarkers(sample.integrated, only.pos = TRUE , min.pct = 0.15, logfc.threshold = 0.25)  #,min.pct = 0.1, logfc.threshold = 0.25
write.csv(sample.integrated.markers, file.path(results_path, "Markers_res.0.5.csv"))

#----细胞注释---------
Idents(sample.integrated)  <- "integrated_snn_res.0.5"
sample.integrated  <- RenameIdents(sample.integrated,
                            '13' =  "D_Mono_T",
                           '3' =  "NK1",
                            '2' =  "NK2",
                           '16' =  "NK3",
                           '9' =  "MAIT",
                           '12' =  "ILC2",
                           '7' =  "ILC",
                           '8' =  "gdT_V9",
                          '14' =  "gdT",
                            '18' =  "gdT",
                           '11' =  "CD4_Tem",
                           '0' =  "CD4_Tn",
                           '19' =  "Tdn",
                           '6' =  "CD4_Tcm",
                                   '5' =  "CD4_Tcm",
                            '10' =  "CD4_Treg",
                           '15' =  "CD8_NELL2",
                           '4' =  "CD8_Tem",
                           '1' =  "CD8_Tn",
                                   '20' =  "D_T_B",
                                   '17' =  "D_T_B"
                           )
sample.integrated@meta.data$ct_Tcell  <- Idents(sample.integrated)
sample.integrated@meta.data$ct_Tcell <- as.character(sample.integrated@meta.data$ct_Tcell)

sample.integrated$ct_Tcell[sample.integrated$integrated_snn_res.1 == 10] <- "CD8_Tcm"
sample.integrated$ct_Tcell[sample.integrated$integrated_snn_res.1 == 29] <- "MAIT"
sample.integrated$ct_Tcell[sample.integrated$integrated_snn_res.0.7 == 24] <- "D_T_Platelet"

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "ct_Tcell"
plan("multicore", workers = 12)
sample.integrated.markers <- FindAllMarkers(sample.integrated, only.pos = TRUE , min.pct = 0.15, logfc.threshold = 0.25)  #,min.pct = 0.1, logfc.threshold = 0.25
write.csv(sample.integrated.markers, file.path(results_path, "Markers_ct_Tcell.csv"))

# 高亮细胞
Idents(sample.integrated) <- "ct_Tcell"
idents_name_list <- unique(sample.integrated@meta.data$ct_Tcell)
pdf(file.path(results_path, "Highlight_ct_Tcell.pdf"), 7.5, 5)
for (idents_name in idents_name_list) {
    p <- DimPlot(sample.integrated,reduction = 'umap',
                 #pt.size =0.3,
            cells.highlight= WhichCells(sample.integrated, idents=idents_name),
            cols.highlight="red",
            #col="grey",
            group.by = "ct_Tcell",
            label = T,
                 pt.size = 0.01,
            raster=FALSE)+ggtitle(idents_name)
    print(p)
}
dev.off()

unique(sample.integrated$ct_Tcell)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("ct_Tcell","integrated_snn_res.0.5","ct_level_2","ct_level_1")
pdf(file.path(results_path, "DimPlot_ct_Tcell.pdf"),7.5,5)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap', cols=allcolors,
               pt.size = 0.01,label = T,label.size = 4,raster = FALSE) 
  print(p)
}
dev.off()

sample.integrated@meta.data %>%
  dplyr::group_by(ct_Tcell) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)

saveRDS(sample.integrated, file = "CCA_Tcell_Batch.rds")

DefaultAssay(sample.integrated) <- "RNA"   
Idents(sample.integrated) <- "ct_Tcell"
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score",
               "VIA.Score.Tcell1","VIB.Score.Tcell1",
               "Stress.Score1","IFN.Score1","scrublet_score","DF.pANN")
pdf(file.path(results_path, "FeaturePlot_QC_VIA.pdf"),10,6)
for (plot_name in plot_list) {
  p <- FeaturePlot(object = sample.integrated,features = plot_name,cols = c("grey", "blue"), pt.size = 0.1,
                   reduction = "umap",raster=FALSE,label= TRUE,label.size = 2)
  print(p)
}
dev.off()


