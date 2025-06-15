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
library(tradeSeq)
library(magrittr)
library(viridis)
library(Startrac)
library(UpSetR)
library(pheatmap)
library(Hmisc)
library(corrplot)


# 设置路径
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/3_Intergrated_CCA"
data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/2_QC_raw_data"
# setwd()指定路径，setwd("~/project/")

#运行环境设置
source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
set.seed(2)
system("hostname") #计算机主机名
use_python("~/miniconda3/envs/SC/bin/python", required = FALSE)
use_condaenv("SC",conda="~/miniconda3/bin/conda")
scr = import("scrublet")

sample.integrated <- readRDS(file.path(data_save_path, "2_QC_raw_umap.rds"))

#筛掉RNA_snn_res.2中51的细胞（HBB）
sample.integrated  <- subset(sample.integrated, RNA_snn_res.2 == 51, invert = TRUE)

sample.integrated@meta.data %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), umi = mean(nCount_RNA), mt = mean(percent.mt))

#干扰基因
hgGenes = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE","JCHAIN",
            "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
IG.gene1 <- grep(pattern = "^IGLV", x = rownames(x = sample.integrated), value = TRUE)
IG.gene2 <- grep(pattern = "^IGKV", x = rownames(x = sample.integrated), value = TRUE)
IG.gene3 <- grep(pattern = "^IGHV", x = rownames(x = sample.integrated), value = TRUE)
IG.genes <- union(IG.gene1,c(IG.gene2, IG.gene3))
TRDV.genes <- grep(pattern = "^TRDV", x = rownames(x = sample.integrated), value = TRUE)
TRAV.genes <- grep(pattern = "^TRAV", x = rownames(x = sample.integrated), value = TRUE)
TRBV.genes <- grep(pattern = "^TRBV", x = rownames(x = sample.integrated), value = TRUE)
TRGV.genes <- grep(pattern = "^TRGV", x = rownames(x = sample.integrated), value = TRUE)
TRV.genes <- union(TRDV.genes,c(TRAV.genes, TRBV.genes,TRGV.genes))
RPS.genes  <- grep(pattern = "^RPS", x = rownames(x = sample.integrated), value = TRUE)
RPL.genes  <- grep(pattern = "^RPL", x = rownames(x = sample.integrated), value = TRUE)
MT.genes  <- grep(pattern = "^MT-", x = rownames(x = sample.integrated), value = TRUE)
HIST.genes  <- grep(pattern = "^HIST", x = rownames(x = sample.integrated), value = TRUE)

DefaultAssay(sample.integrated)  <- "RNA"
sample.list  <- SplitObject(sample.integrated, split.by = "Sample")
var.genes  <- purrr::map(sample.list, function(sample) {
                sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2500)
                VariableFeatures(sample) }) %>%
                 purrr::reduce(union)

var.genes.filtered  <- setdiff(var.genes, c(IG.genes,TRV.genes, MT.genes, RPL.genes, RPS.genes,HIST.genes))

sample.list[c(8,14,34)]

# plan("multicore", workers = 4)
# plan()
# options(future.globals.maxSize = 80 * 1000 * 1024^2)

for (i in 1:length(sample.list)) {
    sample.list[[i]] <- NormalizeData(sample.list[[i]], verbose = FALSE)
    sample.list[[i]] <- FindVariableFeatures(sample.list[[i]], selection.method = "vst",
        nfeatures = 2500, verbose = FALSE)
    sample.list[[i]] <- ScaleData(sample.list[[i]],
                                  vars.to.regress = c("nFeature_RNA","percent.mt",
                                                      "G2M.Score","S.Score"),
                                  features = var.genes )
}

# plan("multicore", workers = 4)
# plan()
# options(future.globals.maxSize = 80 * 1000 * 1024^2)

sample.anchors <- FindIntegrationAnchors(sample.list,
                                         dims = 1:30,
                                         anchor.features = var.genes.filtered,
                                         # reference = c(8,11,13,14,34,35),#D1&D3-A/B/C
                                         reference = c(8,14,34), #D1-A/B/C
                                         reduction ="rpca",
                                         scale = FALSE)


# plan("multicore", workers = 1)
# plan()
# options(future.globals.maxSize = 500 * 1000 * 1024^2,future.seed=TRUE)

sample.integrated <- IntegrateData(anchorset = sample.anchors, dims = 1:30)

# plan("multicore", workers = 4)
# plan()
# options(future.globals.maxSize = 80 * 1000 * 1024^2)

sample.integrated <- FindVariableFeatures(sample.integrated, selection.method = "vst",
        nfeatures = 4000, verbose = FALSE)

var.features  <- VariableFeatures(sample.integrated)
# var.features  <- rownames(sample.integrated)
# plan("multicore", workers = 4)
# plan()
# options(future.globals.maxSize = 80 * 1000 * 1024^2)

sample.integrated <- ScaleData(sample.integrated,
                               # vars.to.regress = c("percent.mt","nFeature_RNA"),
                               features = rownames(sample.integrated))

genes.cc  <- extract_cellcycle(sample.integrated, var.features, cores =10, assay = "integrated", cutoff = 0.05)
genes.stress  <- extract_stress(sample.integrated, var.features, cores = 10, assay = "integrated", cutoff = 0.1)
genes.ifn  <- extract_ifn(sample.integrated, var.features, cores = 10, assay = "integrated", cutoff = 0.1)
genes.hypoxia  <- extract_hypoxia(sample.integrated, var.features, cores = 10, assay = "integrated", cutoff = 0.1)

var.features.filtered  <- var.features[! var.features %in% c(stress.genes,MT.genes, RPS.genes, RPL.genes, HIST.genes,
                                                             cc.genes$g2m.genes, cc.genes$s.genes, hgGenes,
                                                             genes.hypoxia[1:20],
                                                    interferon.genes, genes.ifn[1:30],genes.cc[1:200],genes.newcc,
                                                    genes.stress[1:20])]

# plan("multicore", workers = 10)
# plan()
# options(future.globals.maxSize = 40 * 1000 * 1024^2)

sample.integrated <- RunPCA(sample.integrated, npcs = 30, features = var.features.filtered)
sample.integrated <- RunUMAP(sample.integrated,
                             umap.method = "umap-learn",
                             reduction = "pca",
                             dims = 1:25)
sample.integrated <- FindNeighbors(sample.integrated,reduction = "pca",
                                   dims = 1:25, force.recalc = T)
sample.integrated <- FindClusters(sample.integrated, algorithm = 1,
                                  resolution = c(0.1, 0.2, 0.4, 0.8,1, 1.2, 1.6,2,3,4,8,10,12,16) )

saveRDS(sample.integrated, file.path(results_path, "3_CCA_Integrated_umap.rds"))

sample.integrated@meta.data %>%
  dplyr::group_by(integrated_snn_res.2) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), umi = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(umi) %>%
  print(n = 50)

# sample.integrated <- readRDS(file.path(data_save_path, "3_CCA_Integrated_umap.rds"))
#results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/3_Intergrated_CCA"

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.0.1","integrated_snn_res.0.2",
                      "integrated_snn_res.0.4","integrated_snn_res.0.8",
                      "integrated_snn_res.1","integrated_snn_res.1.2",
                      "integrated_snn_res.1.6","integrated_snn_res.2",
                      "Donor","Condition","Sample","Batch","HTO_classification.global",
                      "scrublet_callrate","DF.classifications")
pdf(file.path(results_path, "CCA_DimPlot.pdf"))
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.3","integrated_snn_res.4",
                      "integrated_snn_res.8","integrated_snn_res.10",
                      "integrated_snn_res.12","integrated_snn_res.16")
pdf(file.path(results_path, "CCA_DimPlot_3_10.pdf"))
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

Idents(sample.integrated) <- "integrated_snn_res.2"
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score","Stress.Score1","IFN.Score1","scrublet_score","DF.pANN")
pdf(file.path(results_path, "FeaturePlot_QC.pdf"))
for (plot_name in plot_list) {
  p <- FeaturePlot(object = sample.integrated,features = plot_name,cols = c("grey", "blue"), pt.size = 0.1,
                   reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
  print(p)
}
dev.off()


idents_name_list <- c("integrated_snn_res.2","Donor","Condition","Sample","Batch")
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score",
               "Stress.Score1","IFN.Score1","scrublet_score","DF.pANN","Hypoxia.Score1")

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

Idents(sample.integrated) <- "integrated_snn_res.2"
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

Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "CCA_FeaturePlot_Genes.pdf"), width = 40, height =60)
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

# 高亮细胞
Idents(sample.integrated) <- "integrated_snn_res.2"
idents_name_list <- unique(sample.integrated@meta.data$integrated_snn_res.2)
pdf(file.path(results_path, "Highlight_integrated_res.2.pdf"), 20, 10)
for (idents_name in idents_name_list) {
    p <- DimPlot(sample.integrated,reduction = 'umap',
                 pt.size =0.3,
            cells.highlight= WhichCells(sample.integrated, idents=idents_name),
            cols.highlight="red",
            #col="grey",
            group.by = "integrated_snn_res.2",
            label = T,
            raster=FALSE)+ggtitle(idents_name)
    print(p)
}
dev.off()






