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
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/4_Cell_type_level_1_2"
data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/3_Intergrated_CCA_Reumap"
# setwd()指定路径，setwd("~/project/")

# #运行环境设置
# source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
# set.seed(2)
# system("hostname") #计算机主机名
# use_python("~/miniconda3/envs/SC/bin/python", required = FALSE)
# use_condaenv("SC",conda="~/miniconda3/bin/conda")
# scr = import("scrublet")

sample.integrated <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results/data/3_Tem_Integrated_celltype.rds")

unique(sample.integrated$CellType)

sample.integrated$CellType[sample.integrated$CellType == "CD8_Tscm"] <- "CD8_Tcm"
sample.integrated$CellType[sample.integrated$CellType == "CD8_Tna"] <- "CD8_Tn"

DefaultAssay(sample.integrated) <- "RNA"
pdf(file.path(results_path, "DimPlot_CellType.pdf"), 15, 10)
  Idents(sample.integrated) <- "CellType"
  p <- DimPlot(object = sample.integrated, reduction = 'umap',
               label = T,ncol=4,label.size = 2,raster = FALSE) 
    print(p)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "FeaturePlot_T_cell.pdf"), width = 40, height=60)
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
pdf(file.path(results_path, "FeaturePlot_B_cell.pdf"), width = 40, height=40)
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
pdf(file.path(results_path, "FeaturePlot_NK_cell.pdf"), width = 40, height=40)
FeaturePlot(object = sample.integrated ,
            features = c("NCAM1","FCGR3A","CD14", "XCL2",
                         "CD3D", "CD4", "CD8A","CCR7", #T_cell
                         "GZMM", "GZMA", "GZMK", "GZMB",#毒性
                         "GZMH", "IL7R",  "TRDC","NKG7",
                         "GNLY", "PRF1", "CTSW", "RGS1",
                          "CREM", "IL32", 
                         "CCL2", "CCL3", "CCL5", "CXCL10", #炎症
                         "CXCL9", "IL6", "IL7", "IL15",
                         "CX3CR1",  "MKI67",
                         "PF4","SLC4A10","JCHAIN", # Doublet_others
                        ),
            cols = c("grey", "blue"), pt.size = 0.1, 
            reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
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

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
pdf(file.path(results_path, "FeaturePlot_Others_cell.pdf"), width = 40, height=30)
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
idents_name_list <- c("integrated_snn_res.2",
                      "integrated_snn_res.0.5",
                      "integrated_snn_res.1.5",
                      "integrated_snn_res.3",
                      "integrated_snn_res.5",
                       "integrated_snn_res.10",
                      "Donor","Condition","Sample","Batch","HTO_classification.global",
                      "scrublet_callrate","DF.classifications")
pdf(file.path(results_path, "DimPlot_all.pdf"),15,15)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

# 高亮细胞
Idents(sample.integrated) <- "integrated_snn_res.2"
idents_name_list <- unique(sample.integrated@meta.data$integrated_snn_res.2)
pdf(file.path(results_path, "Highlight_integrated_res.2.pdf"), 20, 10)
for (idents_name in idents_name_list) {
    p <- DimPlot(sample.integrated,reduction = 'umap',
                 #pt.size =0.3,
            cells.highlight= WhichCells(sample.integrated, idents=idents_name),
            cols.highlight="red",
            #col="grey",
            group.by = "integrated_snn_res.2",
            label = T,
            raster=FALSE)+ggtitle(idents_name)
    print(p)
}
dev.off()

unique(sample.integrated$integrated_snn_res.2)

DefaultAssay(sample.integrated) <- "RNA"                         
idents_name_list <- c("integrated_snn_res.2")
pdf(file.path(results_path, "DimPlot_res.2.pdf"),15,8)
for (idents_name in idents_name_list) {
  Idents(sample.integrated) <- idents_name
  p <- DimPlot(object = sample.integrated, reduction = 'umap',label = T,raster = FALSE) 
  print(p)
}
dev.off()

sample.integrated@meta.data %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(num = n(), gene = mean(nFeature_RNA), count = mean(nCount_RNA), mt = mean(percent.mt)) %>%
  arrange(num)

DefaultAssay(sample.integrated) <- "RNA"
pdf(file.path(results_path, "DimPlot_CellType_color.pdf"), 15, 10)
  Idents(sample.integrated) <- "CellType"
  p <- DimPlot(object = sample.integrated, reduction = 'umap',cols=allcolour,
               label = T,ncol=4,label.size = 2,raster = FALSE) 
    print(p)
dev.off()

allcolour <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
  "#968175"
)

umap = sample.integrated@reductions$umap@cell.embeddings %>%  #坐标信息
      as.data.frame() %>% 
      cbind(cell_type = sample.integrated@meta.data$CellType) # 注释后的label信息 ，改为cell_type
    
head(umap)
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
library(ggrepel)

pdf(file.path(results_path, "CellType_color.pdf"),15,10)
    p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
                    geom_point(size = 0.1 , alpha =1 )  + 
                    scale_color_manual(values = allcolour)
    # 调整umap图 - theme
    p2 <- p  +
      theme(panel.grid.major = element_blank(), #主网格线
            panel.grid.minor = element_blank(), #次网格线
            panel.border = element_blank(), #边框
            axis.title = element_blank(),  #轴标题
            axis.text = element_blank(), # 文本
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = 'white'), #背景色
            plot.background=element_rect(fill="white"))
    #调整umap图 - legend
    p3 <- p2 +         
            theme(
              legend.title = element_blank(), #去掉legend.title 
              legend.key=element_rect(fill='white'), #
            legend.text = element_text(size=10), #设置legend标签的大小
            legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
            guides(color = guide_legend(override.aes = list(size=5))) #设置legend中点的大小 
      #调整umap图 - annotation
    p4 <- p3 + 
      geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                       xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
                   colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
      geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                       xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
                   colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
      annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
               color="black",size = 3, fontface="bold" ) + 
      annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
               color="black",size = 3, fontface="bold" ,angle=90) 
       #p4
    print(p4)
    p5<- p4 +
        geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                       point.padding=unit(0.5, "lines")) +
                       theme(legend.position = "none")  #去掉legend
     print(p5)
dev.off()

DefaultAssay(sample.integrated) <- "RNA"
Idents(sample.integrated) <- "integrated_snn_res.2"
plan("multicore", workers = 4)
sample.integrated.markers <- FindAllMarkers(sample.integrated, only.pos = TRUE , min.pct = 0.15, logfc.threshold = 0.25)  #,min.pct = 0.1, logfc.threshold = 0.25
write.csv(sample.integrated.markers, file.path(results_path, "Markers_res.2.csv"))

#----细胞注释---------
# cells.use <- WhichCells(sample.integrated, expression = integrated_snn_res.2 == 10)
# sample.integrated@meta.data[cells.use, "ct_level_2"] <- "NK1"
Idents(sample.integrated)  <- "CellType"
sample.integrated@meta.data$ct_level_2  <- Idents(sample.integrated)
sample.integrated@meta.data$ct_level_2  <- as.character(sample.integrated@meta.data$ct_level_2)
sample.integrated$ct_level_2[sample.integrated$ct_level_2 == "Bmemory_CD27+/IGHM-"] <- "Bmemory_CD27+_IGHM-"
sample.integrated$ct_level_2[sample.integrated$ct_level_2 == "Mono_CD14/16"] <- "Mono_CD14_CD16"
sample.integrated$ct_level_2[sample.integrated$ct_level_2 == "Bmemory_CD27+/IGHM+"] <- "Bmemory_CD27+_IGHM+"

pdf(file.path(results_path, "ct_level_2_DimPlot.pdf"),12,8)
Idents(sample.integrated) <- "ct_level_2"
p <- DimPlot(object = sample.integrated, reduction = 'umap',cols=allcolour,label.size = 2,pt.size=0.05,label.box = T,
             label = T,raster = FALSE) 
print(p)
dev.off()

unique(sample.integrated@meta.data$ct_level_2)

unique(sample.integrated@meta.data$CellType)

CD4T <- c("CD4_Tn","CD4_Tcm","CD4_Tem","CD4_Treg")
CD8T <- c("CD8_Tn","CD8_Tcm","CD8_Tem","CD8_Temra","CD8_NELL2")
NK <- c("NK1","NK2","NK3")
Monocyte <- c("Mono_CD14","Mono_CD16","Mono_CD14_CD16","Mono_ATG7+")
Bcell <- c("Bmn_IGHD","Bmemory_CD27+_IGHM+","Bmemory_CD27+_IGHM-","Plasma")


Idents(sample.integrated)  <- "ct_level_2"
sample.integrated@meta.data$ct_level_1  <- Idents(sample.integrated)
sample.integrated@meta.data$ct_level_1  <- as.character(sample.integrated@meta.data$ct_level_1)
sample.integrated$ct_level_1[sample.integrated$ct_level_2 %in% CD4T] <- "CD4T"
sample.integrated$ct_level_1[sample.integrated$ct_level_2 %in% CD8T] <- "CD8T"
sample.integrated$ct_level_1[sample.integrated$ct_level_2 %in% NK] <- "NK"
sample.integrated$ct_level_1[sample.integrated$ct_level_2 %in% Bcell] <- "Bcell"
sample.integrated$ct_level_1[sample.integrated$ct_level_2 %in% Monocyte] <- "Monocyte"

# sample.integrated$ct_level_1 <- factor(
#   sample.integrated$ct_level_1,
#   levels = c("CD4T", "CD8T", "NK","Bcell","Monocyte"))
# table(sample.integrated$ct_level_1)

allcolour <- c(
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
  "#968175"
)

pdf(file.path(results_path, "ct_level_1_DimPlot.pdf"),12,8)
Idents(sample.integrated) <- "ct_level_1"
p <- DimPlot(object = sample.integrated, reduction = 'umap',cols=allcolour,label.size = 2,pt.size=0.05,label.box = T,
             label = T,raster = FALSE) 
print(p)
dev.off()

# ---- interferon genes
msigdb <- getBroadSets("/share/reference/sigdb/msigdb_v7.4.xml")
index <- sapply(msigdb, function(gs) {
  bcCategory(collectionType(gs)) == "c2"
})
geneset.c2 <- msigdb[index]
geneset.interferon <- geneset.c2[["BROWNE_INTERFERON_RESPONSIVE_GENES"]]
interferon.genes <- geneset.interferon@geneIds
interferon.genes <- interferon.genes[interferon.genes %in% rownames(sample.integrated)]
cellcycle.genes <- geneset.c2[["KEGG_CELL_CYCLE"]]@geneIds
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(interferon.genes),
  ctrl = length(interferon.genes),
  name = "IFN.Score.sample.integrated",
  assay = "RNA"
)

# ---- 32 stress genes
stress.genes <- read_csv("/share/home/qlab/projects/project_cvv/data/stress.genes.csv") %>%
  dplyr::filter(gene %in% rownames(sample.integrated)) %>%
  dplyr::pull(gene)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(stress.genes),
  ctrl = length(stress.genes),
  name = "Stress.Score.sample.integrated",
  assay = "RNA"
)

# ---- hypoxia genes
genes.hypoxia <- c(
  "VEGFA", "SLC2A1", "PGAM1", "ENO1",
  "LDHA", "TPI1", "P4HA1", "MRPS17",
  "CDKN3", "ADM", "NDRG1", "TUBB6",
  "ALDOA", "MIF", "ACOT7"
)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.hypoxia),
  ctrl = length(genes.hypoxia),
  name = "Hypoxia.Score.sample.integrated",
  assay = "RNA"
)

# ---- variable genes
sample.integrated <- FindVariableFeatures(sample.integrated, selection.method = "vst", nfeatures = 3000)
var.genes <- VariableFeatures(sample.integrated)
# HG genes
hgGenes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
# IG genes
igGenes <- c(
  "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "JCHAIN",
  "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC"
)
IG.gene1 <- grep(pattern = "^IGLV", x = rownames(x = sample.integrated), value = TRUE)
IG.gene2 <- grep(pattern = "^IGKV", x = rownames(x = sample.integrated), value = TRUE)
IG.gene3 <- grep(pattern = "^IGHV", x = rownames(x = sample.integrated), value = TRUE)
IG.genes <- union(IG.gene1, c(IG.gene2, IG.gene3))
# TRV genes
TRDV.genes <- grep(pattern = "^TRDV", x = rownames(x = sample.integrated), value = TRUE)
TRAV.genes <- grep(pattern = "^TRAV", x = rownames(x = sample.integrated), value = TRUE)
TRBV.genes <- grep(pattern = "^TRBV", x = rownames(x = sample.integrated), value = TRUE)
TRGV.genes <- grep(pattern = "^TRGV", x = rownames(x = sample.integrated), value = TRUE)
TRV.genes <- union(TRDV.genes, c(TRAV.genes, TRBV.genes, TRGV.genes))
# RPS & RPL
RPS.genes <- grep(pattern = "^RPS", x = rownames(x = sample.integrated), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(x = sample.integrated), value = TRUE)
# MT genes
MT.genes <- grep(pattern = "^MT-", x = rownames(x = sample.integrated), value = TRUE)
# HIST genes
HIST.genes <- grep(pattern = "^HIST", x = rownames(x = sample.integrated), value = TRUE)
# select variable genes
var.genes <- var.genes[!var.genes %in% c(IG.genes, TRV.genes,MT.genes, RPL.genes, RPS.genes,HIST.genes)]

sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(TRV.genes),
  ctrl = length(TRV.genes),
  name = "TRV.Score.sample.integrated",
  assay = "RNA"
)
genes.via <- c("GMFG","CD27","ACTG1", "LIMD2", "CD74","ARHGDIB",
       "GZMA", "COTL1", "RTRAF5", "WDR1",  "ATP5MC2",
       "ACTB", "PSMB8", "CNN2", "ATP5MG",
       "HIGD2A", "PSMA5", "TMSB4X", "PYCARD",
        "MYL6", "ARPC1B", "HINT", "CCL5","ACTR3", "PSMA4")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.via),
  ctrl = length(genes.via),
  name = "via.Score.sample.integrated",
  assay = "RNA"
)
genes.vib <- c("PKMYT1","E2F7","SPC25", "CCNB2", "CDC45","FAM111B",
       "MYBL2", "DTL", "CLSPN", "MND1",  "TYMS",
       "CDK1", "DLGAP5", "CDCA5", "CDC6",
       "HIST1H3G", "TK1", "CKAP2L", "CDT1",
        "MKI67", "CENPF", "ZWINT", "PCLAF","KNL1", "RRM2",
        "CENPU", "PBK", "MCM10", "CEP55","KLRG1")
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.vib),
  ctrl = length(genes.vib),
  name = "vib.Score.sample.integrated",
  assay = "RNA"
)
# ---- add module score
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(cc.genes$s.genes),
  ctrl = length(cc.genes$s.genes),
  name = "S.Score.sample.integrated",
  search = FALSE,
  assay = "integrated"
)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(cc.genes$g2m.genes),
  ctrl = length(cc.genes$g2m.genes),
  name = "G2M.Score.sample.integrated",
  search = FALSE,
  assay = "integrated"
)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(interferon.genes),
  ctrl = length(interferon.genes),
  name = "IFN.Score.sample.integrated",
  assay = "integrated"
)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(stress.genes),
  ctrl = length(stress.genes),
  name = "Stress.Score.sample.integrated",
  assay = "integrated"
)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = list(genes.hypoxia),
  ctrl = length(genes.hypoxia),
  name = "Hypoxia.Score.sample.integrated",
  assay = "integrated"
)


Idents(sample.integrated) <- "integrated_snn_res.0.5"
plot_list <- c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score",
               "scrublet_score",
               "via.Score.sample.integrated1","vib.Score.sample.integrated1","TRV.Score.sample.integrated1",
               "Stress.Score.sample.integrated1","IFN.Score.sample.integrated1","DF.pANN",
               "Hypoxia.Score.sample.integrated1",
               "S.Score.sample.integrated1","G2M.Score.sample.integrated1")
pdf(file.path(results_path, "CCA_FeaturePlot_QC.pdf"))
for (plot_name in plot_list) {
  p <- FeaturePlot(object = sample.integrated,features = plot_name,cols = c("grey", "blue"), pt.size = 0.1,
                   reduction = "umap",raster=FALSE,label= TRUE,label.size = 3)
  print(p)
}
dev.off()

str(sample.integrated@meta.data)

saveRDS(sample.integrated, file.path(results_path, "4_integrated_level_1_2.rds"))
