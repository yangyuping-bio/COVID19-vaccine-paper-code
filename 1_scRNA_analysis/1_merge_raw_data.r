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
results_path <- "/share/home/qlab/projects/project_cvv/yyp_results/4_Cell_Group/Tcell/CCA"
data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results/data"
# setwd()指定路径，setwd("~/project/")

#运行环境设置
source("/share/home/qlab/projects/project_cvv/scTools.R") #加载到当前的R环境中
set.seed(2)
system("hostname") #计算机主机名
use_python("~/miniconda3/envs/SC/bin/python", required = FALSE)
use_condaenv("SC",conda="~/miniconda3/bin/conda")
scr = import("scrublet")

#利用Scrublet对给定的数据进行双峰细胞检测，并将检测结果添加到数据对象的属性中
run_scrublet <- function(seu, npcs = 20, doublet_rate = 0.05){
  cells      = rownames(seu@meta.data) 
  mat        = t(as.matrix(seu[["RNA"]]@counts[,cells]))
  scrub      = scr$scrublet$Scrublet(mat, expected_doublet_rate= doublet_rate, sim_doublet_ratio = 2)
  res        = scrub$scrub_doublets(min_counts=2,
                                    min_cells=3,
                                    min_gene_variability_pctl=85,
                                    n_prin_comps=as.integer(npcs),
                                    log_transform = TRUE)
  names(res) = c("scrublet_score","scrublet_call ")
  print(paste0("Doublet number:",length(which(res$scrublet_call))))
  seu@meta.data$scrublet_score = res$scrublet_score
  seu@meta.data$scrublet_call  = res$scrublet_call
  return(seu)
}

#样本数据的循环处理
sample.combined  <- NULL
for(batch_id in unique(sample.df$Batch)){

  print(paste0("loading data ",batch_id))
  sample.raw <- Read10X_h5(paste0("./matrix/GEX/sc5r",batch_id,"gexhash/raw_feature_bc_matrix.h5"))
  sample.gex  <- sample.raw$`Gene Expression`
  sample.hash  <- sample.raw$`Antibody Capture`
  
  newnames  <-  paste0(batch_id,"_",colnames(sample.gex))
  colnames(sample.gex)  <- newnames
  colnames(sample.hash)  <- newnames
  
  print("calculate empty drop and filter..")
  
  empty.out <- emptyDrops(sample.gex, lower = 100, ignore = 200,BPPARAM  = DoparParam())
  is.cell <- empty.out$FDR <= 0.01
  sum(is.cell, na.rm=TRUE)
  table(Limited=empty.out$Limited, Significant=is.cell)
  cells.use  <- empty.out %>% as.data.frame() %>%
    rownames_to_column(var = "CellName") %>%
    dplyr::filter(FDR < 0.01, !is.na(FDR)) %>%
    pull(CellName)
  
  print("create..")
  sample <- CreateSeuratObject(counts = sample.gex[,cells.use], project = batch_id)
  sample  <- NormalizeData(object = sample, normalization.method = "LogNormalize",
                           scale.factor = 10000)
  hashtags  <- sample.df %>%
    dplyr::filter(Batch == batch_id) %>%
    dplyr::pull(Hashtag) %>%
    unique()
  
  sample[["HTO"]] <- CreateAssayObject(counts = sample.hash[hashtags,cells.use])
  sample <- NormalizeData(sample, assay = "HTO", normalization.method = "CLR")
  sample <- HTODemux(sample, assay = "HTO", positive.quantile = 0.95, nsamples = 1000)
  
  #load souporcell results
  soup.df  <- read_tsv(paste0("./matrix/GEX/sc5r",batch_id,"gexhash/clusters.tsv")) %>%
    dplyr::mutate(barcode = paste0(batch_id,"_",barcode)) %>%
    dplyr::select(barcode, status,assignment) %>%
    dplyr::filter(barcode %in% cells.use) %>%
    column_to_rownames(var ="barcode")
  
  sample  <- AddMetaData(object = sample, metadata = soup.df)
  
  metadata.df  <- sample@meta.data %>%
    dplyr::mutate(hash.ID = as.character(hash.ID)) %>%
    dplyr::filter(status == "singlet", HTO_classification.global == "Singlet")
  
  #match hashtags with soup results
  overlap.df  <- table(metadata.df$assignment, metadata.df$hash.ID) %>%
    as.data.frame.matrix() %>%
    tibble::rownames_to_column(var = "soup_class") %>%
    tidyr::gather(key = "hashtag", value = "num", -soup_class) %>%
    dplyr::group_by(hashtag) %>%
    dplyr::filter(num == max(num)) %>%
    dplyr::ungroup()
  
  match.df  <- sample.df %>%
    dplyr::filter(Batch == batch_id) %>%
    dplyr::left_join(overlap.df, by = c("Hashtag"="hashtag")) %>%
    dplyr::select(-num)
  
  dbl_rate  <- sample.df %>%
    dplyr::filter(Batch == batch_id) %>%
    dplyr::slice(1:1) %>%
    dplyr::pull(DblRate)
  
  sample  <- run_scrublet(sample,npcs = 20, doublet_rate = dbl_rate)
  sample@meta.data
  cells.doublet  <- sample@meta.data %>%
    tibble::rownames_to_column(var = "Cell") %>%
    top_frac(dbl_rate, scrublet_score) %>%
    pull(Cell)
  
  sample@meta.data$scrublet_callrate  <- "Singlet"
  sample@meta.data[cells.doublet, "scrublet_callrate"]  <- "Doublet"
  
  pk.all  <- 0.16
  expected_doublet_rate  <- dbl_rate
  sample <- NormalizeData(sample)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  sample <- ScaleData(sample, vars.to.regress = c("nCount_RNA"))
  sample  <- RunPCA(sample)
  nExp_poi <- round(expected_doublet_rate * length(rownames(sample@meta.data)))
  sample <- doubletFinder_v3(sample, PCs = 1:20, nExp = nExp_poi, pN = 0.25, pK = pk.all,
                             reuse.pANN = FALSE, sct = FALSE)
  num  <- ncol(sample@meta.data)
  sample@meta.data  <- sample@meta.data %>%
    dplyr::rename("DF.classifications" = num) %>%
    dplyr::rename("DF.pANN" = num -1)
  
  cells.new  <- sample@meta.data %>%
    tibble::rownames_to_column(var = "Cell") %>%
    dplyr::filter(assignment %in% c(0,1,2,3,4)) %>%
    dplyr::pull(Cell)
  
  sample.new  <- subset(sample, cells = cells.new)
  metadata.new  <- sample.new@meta.data %>%
    tibble::rownames_to_column(var = "Cell") %>%
    dplyr::left_join(match.df, by = c("assignment"="soup_class")) %>%
    tibble::column_to_rownames(var = "Cell")
  
  sample.new@meta.data  <- metadata.new
  
  if(is.null(sample.combined)) {
    sample.combined <- sample.new
  } else {
    sample.combined <- merge(sample.combined, sample.new)
    print("merge..")
  }
  print(paste0("done: batch ",batch_id))
  
}

sample.combined@meta.data$orig.ident%>% unique()
class(sample.combined)
saveRDS(sample.combined, file.path(data_save_path, "sample_raw_soup.rds"))
















