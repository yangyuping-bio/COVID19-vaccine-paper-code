library("tidyverse")
library("data.table")
library(ggplot2)
library("Matrix")
library(Seurat)
library(viridis)
library(RColorBrewer)

# 文件路径列表
bcr_files <- c(
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX10bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX11bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX12bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX50bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX51bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX52bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX53bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX5bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX6bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX7bcr/outs/filtered_contig_annotations.csv',
    '/share/home/qlab/project_10x/3_data/project_cvv//sc5rMIX8bcr/outs/filtered_contig_annotations.csv'
)

# add function
add_clonotype <- function(file_dir, type) {
  file_dir <- gsub("/outs/filtered_contig_annotations.csv", "", file_dir) #字符串替换为空字符串
  mix_name <- gsub("/share/home/qlab/project_10x/3_data/project_cvv/", "", file_dir)
  mix_name <- gsub("sc5r", "", mix_name)
  mix_name <- gsub(type, "", mix_name)

  tcr <- read.csv(paste0(file_dir, "/outs/filtered_contig_annotations.csv"), sep = ",")
  tcr$barcode <- paste(mix_name, tcr$barcode, sep = "_") #合并为新列
  print(paste0(mix_name,
               "  unique:  ",length(unique(tcr$barcode)),
               "  all:  ",length(tcr$barcode),
               "  duplicated:  ",length(tcr$barcode)-length(unique(tcr$barcode)),
               "  pro:  ",length(unique(tcr$barcode))/length(tcr$barcode)
              ))
  tcr <- tcr[!duplicated(tcr$barcode), ] #筛选出不重复的barcode
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

  # Clonotype-centric info
  clono <- read.csv(paste0(file_dir, "/outs/clonotypes.csv"))

  # Slap the AA sequences onto our original table by clonotype_id.
  #tcr <- merge(tcr[, c("barcode", "clonotype_id")], clono)
  tcr <- merge(tcr, clono)

  # Reorder so barcodes are first column and set them as rownames.
  tcr$chain_cdr3 <- paste(tcr$chain, tcr$cdr3, sep = "_")
  tcr$clonotype_mix <- paste(mix_name, tcr$clonotype_id, sep = "_")
  colnames(tcr) <- paste(type, colnames(tcr), sep = "_")
  return(tcr)
}


# ---- add BCR metadata
bcr_metadata <- lapply(1:length(bcr_files), function(x) {
  add_clonotype(bcr_files[x], type = "bcr")
})
bcr_metadata <- rbindlist(bcr_metadata)
bcr_metadata$bcr_barcode <- gsub("/", "", bcr_metadata$bcr_barcode)
bcr_metadata$bcr_clonotype_mix <- gsub("/", "", bcr_metadata$bcr_clonotype_mix)
rownames(bcr_metadata) <- bcr_metadata$bcr_barcode
#bcr_metadata[, "bcr_barcode"] <- NULL

colnames(bcr_metadata)

data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3"
Bcell <- readRDS(file.path(data_save_path, "2_Bcell/2_Reumap_with_Plasma/Reumap_Bcell_with_Plasma.rds"))

# ---- add TCR clonotype, rename bcr_clonotype_id based on cdr3s_aa
Bcell <- AddMetaData(object = Bcell, metadata = bcr_metadata)
table(!is.na(Bcell$bcr_clonotype_id))
change_id <- data.frame(raw = na.omit(unique(Bcell$bcr_cdr3s_aa)))
change_id$new <- paste0("clonotype", 1:nrow(change_id))
Bcell$bcr_clonotype_id_new <- change_id$new[match(Bcell$bcr_cdr3s_aa, change_id$raw)]
table(Bcell$ct_level_4, Bcell$bcr_clonotype_id_new > 0)
head(change_id)

saveRDS(bcr_metadata, "data/1_bcr_metadata.rds")
write.csv(bcr_metadata, "data/1_bcr_metadata.csv")
saveRDS(Bcell, "data/1_Bcell_BCR.rds")

getwd()

Bcell <- readRDS("data/1_Bcell_BCR.rds")

unique(Bcell$ct_level_4)
unique(Bcell$Condition)

Bcell$Condition <- factor(Bcell$Condition,
                         levels=c('A0','A1','A2',
                                  'B0','B1','B2',
                                  'C0','C1','C2'))
unique(Bcell$Condition)

# 定义函数来处理细胞对象的克隆型数据
process_clonotype_data <- function(cell_type_obj,ct_seq) {
  # 设置 ct_level_4 因子顺序
  cell_type_obj$ct_level_4 <- factor(
    cell_type_obj$ct_level_4,
    levels = ct_seq
  )
  
  # 设置 Donor 因子顺序
  cell_type_obj$Donor <- factor(
    cell_type_obj$Donor,
    levels = unique(cell_type_obj$Donor)
  )

  # 初始化 Clonotype_num_ld 列
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

nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
ct_seq <- c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-','Plasma')
clonotype_num <- as.data.frame(table(Bcell$bcr_clonotype_id_new))
Bcell$Clonotype_num <- clonotype_num$Freq[match(Bcell$bcr_clonotype_id_new, clonotype_num$Var1)]
Bcell <- process_clonotype_data(Bcell,ct_seq)

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(Bcell$ct_level_4, Bcell$Clonotype_num_ld))
colnames(data) <- c("ct_level_4", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/4_Bcell_Bar_Clonotype_Levels_all.pdf",3.5,2.8)
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
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
      plot.margin = unit(c(2, 2, 2, 2), "char"),
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1,color='black')# 设置 X 轴文本
  ) 
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
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
      plot.margin = unit(c(2, 2, 2, 2), "char"),
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1,color='black')# 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)
Bcell_f <- Bcell@meta.data%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
Bcell_f$Clonotype_num_ld <- factor(Bcell_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(Bcell_f$ct_level_4, Bcell_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("ct_level_4", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/4_Bcell_Bar_Clonotype_Levels.pdf",3.5,2.8)
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
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
       plot.margin = unit(c(2, 2, 2, 2), "char"),
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1,color='black'), # 设置 X 轴文本
  ) 
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
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
      plot.margin = unit(c(2, 2, 2, 2), "char"),
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1,color='black')# 设置 X 轴文本
  ) 
dev.off()

color_mapping <-  c("Not detected" ="lightgrey", 
                      "n = 1" ="darkgrey", 
                      "1 < n <= 3" ="#F3B1A0", 
                      "3 < n <= 5" ="#B22C2CFF",
                      "n > 5" ="#573333FF")
data <- as.data.frame(table(Bcell$Condition, Bcell$Clonotype_num_ld))
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")
head(data)

pdf("figures/4_Bcell_Bar_Clonotype_Levels_Condition_all.pdf",3.6,2.5)
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
    axis.text.x = element_text(size = 8,color='black')# 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)
Bcell_f <- Bcell@meta.data%>% filter(Clonotype_num_ld %in% c("1 < n <= 3", "3 < n <= 5",  "n > 5"))
Bcell_f$Clonotype_num_ld <- factor(Bcell_f$Clonotype_num_ld,
                     levels = rev(c( "1 < n <= 3", "3 < n <= 5", "n > 5")))
data <- as.data.frame(table(Bcell_f$Condition, Bcell_f$Clonotype_num_ld))

# 将列重命名
colnames(data) <- c("Condition", "Clonotype_num_ld", "Count")

# 绘制堆叠条形图
pdf("figures/4_Bcell_Bar_Clonotype_Levels_Condition.pdf",3.6,2.5)
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
    axis.text.x = element_text(size = 8,color='black'), # 设置 X 轴文本
  ) 
dev.off()





# 定义函数来处理细胞对象的克隆型数据
process_clonotype_data <- function(cell_type_obj,ct_seq) {
  # 设置 ct_level_4 因子顺序
  cell_type_obj$ct_level_4 <- factor(
    cell_type_obj$ct_level_4,
    levels = ct_seq
  )
  
  # 设置 Donor 因子顺序
  cell_type_obj$Donor <- factor(
    cell_type_obj$Donor,
    levels = c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5", "Donor6", "Donor7", "Donor8",
               "Donor9", "Donor10", "Donor11", "Donor12")
  )
  
  # 设置 Batch 因子顺序
  cell_type_obj$Batch <- factor(
    cell_type_obj$Batch,
    levels = c("MIX5", "MIX6", "MIX7", "MIX8", "MIX10", "MIX11", "MIX12", "MIX50", "MIX51",
               "MIX52", "MIX53")
  )
  
  # 初始化 Clonotype_num_ld 列
  cell_type_obj$Clonotype_num_ld <- "Not detected"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num == 1] <- "n = 1"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num == 2] <- "n = 2"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num == 3] <- "n = 3"
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 3] <- "n > 3"
  
  # 将 Clonotype_num_ld 转换为因子，并设定顺序
  cell_type_obj$Clonotype_num_ld <- factor(
    cell_type_obj$Clonotype_num_ld,
    levels = c("Not detected", "n = 1", "n = 2", "n = 3",  "n > 3")
  )
  
  return(cell_type_obj)
}

# 定义通用绘图函数
plot_clonotype_umap <- function(cell_type_obj, file_name, nclos_clonotype,w,h,pt) {
  pdf(file_name, w,h)
  print(
    DimPlot(
      object = cell_type_obj, reduction = "umap", label = FALSE, pt.size = pt, order = TRUE,
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
      pt.size = 0.01, 
      order = TRUE,
      group.by = "Clonotype_num_ld", 
      split.by = split_by,
      cols = nclos_clonotype
    )
  )
}


# 函数 调用：
nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
ct_seq <-c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-','Plasma')
clonotype_num <- as.data.frame(table(Bcell$bcr_clonotype_id_new))
Bcell$Clonotype_num <- clonotype_num$Freq[match(Bcell$bcr_clonotype_id_new, clonotype_num$Var1)]
Bcell <- process_clonotype_data(Bcell,ct_seq)

plot_clonotype_umap(Bcell, "figures/2_Bcell_clonotype_umap.pdf", nclos_clonotype,5,3.5,0.0001)

pdf("figures/2_Bcell_clonotype_umap_split.pdf", 18, 4)
split_conditions <- c("Clonotype_num_ld", "Condition", "ct_level_4", "Donor", "Batch")
lapply(split_conditions, function(split_by) plot_dim(Bcell, split_by))
dev.off()

nclos_clonotype <-  c("lightgrey", "darkgrey", 
                      "#F3B1A0",  "#B22C2CFF","#573333FF")
ct_seq <-c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-','Plasma')
clonotype_num <- as.data.frame(table(Bcell$bcr_clonotype_id_new))
Bcell$Clonotype_num <- clonotype_num$Freq[match(Bcell$bcr_clonotype_id_new, clonotype_num$Var1)]
Bcell <- process_clonotype_data(Bcell,ct_seq)

table(Bcell$Clonotype_num)

table(Bcell$Clonotype_num_ld)

table(Bcell$ct_level_4,Bcell$Clonotype_num_ld)

table(Bcell$Condition,Bcell$Clonotype_num_ld)



# 载入 ggplot2 包
library(ggplot2)

# 创建一个数据框用于绘制图形
Bcell$ct_level_4 <- factor(Bcell$ct_level_4,
                     levels = c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-','Plasma'))
Bcell$Clonotype_num_ld <- factor(Bcell$Clonotype_num_ld,
                     levels = rev(c("Not detected", "n = 1", "n = 2", "n = 3",  "n > 3")))
data <- as.data.frame(table(Bcell$ct_level_4, Bcell$Clonotype_num_ld))

# 将列重命名为更友好的名称
colnames(data) <- c("ct_level_4", "Clonotype_num_ld", "Count")

# 设置颜色对应
# color_mapping <- c("Not detected"="#57C3F3",
#                    "n = 1"="#CCE0F5", 
#                    "n = 2" ="#DDC0DC", 
#                    "n = 3" ="#BE84BC", 
#                    "n > 3" ="#810F7C")
color_mapping <- c("Not detected"="lightgrey",
                   "n = 1"="darkgrey", 
                   "n = 2" ="#F3B1A0", 
                   "n = 3" ="#B22C2CFF", 
                   "n > 3" ="#573333FF")
# # 绘制堆叠条形图
# ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = color_mapping) +
#   labs(title = "Stacked Bar Plot of Clonotype Levels", 
#        x = "CT Level 4", 
#        y = "Count", 
#        fill = "Clonotype Number Level") +
#   theme_minimal()
# 绘制堆叠条形图
pdf("figures/2_Bar_Plot_Clonotype_Levels_all.pdf",3,3)
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
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
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1), # 设置 X 轴文本
  ) 
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
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
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1), # 设置 X 轴文本
  ) 
dev.off()

# 载入 ggplot2 包
library(ggplot2)

# 创建一个数据框用于绘制图形
Bcell$ct_level_4 <- factor(Bcell$ct_level_4,
                     levels = c('Bmn_IGHD','Bm_CD27+_IGHM+',
                                'Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-','Plasma'))
Bcell_f <- Bcell@meta.data%>% filter(Clonotype_num_ld %in% c("n = 2", "n = 3",  "n > 3"))
Bcell_f$Clonotype_num_ld <- factor(Bcell_f$Clonotype_num_ld,
                     levels = rev(c( "n = 2", "n = 3",  "n > 3")))
data <- as.data.frame(table(Bcell_f$ct_level_4, Bcell_f$Clonotype_num_ld))

# 将列重命名为更友好的名称
colnames(data) <- c("ct_level_4", "Clonotype_num_ld", "Count")

# 设置颜色对应
# color_mapping <- c(
#                    "n = 2" ="#DDC0DC", 
#                    "n = 3" ="#BE84BC", 
#                    "n > 3" ="#810F7C")
color_mapping <- c(
                   "n = 2" ="#F3B1A0", 
                   "n = 3" ="#B22C2CFF", 
                   "n > 3" ="#573333FF")
# 绘制堆叠条形图
pdf("figures/2_Bar_Plot_Clonotype_Levels.pdf",3,3)
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
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
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1), # 设置 X 轴文本
  ) 
dev.off()

# 绘制堆叠条形图
pdf("figures/2_Bar_Plot_Clonotype_Levels_count.pdf",3,3)
ggplot(data, aes(x = ct_level_4, y = Count, fill = Clonotype_num_ld)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Clonotype Number", 
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
    axis.text.x = element_text(size = 8, angle = 45,hjust = 1), # 设置 X 轴文本
  ) 
dev.off()



clonotypes <- data.frame(
  id = unique(Bcell$bcr_clonotype_id_new),
  Bmn_IGHD = 0,
  Bm_CD27_pos_IGHM_pos = 0,
  Bm_CD27_pos_IGHM_pos_SOX5 = 0,
  Bm_CD27_pos_IGHM_neg = 0,
  Plasma = 0
)
clonotypes <- clonotypes[!is.na(clonotypes$id), ]
for (i in 1:nrow(Bcell)) {
  id <- Bcell$bcr_clonotype_id_new[i]
  if (is.na(id)) {
    next
  }
  CellType <- as.character(Bcell$ct_level_4[i])
  # 修改 CellType 以匹配 clonotypes 中的列名
  CellType <- gsub("\\+", "_pos", CellType)
  CellType <- gsub("-", "_neg", CellType)
  clonotypes[clonotypes$id == id, CellType] <- clonotypes[clonotypes$id == id, CellType] + 1
}
clonotypes[clonotypes > 0] <- 1
head(clonotypes)

pdf("figures/3_Bcell_upset_ct_level_4.pdf", width = 8, height = 6)
upset(clonotypes,sets = c('Bmn_IGHD','Bm_CD27_pos_IGHM_pos',
                          'Bm_CD27_pos_IGHM_pos_SOX5','Bm_CD27_pos_IGHM_neg','Plasma'),
            keep.order = TRUE, text.scale = 1,
            mb.ratio = c(0.7, 0.3),
            matrix.color = "#B22C2CFF", main.bar.color = "#B22C2CFF",
           shade.color = "#F3B1A0",sets.bar.color = "#B22C2CFF", 
            point.size = 2, line.size = 0.8)
dev.off()



# ---- upset plot
library(UpSetR)
clonotypes <- data.frame(
  id = unique(Bcell$bcr_clonotype_id_new),
  Donor1 = 0,
  Donor2 = 0,
  Donor3 = 0,
  Donor4 = 0,
  Donor5 = 0,
  Donor6 = 0,
  Donor7 = 0,
  Donor8 = 0,
  Donor9 = 0,
  Donor10 = 0,
  Donor11 = 0,
  Donor12 = 0
)
clonotypes <- clonotypes[!is.na(clonotypes$id), ]

for (i in c(1:ncol(Bcell))) {
  #print(i)
  id <- Bcell$bcr_clonotype_id_new[i]
  if (is.na(id)) {
    next
  }
  donor_n <- as.character(Bcell$Donor[i])
  clonotypes[clonotypes$id == id, donor_n] <- clonotypes[clonotypes$id == id, donor_n] + 1
}
clonotypes[clonotypes > 0] <- 1

pdf("figures/3_Bcell_upset_donor.pdf", width = 10, height = 8)
upset(clonotypes,sets = c("Donor1" ,"Donor2","Donor3" ,"Donor4","Donor5" ,"Donor6","Donor7" ,"Donor8",
                          "Donor9" ,"Donor10","Donor11" ,"Donor12"),
            keep.order = TRUE, text.scale = 1,
            matrix.color = "#B22C2CFF", main.bar.color = "#B22C2CFF",
           shade.color = "#F3B1A0",sets.bar.color = "#B22C2CFF", 
            mb.ratio = c(0.6, 0.4),
            point.size = 2, line.size = 1)
dev.off()


