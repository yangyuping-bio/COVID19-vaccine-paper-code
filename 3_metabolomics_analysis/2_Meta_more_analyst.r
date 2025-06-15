library(Seurat)
library(dplyr)
library(tidyr)
library(jsonlite)

getwd()

# 读取通路和基因数据
df <- read.csv("data/kegg_metabolism_pathways.csv")
head(df,n=1)

calculate_pathway_activity <- function(expr_matrix, pathways, cell_types) {
  if (ncol(expr_matrix) != length(cell_types)) {
    stop("Number of columns in expr_matrix does not match the length of cell_types.")
  }
  pathway_scores <- list()
  
  for (pathway in pathways$Pathway) {
    #cat("Processing pathway:", pathway, "\n")
    # 提取通路基因
    pathway_genes <- unlist(strsplit(pathways$Genes[pathways$Pathway == pathway], ", "))
    matching_genes <- pathway_genes[pathway_genes %in% rownames(expr_matrix)]
    #cat("Matched genes:", length(matching_genes), "\n")
    
    if (length(matching_genes) == 0) {
      message("No matching genes found for pathway: ", pathway)
      next
    }
    
    # 子集化表达矩阵
    expr_subset <- expr_matrix[matching_genes, , drop = FALSE]
    #cat("Expression subset dimensions:", dim(expr_subset), "\n")
    
    # 按细胞类型计算基因的均值
    indices_by_cell_type <- split(1:ncol(expr_subset), cell_types)
    mean_expr_by_cell <- lapply(indices_by_cell_type, function(indices) {
      rowMeans(expr_subset[, indices, drop = FALSE], na.rm = TRUE)
    })
    # 转换为基因 × 细胞类型矩阵
    mean_expr_by_cell <- do.call(cbind, mean_expr_by_cell)
    # 计算相对表达
    relative_expr <- sweep(mean_expr_by_cell, 1, rowMeans(mean_expr_by_cell, na.rm = TRUE), FUN = "/")
    # 过滤异常值
    upper_bound <- 3 * apply(relative_expr, 1, quantile, probs = 0.9, na.rm = TRUE) #0.75换成0.9
    lower_bound <- (1 / 3) * apply(relative_expr, 1, quantile, probs = 0.1, na.rm = TRUE) #0.25换成0.1
    valid_genes <- sweep(relative_expr, 1, upper_bound, FUN = "<=") &
                   sweep(relative_expr, 1, lower_bound, FUN = ">=")
    relative_expr[!valid_genes] <- NA
    # 计算权重
    weights <- sapply(matching_genes, function(gene) sum(grepl(gene, pathways$Genes)))
    weights <- 1 / weights
    
    # 计算路径活性
    tryCatch({
      pathway_activity <- colSums(relative_expr * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
      pathway_scores[[pathway]] <- pathway_activity
    }, error = function(e) {
      message("Error in pathway activity calculation for pathway: ", pathway, " - ", e$message)
    })
  }
  
  pathway_scores_df <- do.call(cbind, pathway_scores)
  colnames(pathway_scores_df) <- pathways$Pathway
  return(as.data.frame(pathway_scores_df))
}

#无合并数据的操作

calculate_and_plot_old <- function(seurat_obj, 
                               df, 
                               group_var1, 
                               group_var2, 
                               output_prefix, 
                               column_order1, 
                               column_order2, 
                               heatmap_width = 8, 
                               heatmap_height = 12) {
  # 加载必要的包
  library(Seurat)
  library(pheatmap)
  
  # Helper Function: 检查和处理数据
  preprocess_matrix <- function(matrix, column_order) {
    # 替换 NA
    matrix[is.na(matrix)] <- 0
    
    # 去除标准差为 0 的行
    row_sd <- apply(matrix, 1, sd)
    matrix <- matrix[row_sd != 0, ]
    
    # 补全缺失的列
    missing_cols <- setdiff(column_order, colnames(matrix))
    for (col in missing_cols) {
      matrix[, col] <- 0
    }
    matrix <- matrix[, column_order, drop = FALSE]
    
    return(matrix)
  }
  
  # Step 1: 对 group_var1 分簇计算
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")  # 保持稀疏矩阵格式
  cell_types <- seurat_obj[[group_var1]][, 1]  # 如果是 data.frame，提取第一列
  
  if (ncol(expr_matrix) != length(cell_types)) {
    stop("Number of columns in expr_matrix does not match the length of cell_types.")
  }
  
  pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, cell_types)
  write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", group_var1, "_pathway_activity.csv"))
  
  pathway_activity_transposed <- t(pathway_activity_scores)
  pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order1)
  
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE, #FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", group_var1),
    filename = paste0("figures/", output_prefix, "_", group_var1, "_pathway_activity.pdf"),
    width = heatmap_width,
    height = heatmap_height,
    show_colnames = TRUE
  )
  
  # Step 2: 对 group_var2 分簇计算
  conditions <- seurat_obj[[group_var2]][, 1]
  
  if (ncol(expr_matrix) != length(conditions)) {
    stop("Number of columns in expr_matrix does not match the length of conditions.")
  }
  
  pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, conditions)
  write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", group_var2, "_pathway_activity.csv"))
  
  pathway_activity_transposed <- t(pathway_activity_scores)
  pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order2)
  
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE, #FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", group_var2),
    filename = paste0("figures/", output_prefix, "_", group_var2, "_pathway_activity.pdf"),
    width = heatmap_width,
    height = heatmap_height,
    show_colnames = TRUE
  )
  
  # Step 3: 对 group_var1 下每个类别的 group_var2 分簇计算
  categories <- unique(seurat_obj[[group_var1]][, 1])
  
  for (category in categories) {
    subset_obj <- subset(seurat_obj, subset = !!rlang::sym(group_var1) == category)
    expr_matrix <- GetAssayData(subset_obj, assay = "RNA", slot = "data")
    if (ncol(expr_matrix) == 0) next
    
    conditions <- subset_obj[[group_var2]][, 1]
    
    if (ncol(expr_matrix) != length(conditions)) next
    
    pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, conditions)
    pathway_activity_scores[is.na(pathway_activity_scores)] <- 0  # 替换 NA
    
    pathway_activity_transposed <- t(pathway_activity_scores)
    pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order2)
    
    # 保存结果为 CSV
    write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", category, "_", group_var2, "_pathway_activity.csv"))
    
    # 绘制热图
    pheatmap(
      pathway_activity_transposed,
      scale = "row",
      cluster_rows = TRUE, #FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(c("blue", "white", "red"))(50),
      main = paste("Pathway Activity Heatmap -", category),
      filename = paste0("figures/", output_prefix, "_", category, "_", group_var2, "_pathway_activity.pdf"),
      width = heatmap_width,
      height = heatmap_height
    )
  }
}

calculate_and_plot <- function(seurat_obj, 
                               df, 
                               group_var1, 
                               group_var2, 
                               output_prefix, 
                               column_order1, 
                               column_order2, 
                               heatmap_width = 8, 
                               heatmap_height = 12) {
  # 加载必要的包
  library(Seurat)
  library(pheatmap)
  library(dplyr)
  
  # Helper Function: 检查和处理数据
  preprocess_matrix <- function(matrix, column_order) {
    matrix[is.na(matrix)] <- 0
    row_sd <- apply(matrix, 1, sd)
    matrix <- matrix[row_sd != 0, ]
    missing_cols <- setdiff(column_order, colnames(matrix))
    for (col in missing_cols) {
      matrix[, col] <- 0
    }
    matrix <- matrix[, column_order, drop = FALSE]
    return(matrix)
  }
  
  # 合并所有 pathway_activity_scores 的结果
  merge_data <- list()
  
  # Step 1: 对 group_var1 分簇计算
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  cell_types <- seurat_obj[[group_var1]][, 1]
  
  if (ncol(expr_matrix) != length(cell_types)) {
    stop("Number of columns in expr_matrix does not match the length of cell_types.")
  }
  
  pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, cell_types)
  write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", group_var1, "_pathway_activity.csv"))
  
  pathway_activity_transposed <- t(pathway_activity_scores)
  pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order1)
  
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE, 
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", group_var1),
    filename = paste0("figures/", output_prefix, "_", group_var1, "_pathway_activity.pdf"),
    width = heatmap_width,
    height = heatmap_height,
    show_colnames = TRUE
  )
    
  # Step 2: 对 group_var2 分簇计算
  conditions <- seurat_obj[[group_var2]][, 1]
  
  if (ncol(expr_matrix) != length(conditions)) {
    stop("Number of columns in expr_matrix does not match the length of conditions.")
  }
  
  pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, conditions)
  write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", group_var2, "_pathway_activity.csv"))
  
  pathway_activity_transposed <- t(pathway_activity_scores)
  pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order2)
  
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE, #FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", group_var2),
    filename = paste0("figures/", output_prefix, "_", group_var2, "_pathway_activity.pdf"),
    width = heatmap_width,
    height = heatmap_height,
    show_colnames = TRUE
  )
    
  # Step 3: 对 group_var1 下每个类别的 group_var2 分簇计算
  categories <- unique(seurat_obj[[group_var1]][, 1])
  
  for (category in categories) {
    subset_obj <- subset(seurat_obj, subset = !!rlang::sym(group_var1) == category)
    expr_matrix <- GetAssayData(subset_obj, assay = "RNA", slot = "data")
    if (ncol(expr_matrix) == 0) next
    
    conditions <- subset_obj[[group_var2]][, 1]
    if (ncol(expr_matrix) != length(conditions)) next
    
    pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, conditions)
    pathway_activity_scores[is.na(pathway_activity_scores)] <- 0  # 替换 NA
    
    # 添加 cell_type 列
    pathway_activity_scores <- as.data.frame(pathway_activity_scores)
    pathway_activity_scores$cell_type <- category
    
    # 合并到 merge_data
    merge_data[[category]] <- pathway_activity_scores
    
    # 热图
    pathway_activity_transposed <- t(pathway_activity_scores[ , -ncol(pathway_activity_scores)])  # 除去 cell_type 列
    pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order2)
    
    pheatmap(
      pathway_activity_transposed,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = colorRampPalette(c("blue", "white", "red"))(50),
      main = paste("Pathway Activity Heatmap -", category),
      filename = paste0("figures/", output_prefix, "_", category, "_", group_var2, "_pathway_activity.pdf"),
      width = heatmap_width,
      height = heatmap_height
    )
  }
  
  # 合并所有 pathway_activity_scores 并保存
  final_merge_data <- do.call(rbind, merge_data)
  write.csv(final_merge_data, file = paste0("data/", output_prefix, "_merged_pathway_activity.csv"), row.names = TRUE)
}

calculate_and_plot_only_1 <- function(seurat_obj, 
                               df, 
                               group_var1, 
                               group_var2, 
                               output_prefix, 
                               column_order1, 
                               column_order2, 
                               heatmap_width = 8, 
                               heatmap_height = 12) {
  # 加载必要的包
  library(Seurat)
  library(pheatmap)
  library(dplyr)
  
  # Helper Function: 检查和处理数据
  preprocess_matrix <- function(matrix, column_order) {
    matrix[is.na(matrix)] <- 0
    row_sd <- apply(matrix, 1, sd)
    matrix <- matrix[row_sd != 0, ]
    missing_cols <- setdiff(column_order, colnames(matrix))
    for (col in missing_cols) {
      matrix[, col] <- 0
    }
    matrix <- matrix[, column_order, drop = FALSE]
    return(matrix)
  }
  
  # 合并所有 pathway_activity_scores 的结果
  merge_data <- list()
  
  # Step 1: 对 group_var1 分簇计算
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  cell_types <- seurat_obj[[group_var1]][, 1]
  
  if (ncol(expr_matrix) != length(cell_types)) {
    stop("Number of columns in expr_matrix does not match the length of cell_types.")
  }
  
  pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, cell_types)
  write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", group_var1, "_pathway_activity.csv"))
  
  pathway_activity_transposed <- t(pathway_activity_scores)
  pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order1)
  
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE, 
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", group_var1),
    filename = paste0("figures/", output_prefix, "_", group_var1, "_pathway_activity.pdf"),
    width = heatmap_width,
    height = heatmap_height,
    show_colnames = TRUE
  )
  
}

sample.integrated <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/5_new_level_1_2/3_All_Cell_CCA_integrated_by_Sample.rds")

# All
calculate_and_plot(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel1",
  group_var2 = "Condition",
  output_prefix = "1_All_ct1",
  column_order1 = c('CD4T','CD8T','Other_T','NK','Bcell','Monocyte','DC','Platelet','ILC','HPSC'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

# All
calculate_and_plot(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel2",
  group_var2 = "Condition",
  output_prefix = "2_All_ct2",
  column_order1 = c('CD4T','CD8T','MAIT','gdT','Tdn','NK','Bcell','Monocyte','cDC','pDC','Platelet','ILC','HPSC'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

# All
calculate_and_plot(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel3",
  group_var2 = "Condition",
  output_prefix = "3_All_ct3",
  column_order1 = unique(sample.integrated$celltypel3),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

# All
calculate_and_plot_only_1(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel3",
  group_var2 = "Condition",
  output_prefix = "3_All_ct3",
  column_order1 = sort(unique(sample.integrated$celltypel3)),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2'),
   heatmap_width = 10.5, heatmap_height = 12
)

data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3"
Monocyte <- readRDS(file.path(data_save_path, "2_Mono/2_Reumap/Trajectory_Mono.rds"))

# CD4T
calculate_and_plot(
  seurat_obj = Monocyte,
  df = df,
  group_var1 = "ct_level_4",
  group_var2 = "Condition",
  output_prefix = "Monocyte",
  column_order1 = c('Mono_CD14','Mono_CD14_ATG7+','Mono_CD14_CD16','Mono_CD16','Mono_CD16_ATG7+'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)



calculate_and_plot_only12 <- function(seurat_obj, 
                               df, 
                               group_var1, 
                               group_var2, 
                               output_prefix, 
                               column_order1, 
                               column_order2, 
                               heatmap_width = 8, 
                               heatmap_height = 12) {
  
  # 加载必要的包
  library(Seurat)
  library(pheatmap)
  
  # Helper Function: 检查和处理数据
  preprocess_matrix <- function(matrix, column_order) {
    # 替换 NA
    matrix[is.na(matrix)] <- 0
    
    # 去除标准差为 0 的行
    row_sd <- apply(matrix, 1, sd)
    matrix <- matrix[row_sd != 0, ]
    
    # 补全缺失的列
    missing_cols <- setdiff(column_order, colnames(matrix))
    for (col in missing_cols) {
      matrix[, col] <- 0
    }
    matrix <- matrix[, column_order, drop = FALSE]
    
    return(matrix)
  }
  
  # Step 1: 对 group_var1 分簇计算
  expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")  # 保持稀疏矩阵格式
  cell_types <- seurat_obj[[group_var1]][, 1]  # 如果是 data.frame，提取第一列
  
  if (ncol(expr_matrix) != length(cell_types)) {
    stop("Number of columns in expr_matrix does not match the length of cell_types.")
  }
  
  pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, cell_types)
  write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", group_var1, "_pathway_activity.csv"))
  
  pathway_activity_transposed <- t(pathway_activity_scores)
  pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order1)
  
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE, #FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", group_var1),
    filename = paste0("figures/", output_prefix, "_", group_var1, "_pathway_activity.pdf"),
    width = heatmap_width,
    height = heatmap_height,
    show_colnames = TRUE
  )
  
  # Step 2: 对 group_var2 分簇计算
  conditions <- seurat_obj[[group_var2]][, 1]
  
  if (ncol(expr_matrix) != length(conditions)) {
    stop("Number of columns in expr_matrix does not match the length of conditions.")
  }
  
  pathway_activity_scores <- calculate_pathway_activity(as.matrix(expr_matrix), df, conditions)
  write.csv(pathway_activity_scores, file = paste0("data/", output_prefix, "_", group_var2, "_pathway_activity.csv"))
  
  pathway_activity_transposed <- t(pathway_activity_scores)
  pathway_activity_transposed <- preprocess_matrix(pathway_activity_transposed, column_order2)
  
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE, #FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", group_var2),
    filename = paste0("figures/", output_prefix, "_", group_var2, "_pathway_activity.pdf"),
    width = heatmap_width,
    height = heatmap_height,
    show_colnames = TRUE
  )
}

data_save_path <- "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3"
CD4T <- readRDS(file.path(data_save_path, "4_Tcell_CD4T/3_Reumap_with_Treg/CD4T_Reumap_with_Treg.rds"))
CD8T <- readRDS(file.path(data_save_path, "4_Tcell_CD8T/2_Reumap/Trajectory_CD8T.rds"))
NK <- readRDS(file.path(data_save_path, "4_Tcell_NK/1_CCA/Trajectory_NK.rds"))
Bcell<- readRDS(file.path(data_save_path, "2_Bcell/2_Reumap_with_Plasma/Reumap_Bcell_with_Plasma.rds"))
Monocyte <- readRDS(file.path(data_save_path, "2_Mono/2_Reumap/Trajectory_Mono.rds"))


# Bcell
calculate_and_plot(
  seurat_obj = Bcell,
  df = df,
  group_var1 = "ct_level_4",
  group_var2 = "Condition",
  output_prefix = "Bcell",
  column_order1 = c('Bmn_IGHD', 'Bm_CD27+_IGHM+', 'Bm_CD27+_IGHM+_SOX5', 'Bm_CD27+_IGHM-', 'Plasma'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

unique(CD4T$ct_level_4)

# CD4T
calculate_and_plot(
  seurat_obj = CD4T,
  df = df,
  group_var1 = "ct_level_4",
  group_var2 = "Condition",
  output_prefix = "CD4T",
  column_order1 = c('CD4_Tn','CD4_Tcm','CD4_T_CCR6','CD4_Th2','CD4_Treg','CD4_Tem','CD4_Temra'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

unique(CD8T$ct_level_4)

# CD8T
calculate_and_plot(
  seurat_obj = CD8T,
  df = df,
  group_var1 = "ct_level_4",
  group_var2 = "Condition",
  output_prefix = "CD8T",
  column_order1 = c('CD8_Tn','CD8_NELL2','CD8_Tem','CD8_Temra'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

unique(Monocyte$ct_level_4)

# CD4T
calculate_and_plot(
  seurat_obj = Monocyte,
  df = df,
  group_var1 = "ct_level_4",
  group_var2 = "Condition",
  output_prefix = "Monocyte",
  column_order1 = c('Mono_CD14','Mono_CD14_ATG7+','Mono_CD14_CD16','Mono_CD16','Mono_CD16_ATG7+'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

unique(NK$ct_level_3)

# NK
calculate_and_plot(
  seurat_obj = NK,
  df = df,
  group_var1 = "ct_level_3",
  group_var2 = "Condition",
  output_prefix = "NK",
  column_order1 = c('NK_cyto_KLRC2','NK_cyto_FCER1G','NK_rest'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)



Platelet <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/1_Platelet/1_Reumap/Platelet_Reumap.rds")

table(Platelet$Condition)

# 对Platelet,Conidtion分簇计算
expr_matrix <- as.matrix(GetAssayData(Platelet, assay = "RNA", slot = "data"))
pathway_activity_scores <- calculate_pathway_activity(expr_matrix, df, Platelet$Condition)

write.csv(pathway_activity_scores,'data/1_Platelet_Condition_pathway_activity.csv')

# 转置 pathway 活性得分矩阵
pathway_activity_transposed <- t(pathway_activity_scores)

# 替换 NA 和 NaN 为 0，Inf 替换为有限值
pathway_activity_transposed[is.na(pathway_activity_transposed)] <- 0
pathway_activity_transposed[is.nan(pathway_activity_transposed)] <- 0
pathway_activity_transposed[is.infinite(pathway_activity_transposed)] <- 0

# 检查标准差，删除全为 0 或无变化的行
row_sd <- apply(pathway_activity_transposed, 1, sd, na.rm = TRUE)
pathway_activity_transposed <- pathway_activity_transposed[row_sd != 0, , drop = FALSE]

# 确保路径活性矩阵不为空
if (nrow(pathway_activity_transposed) == 0) {
  stop("No valid rows left after filtering. Please check input data.")
}

# 加载必要的包
library(pheatmap)
# 假设 pathway_activity_transposed 的列名是通路名称
column_order <- c('A0','A1','A2','B0','B1','B2','C0','C1','C2')

# 如果矩阵只有一行，禁用行聚类
cluster_rows <- nrow(pathway_activity_transposed) > 1

pheatmap(
  pathway_activity_transposed,
  scale = "row",
  cluster_rows = cluster_rows,#FALSE,  # 仅多行时启用行聚类
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Pathway Activity Heatmap",
  filename = "figures/1_Platelet_Condition_pathway_activity.pdf",
  width = 8,
  height = 12,
    column_order = column_order
)



sample.integrated <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/5_new_level_1_2/3_All_Cell_CCA_integrated_by_Sample.rds")

unique(sample.integrated$celltypel2)

# All
calculate_and_plot(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel2",
  group_var2 = "Condition",
  output_prefix = "2_All_ct2",
  column_order1 = c('CD4T','CD8T','MAIT','gdT','Tdn','NK','Bcell','Monocyte','cDC','pDC','Platelet','ILC','HPSC'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

unique(sample.integrated$celltypel1)

# All
calculate_and_plot(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel1",
  group_var2 = "Condition",
  output_prefix = "1_All_ct1",
  column_order1 = c('CD4T','CD8T','Other_T','NK','Bcell','Monocyte','DC','Platelet','ILC','HPSC'),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

unique(sample.integrated$celltypel3)

# All
calculate_and_plot(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel3",
  group_var2 = "Condition",
  output_prefix = "3_All_ct3",
  column_order1 = unique(sample.integrated$celltypel3),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')
)

sort(unique(sample.integrated$celltypel3))

# All
calculate_and_plot_only12(
  seurat_obj = sample.integrated ,
  df = df,
  group_var1 = "celltypel3",
  group_var2 = "Condition",
  output_prefix = "3_All_ct3",
  column_order1 = sort(unique(sample.integrated$celltypel3)),
  column_order2 = c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2'),
   heatmap_width = 11, 
   heatmap_height = 12
    
)

calculate_pathway_activity <- function(expr_matrix, pathways, cell_types) {
  if (ncol(expr_matrix) != length(cell_types)) {
    stop("Number of columns in expr_matrix does not match the length of cell_types.")
  }
  
  pathway_scores <- list()
  
  for (pathway in pathways$Pathway) {
    #cat("Processing pathway:", pathway, "\n")
    
    # 提取通路基因
    pathway_genes <- unlist(strsplit(pathways$Genes[pathways$Pathway == pathway], ", "))
    matching_genes <- pathway_genes[pathway_genes %in% rownames(expr_matrix)]
    #cat("Matched genes:", length(matching_genes), "\n")
    
    if (length(matching_genes) == 0) {
      message("No matching genes found for pathway: ", pathway)
      next
    }
    
    # 子集化表达矩阵
    expr_subset <- expr_matrix[matching_genes, , drop = FALSE]
    #cat("Expression subset dimensions:", dim(expr_subset), "\n")
    
    # 按细胞类型计算基因的均值
    indices_by_cell_type <- split(1:ncol(expr_subset), cell_types)
    mean_expr_by_cell <- lapply(indices_by_cell_type, function(indices) {
      rowMeans(expr_subset[, indices, drop = FALSE], na.rm = TRUE)
    })
    
    # 转换为基因 × 细胞类型矩阵
    mean_expr_by_cell <- do.call(cbind, mean_expr_by_cell)
    
    # 计算相对表达
    relative_expr <- sweep(mean_expr_by_cell, 1, rowMeans(mean_expr_by_cell, na.rm = TRUE), FUN = "/")
    
    # 过滤异常值
    upper_bound <- 3 * apply(relative_expr, 1, quantile, probs = 0.75, na.rm = TRUE)
    lower_bound <- (1 / 3) * apply(relative_expr, 1, quantile, probs = 0.25, na.rm = TRUE)
    valid_genes <- sweep(relative_expr, 1, upper_bound, FUN = "<=") &
                   sweep(relative_expr, 1, lower_bound, FUN = ">=")
    relative_expr[!valid_genes] <- NA
    
    # 计算权重
    weights <- sapply(matching_genes, function(gene) sum(grepl(gene, pathways$Genes)))
    weights <- 1 / weights
    
    # 计算路径活性
    tryCatch({
      pathway_activity <- colSums(relative_expr * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
      pathway_scores[[pathway]] <- pathway_activity
    }, error = function(e) {
      message("Error in pathway activity calculation for pathway: ", pathway, " - ", e$message)
    })
  }
  
  pathway_scores_df <- do.call(cbind, pathway_scores)
  colnames(pathway_scores_df) <- pathways$Pathway
  return(as.data.frame(pathway_scores_df))
}

#对celltypes分簇计算
expr_matrix <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))
pathway_activity_scores <- calculate_pathway_activity(expr_matrix, df, seurat_obj$ct_level_4)

# 转置数据框，使行代表细胞类型，列代表通路
pathway_activity_transposed <- t(pathway_activity_scores)
write.csv(pathway_activity_scores,'data/1_Bcell_ct_level_4_pathway_activity.csv')
# 加载必要的包
library(pheatmap)
# 假设 pathway_activity_transposed 的列名是通路名称
column_order <- c('Bmn_IGHD','Bm_CD27+_IGHM+','Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-','Plasma')
pheatmap(pathway_activity_transposed,
         scale = "row",             # 按行标准化
         cluster_rows = TRUE,       # 对行（细胞类型）聚类
         cluster_cols = FALSE,       # 对列（通路）聚类
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Pathway Activity Heatmap",
        filename='figures/1_Bcell_ct_level_4_pathway_activity.pdf',
         width = 10, height = 15,
        show_colnames = TRUE,     # 显示列名
         column_order = column_order)


# 对Conidtion分簇计算
pathway_activity_scores <- calculate_pathway_activity(expr_matrix, df, seurat_obj$Condition)

write.csv(pathway_activity_scores,'data/1_Bcell_Condition_pathway_activity.csv')

# 转置数据框，使行代表细胞类型，列代表通路
pathway_activity_transposed <- t(pathway_activity_scores)

# 加载必要的包
library(pheatmap)
# 假设 pathway_activity_transposed 的列名是通路名称
column_order <- c('A0','A1','A2','B0','B1','B2','C0','C1','C2')
pheatmap(pathway_activity_transposed,
         scale = "row",             # 按行标准化
         cluster_rows = TRUE,       # 对行（细胞类型）聚类
         cluster_cols = FALSE,       # 对列（通路）聚类
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Pathway Activity Heatmap",
        filename='figures/1_Bcell_Condition_pathway_activity.pdf',
         width = 8, height = 12,
        show_colnames = TRUE,     # 显示列名
         column_order = column_order)


#对每个celltype下每个Conidtion分簇计算
for (category in categories) {
  # 子集化数据
  subset_obj <- subset(seurat_obj, subset = ct_level_4 == category)
  expr_matrix <- as.matrix(GetAssayData(subset_obj, assay = "RNA", slot = "data"))
  if (nrow(expr_matrix) == 0) next # 跳过空数据

  # 计算路径活动评分
  pathway_activity_scores <- calculate_pathway_activity(expr_matrix, df, subset_obj$Condition)
  
  # 转置前检查并处理 NA
  pathway_activity_scores[is.na(pathway_activity_scores)] <- 0
  
  # 转置并检查标准化问题
  pathway_activity_transposed <- t(pathway_activity_scores)
  row_sd <- apply(pathway_activity_transposed, 1, sd)
  pathway_activity_transposed <- pathway_activity_transposed[row_sd != 0, ]
  
  # 确保列名完整
  missing_cols <- setdiff(column_order, colnames(pathway_activity_transposed))
  for (col in missing_cols) {
    pathway_activity_transposed[, col] <- 0
  }
  pathway_activity_transposed <- pathway_activity_transposed[, column_order]

  # 绘图
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", category),
    filename = paste0("figures/", category, "_Condition_pathway_activity.pdf"),
    width = 8, height = 12
  )
}

for (category in categories) {
  # 子集化数据
  subset_obj <- subset(seurat_obj, subset = ct_level_4 == category)
  expr_matrix <- as.matrix(GetAssayData(subset_obj, assay = "RNA", slot = "data"))
  if (nrow(expr_matrix) == 0) next # 跳过空数据

  # 计算路径活动评分
  pathway_activity_scores <- calculate_pathway_activity(expr_matrix, df, subset_obj$Condition)
  
  # 转置前检查并处理 NA
  pathway_activity_scores[is.na(pathway_activity_scores)] <- 0
  
  # 转置并检查标准化问题
  pathway_activity_transposed <- t(pathway_activity_scores)
  row_sd <- apply(pathway_activity_transposed, 1, sd)
  pathway_activity_transposed <- pathway_activity_transposed[row_sd != 0, ]
  
  # 确保列名完整
  missing_cols <- setdiff(column_order, colnames(pathway_activity_transposed))
  for (col in missing_cols) {
    pathway_activity_transposed[, col] <- 0
  }
  pathway_activity_transposed <- pathway_activity_transposed[, column_order]

  # 绘图
  pheatmap(
    pathway_activity_transposed,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", category),
    filename = paste0("figures/", category, "_Condition_pathway_activity.pdf"),
    width = 8, height = 12
  )
}

# 加载必要的包
library(Seurat)
library(pheatmap)

# 获取所有唯一的 ct_level_4 值
categories <- unique(seurat_obj$ct_level_4)

# 遍历每个类别进行分析
for (category in categories) {
  # 1. 子集化 Seurat 对象
  subset_obj <- subset(seurat_obj, subset = ct_level_4 == category)
  
  # 2. 提取表达矩阵
  expr_matrix <- as.matrix(GetAssayData(subset_obj, assay = "RNA", slot = "data"))
  
  # 3. 计算 pathway activity scores
  pathway_activity_scores <- calculate_pathway_activity(expr_matrix, df, subset_obj$Condition)
  
  # 4. 保存 pathway activity 数据
  output_csv <- paste0("data/", category, "_Condition_pathway_activity.csv")
  write.csv(pathway_activity_scores, output_csv)
  
  # 5. 转置数据以供绘图
  pathway_activity_transposed <- t(pathway_activity_scores)
  
  # 6. 设置列顺序 (根据具体需求修改)
  column_order <- c('A0','A1','A2','B0','B1','B2','C0','C1','C2')
  
  # 7. 绘制热图并保存
  output_pdf <- paste0("figures/", category, "_Condition_pathway_activity.pdf")
  pheatmap(
    pathway_activity_transposed,
    scale = "row",             # 按行标准化
    cluster_rows = TRUE,       # 对行（细胞类型）聚类
    cluster_cols = FALSE,      # 不对列（通路）聚类
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste("Pathway Activity Heatmap -", category),
    filename = output_pdf,     # 输出文件名
    width = 10, height = 15,
    show_colnames = TRUE,      # 显示列名
    column_order = column_order
  )
}

unique(seurat_obj$ct_level_4)

# 转换数据框为长格式以便于 ggplot2 可视化
library(reshape2)
library(ggplot2)

pathway_long <- melt(pathway_activity_scores, 
                     variable.name = "Pathway", 
                     value.name = "Activity")
pdf('figures/2_boxplot_Bcell.pdf',35,10)
# 绘制箱线图
ggplot(pathway_long, aes(x = Pathway, y = Activity, fill = Pathway)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 8) +
  labs(title = "Pathway Activity Distribution by Cell Types",
       x = "Pathway",
       y = "Activity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 选择特定通路
selected_pathway <- "Glycolysis / Gluconeogenesis"

# 提取对应数据
barplot_data <- pathway_activity_scores[, selected_pathway, drop = FALSE]
barplot_data$CellType <- rownames(barplot_data)
pdf('figures/3_barplot_Bcell.pdf',5,5)
# 绘制条形图
ggplot(barplot_data, aes(x = CellType, y = get(selected_pathway), fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(title = paste("Activity of", selected_pathway, "Pathway"),
       x = "Cell Type",
       y = "Activity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 绘制小提琴图
pdf('figures/4_violin_Bcell.pdf',35,10)
ggplot(pathway_long, aes(x = Pathway, y = Activity, fill = Pathway)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.color = "red") +
  labs(title = "Violin Plot of Pathway Activity by Cell Types",
       x = "Pathway",
       y = "Activity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()










