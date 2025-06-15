library(dplyr)
library(tidyr)

dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/8_Others/2_Metabolism/8_More_Meta_analyst/2_Pathway_activity_analysis'

Pathway_activity_res <- read.csv(paste0(dir_path,'/data/',"3_All_ct3_merged_pathway_activity.csv"),check.names = FALSE,row.names = 1)
colnames(Pathway_activity_res)

unique(Pathway_activity_res$cell_type)



celltype <- 'CD4_Temra'
selected_res <- Pathway_activity_res %>% 
  filter(cell_type == celltype) 
rownames(selected_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
selected_res$cell_type <- NULL

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    min_values <- sort(., decreasing = FALSE, na.last = TRUE)[1:3]
    # 检查第1行是否是该列最小值且不等于0
    Con1 <- (.[1] %in% min_values && .[1] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:4]
    Con2 <- !.[4] %in% top_values
    Con3 <- .[8] %in% top_values
    Con4 <- .[9] %in% top_values
    
    # 满足所有条件
    Con1  && Con3 && Con4 && Con2  
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = paste0('figures/1_',celltype,'_pheatmap.pdf'),
        width = 6.5, height = 3.5)      # Y轴标签

celltype <- 'CD4_Treg'
selected_res <- Pathway_activity_res %>% 
  filter(cell_type == celltype) 
rownames(selected_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
selected_res$cell_type <- NULL

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    min_values <- sort(., decreasing = FALSE, na.last = TRUE)[1:3]
    # 检查第1行是否是该列最小值且不等于0
    Con1 <- (.[1] %in% min_values && .[1] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:1]
    top_values_2 <- sort(., decreasing = TRUE, na.last = TRUE)[1:4]
    Con2 <- .[2] %in% top_values
    #Con3 <- .[7] %in% top_values_2
    #Con4 <- .[6] %in% top_values_2
    
    # 满足所有条件
    Con1  && Con2 #&& Con3  && Con4  
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = paste0('figures/1_',celltype,'_pheatmap.pdf'),
        width = 6.5, height = 3.5)      # Y轴标签

celltype <- 'CD4_Treg'
selected_res <- Pathway_activity_res %>% 
  filter(cell_type == celltype) 
rownames(selected_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
selected_res$cell_type <- NULL

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    min_values <- sort(., decreasing = FALSE, na.last = TRUE)[1:3]
    # 检查第1行是否是该列最小值且不等于0
    Con1 <- (.[1] %in% min_values && .[1] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:1]
    top_values_2 <- sort(., decreasing = TRUE, na.last = TRUE)[1:3]
    Con2 <- .[2] %in% top_values
    Con3 <- .[7] %in% top_values_2
    #Con4 <- .[6] %in% top_values_2
    
    # 满足所有条件
    Con1  && Con2 && Con3  #&& Con4  
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = paste0('figures/1_',celltype,'_pheatmap_less.pdf'),
        width = 6.5, height = 3.5)      # Y轴标签





Mono_res <- Pathway_activity_res %>% 
  filter(cell_type == 'Mono_CD14_ATG7+') 
rownames(Mono_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")

# 筛选满足条件的列
selected_columns <- Mono_res %>%
  select(where(~ {
    # 检查第二行是否是该列最大值
    is_second_row_max <- .[2] == max(., na.rm = TRUE)
    
    # 检查第四行和第七行的值是否不在列的前3/2大值中
    top_3_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:3]
    top_2_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:2]
    is_fourth_row_not_top3 <- !.[4] %in% top_3_values
    is_seventh_row_not_top3 <- !.[7] %in% top_2_values
    
    # 满足所有条件
    is_second_row_max && is_fourth_row_not_top3 && is_seventh_row_not_top3
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = "Mono_CD14_ATG7", # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = 'figures/1_Mono_CD14_ATG7_pheatmap.pdf',
        width = 6.5, height = 3.5)      # Y轴标签





celltype <- 'CD8_Temra'
selected_res <- Pathway_activity_res %>% 
  filter(cell_type == celltype) 
rownames(selected_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
selected_res$cell_type <- NULL

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    # 检查第1行是否是该列最小值且不等于0
    Con1 <- (.[1] == min(., na.rm = TRUE) && .[1] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_3_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:3]
    #Con2 <- .[7] %in% top_3_values
    Con3 <- .[8] %in% top_3_values
    Con4 <- .[9] %in% top_3_values
    
    # 满足所有条件
    Con1  && Con3 && Con4 #&& Con2
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = 'figures/1_CD8_Temra_pheatmap_more.pdf',
        width = 6.5, height = 3.5)      # Y轴标签

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    # 检查第1行是否是该列最小值且不等于0
    Con1 <- (.[1] == min(., na.rm = TRUE) && .[1] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_3_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:3]
    Con2 <- .[7] %in% top_3_values
    Con3 <- .[8] %in% top_3_values
    Con4 <- .[9] %in% top_3_values
    
    # 满足所有条件
    Con1  && Con3 && Con4 && Con2
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = 'figures/1_CD8_Temra_pheatmap.pdf',
        width = 6.5, height = 3.5)      # Y轴标签

celltype <- 'CD8_NELL2'
selected_res <- Pathway_activity_res %>% 
  filter(cell_type == celltype) 
rownames(selected_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
selected_res$cell_type <- NULL

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    # 检查第2行是否是该列最大值且不等于0
    Con1 <- (.[2] == max(., na.rm = TRUE) && .[1] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_3_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:3]
    Con2 <- !.[1] %in% top_3_values
    # Con3 <- .[8] %in% top_3_values
    # Con4 <- .[9] %in% top_3_values
    
    # 满足所有条件
    Con1  && Con2 #&& Con3 && Con4
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = paste0('figures/1_',celltype,'_pheatmap.pdf'),
        width = 6.5, height = 3.5)      # Y轴标签

celltype <- 'Plasma'
selected_res <- Pathway_activity_res %>% 
  filter(cell_type == celltype) 
rownames(selected_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
selected_res$cell_type <- NULL

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    max_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:1]
    # 检查第1行是否是该列最大值且不等于0
    Con1 <- (.[2] %in% max_values && .[2] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:6]
    #top_values_2 <- sort(., decreasing = TRUE, na.last = TRUE)[1:3]
    Con2 <- .[5] %in% top_values
    #Con3 <- .[7] %in% top_values_2
    #Con4 <- .[6] %in% top_values_2
    
    # 满足所有条件
    Con1  #&& Con2 #&& Con3  #&& Con4  
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = paste0('figures/1_',celltype,'_pheatmap.pdf'),
        width = 6.5, height = 3.5)      # Y轴标签


celltype <- 'Platelet'
selected_res <- Pathway_activity_res %>% 
  filter(cell_type == celltype) 
rownames(selected_res) <-  c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
selected_res$cell_type <- NULL

# 筛选满足条件的列
selected_columns <-selected_res %>%
  select(where(~ {
    max_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:1]
    # 检查第1行是否是该列最大值且不等于0
    Con1 <- (.[5] %in% max_values && .[5] != 0)
      
    # 检查第7行到第9行的值是否在列的前3/2大值中
    top_values <- sort(., decreasing = TRUE, na.last = TRUE)[1:4]
    #top_values_2 <- sort(., decreasing = TRUE, na.last = TRUE)[1:3]
    Con2 <- .[4] %in% top_values
    #Con3 <- .[7] %in% top_values_2
    #Con4 <- .[6] %in% top_values_2
    
    # 满足所有条件
    Con1  #&& Con2 #&& Con3  #&& Con4  
  }))

# 输出结果
colnames(selected_columns)

# 加载必要的包
library(pheatmap)
library(RColorBrewer)
heatmap_data <- t(selected_columns)
# 创建热图
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
         border_color = "white",
         cellwidth = 10, 
         cellheight = 10,
         cluster_rows = TRUE,   # 不聚类行
         cluster_cols = FALSE,   # 不聚类列
         scale = "row",         # 数据已经标准化
         display_numbers = FALSE, # 显示具体数值
         main = celltype, # 热图标题
         xlab = "Conidtion",         # X轴标签
         ylab = "Pathway",
        filename = paste0('figures/1_',celltype,'_pheatmap.pdf'),
        width = 6.5, height = 3.5)      # Y轴标签


dev.list()
dev.off()  # 关闭当前设备
while (!is.null(dev.list())) dev.off()

Pathway_activity_res <- read.csv(paste0(dir_path,'/data/',"2_All_ct2_merged_pathway_activity.csv"),check.names = FALSE,row.names = 1)
#colnames(Pathway_activity_res)

library(dplyr)
library(tibble)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
generate_heatmap <- function(data, pathway_name, output_file,
                             heatmap_title = "Heatmap", 
                             scale_rows = TRUE) {
  # 检查必要的包是否加载
  if (!requireNamespace("tibble", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("pheatmap", quietly = TRUE) ||
      !requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please ensure all required libraries are installed: tibble, tidyr, dplyr, pheatmap, RColorBrewer")
  }
  
  # 数据处理
  selected_res <- data %>%
    select({{ pathway_name }}) %>%
    rownames_to_column(var = "rownames") %>%
    mutate(
      Condition = sub(".*\\.", "", rownames),  # 提取 '.' 后部分为 Condition
      celltype = sub("\\..*", "", rownames)   # 提取 '.' 前部分为 celltype
    ) %>%
    column_to_rownames(var = "rownames")
  
  # 转换为宽表
  wide_table <- selected_res %>%
    pivot_wider(
      names_from = Condition,
      values_from = {{ pathway_name }}
    )
  
  # 设置行名并移除 celltype 列
  wide_table <- data.frame(wide_table)
  rownames(wide_table) <- wide_table$celltype
  wide_table$celltype <- NULL
  
  # 设置 scale 参数
  scale_option <- if (scale_rows) "row" else "none"
  
  # 创建热图
  pheatmap(
    wide_table,
    color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
    border_color = "white",
    cellwidth = 10,
    cellheight = 10,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = scale_option,         # 动态控制是否标准化
    display_numbers = FALSE,
    main = heatmap_title,
    xlab = "Condition",
    ylab = "Pathway",
    filename = output_file,
    width = 5,
    height = 5
  )
}

# 调用函数示例
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Primary bile acid biosynthesis`,
  output_file = "figures/2_Bile_acid_pheatmap_celltypel2.pdf",
  heatmap_title = "Primary bile acid biosynthesis",
    scale_rows = FALSE #TRUE
)

# 调用函数示例
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Tryptophan metabolism`,
  output_file = "figures/2_Tryptophan_metabolism_pheatmap_celltypel2.pdf",
  heatmap_title = "Tryptophan metabolism",
    scale_rows = FALSE #TRUE
)

# 调用函数示例
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Phenylalanine, tyrosine and tryptophan biosynthesis`,
  output_file = "figures/2_Tryptophan_metabolism_2_pheatmap_celltypel2.pdf",
  heatmap_title = "Phenylalanine, tyrosine and tryptophan biosynthesis",
    scale_rows = FALSE #TRUE
)

Phenylalanine, tyrosine and tryptophan biosynthesis

generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Arginine biosynthesis`,
  output_file = "figures/2_Arginine_biosynthesis_pheatmap_celltypel2.pdf",
  heatmap_title = "Arginine biosynthesis",
    scale_rows = FALSE #TRUE
)
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Arginine biosynthesis`,
  output_file = 'figures/2_Arginine_biosynthesis_pheatmap_celltypel2_nom_row.pdf',
  heatmap_title = "Arginine biosynthesis",
    scale_rows = TRUE
)


generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Arginine and proline metabolism`,
  output_file = "figures/2_Arginine_and_proline_pheatmap_celltypel2.pdf",
  heatmap_title = "Arginine and proline metabolism",
    scale_rows = FALSE #TRUE
)
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Arginine and proline metabolism`,
  output_file = 'figures/2_Arginine_and_proline_pheatmap_celltypel2_nom_row.pdf',
  heatmap_title = "Arginine and proline metabolism",
    scale_rows = TRUE
)




generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Sphingolipid metabolism`,
  output_file = "figures/2_Sphingolipid_metabolism_pheatmap_celltypel2.pdf",
  heatmap_title = "Sphingolipid metabolism",
    scale_rows = FALSE #TRUE
)
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Sphingolipid metabolism`,
  output_file = 'figures/2_Sphingolipid_metabolism_pheatmap_celltypel2_nom_row.pdf',
  heatmap_title = "Sphingolipid metabolism",
    scale_rows = TRUE
)


generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Glycosphingolipid biosynthesis - ganglio series`,
  output_file = "figures/2_Glycosphingolipid_biosynthesis_pheatmap_celltypel2.pdf",
  heatmap_title = "Glycosphingolipid biosynthesis - ganglio series",
    scale_rows = FALSE #TRUE
)
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Glycosphingolipid biosynthesis - ganglio series`,
  output_file = "figures/2_Glycosphingolipid_biosynthesis_pheatmap_celltypel2_nom_row.pdf",
  heatmap_title = "Glycosphingolipid biosynthesis - ganglio series",
    scale_rows = TRUE
)


generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Glycosphingolipid biosynthesis - globo and isoglobo series`,
  output_file = "figures/2_Glycosphingolipid_biosynthesis_2_pheatmap_celltypel2.pdf",
  heatmap_title = "Glycosphingolipid biosynthesis - globo and isoglobo series",
    scale_rows = FALSE #TRUE
)
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Glycosphingolipid biosynthesis - globo and isoglobo series`,
  output_file = "figures/2_Glycosphingolipid_biosynthesis_2_pheatmap_celltypel2_nom_row.pdf",
  heatmap_title = "Glycosphingolipid biosynthesis - globo and isoglobo series",
    scale_rows = TRUE
)



generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Steroid biosynthesis`,
  output_file = "figures/2_Steroid_biosynthesis_pheatmap_celltypel2.pdf",
  heatmap_title = "Steroid biosynthesis",
    scale_rows = FALSE #TRUE
)
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Steroid biosynthesis`,
  output_file = 'figures/2_Steroid_biosynthesis_pheatmap_celltypel2_nom_row.pdf',
  heatmap_title = "Steroid biosynthesis",
    scale_rows = TRUE
)


generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Steroid hormone biosynthesis`,
  output_file = "figures/2_Steroid_hormone_biosynthesis_pheatmap_celltypel2.pdf",
  heatmap_title = "Steroid hormone biosynthesis",
    scale_rows = FALSE #TRUE
)
generate_heatmap(
  data = Pathway_activity_res,
  pathway_name = `Steroid hormone biosynthesis`,
  output_file = 'figures/2_Steroid_hormone_biosynthesis_pheatmap_celltypel2_nom_row.pdf',
  heatmap_title = "Steroid hormone biosynthesis",
    scale_rows = TRUE
)










