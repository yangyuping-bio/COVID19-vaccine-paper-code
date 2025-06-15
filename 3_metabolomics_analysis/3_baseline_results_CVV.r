library(dplyr)
library(tidyr)

dir_path <- '/share/home/qlab/projects/project_cvv/yyp_results_3/8_Others/2_Metabolism/8_More_Meta_analyst/2_Pathway_activity_analysis'

Pathway_activity_res <- read.csv(paste0(dir_path,'/data/',"2_All_ct2_merged_pathway_activity.csv"),check.names = FALSE,row.names = 1)
#colnames(Pathway_activity_res)

unique(Pathway_activity_res$cell_type)

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
  
  # 定义行顺序
  row_order <- c('CD4T', 'CD8T', 'Bcell', 'Monocyte', 'NK', 'ILC',
                     'MAIT', 'gdT', 'cDC', 'Tdn', 'pDC', 'HPSC')
  
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
  
  # 按照指定顺序对行进行排序
  if (!is.null(row_order)) {
    wide_table <- wide_table[row_order, , drop = FALSE]
  }
  
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


library(dplyr)
library(tibble)
library(pheatmap)
library(tidyr)
library(RColorBrewer)

generate_heatmap_baseline <- function(data, pathway_name, output_file,pdf_file=NA,
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
  
  # 定义行顺序
  row_order <- c('CD4T', 'CD8T', 'Bcell', 'Monocyte', 'NK', 'ILC',
                'MAIT', 'gdT', 'cDC', 'Tdn', 'pDC','Platelet', 'HPSC')
  
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
  
   # 计算相对于 A0 的 fold change，保留正负趋势
  # wide_table <- wide_table %>%
  # mutate(across(-celltype, ~ (. - A0) / abs(A0 + 1e-6)))  # 添加平滑值避免除以零
    # 计算相对于 A0 的 fold change，并映射到 [-1, 1]
   wide_table <- wide_table %>%
      mutate(across(-celltype, ~ {
        fc <- (. - A0) / abs(A0 + 1e-6)  # 计算原始 fold change
        log_fc <- sign(fc) * log1p(abs(fc))  # 对数变换保留正负号
        max_abs <- max(abs(log_fc), na.rm = TRUE)
        if (max_abs == 0) log_fc else log_fc / max_abs  # 标准化到 [-1, 1]
      }))
      
  # 设置行名并移除 celltype 列
  wide_table <- data.frame(wide_table)
  rownames(wide_table) <- wide_table$celltype
  wide_table$celltype <- NULL
  
  # 按照指定顺序对行进行排序
  if (!is.null(row_order)) {
    wide_table <- wide_table[row_order, , drop = FALSE]
    wide_table <- na.omit(wide_table)  # 移除可能出现的NA行
  }
    
  # 保存处理后的数据到CSV文件
  # 创建data文件夹（如果不存在）
  if (!dir.exists("data")) {
    dir.create("data")
  }
  
  # 生成CSV文件名（基于output_file的名称）
  csv_file <- sub("\\..*$", "", basename(output_file))  # 提取文件名（不含扩展名）
  csv_path <- file.path("data", paste0(csv_file, "_data.csv"))
  
  # 保存wide_table为CSV文件
  write.csv(wide_table, file = csv_path, row.names = TRUE)
    
  # 设置 scale 参数
  scale_option <- if (scale_rows) "row" else "none"
  # 定义legend的breaks，从0到3
  breaks <- seq(-1, 1, length.out = 51)  # 创建0到3的等间隔序列，与颜色数量匹配

  
  # 创建热图
  pheatmap(
    wide_table,
    color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
    breaks = breaks,  # 添加breaks参数控制legend范围
    border_color = "white",
    cellwidth = 10,
    cellheight = 10,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = scale_option,
    display_numbers = FALSE,
    main = paste(heatmap_title, "(fold change over baseline)"),
    xlab = "Condition",
    ylab = "Cell Type",
    filename = pdf_file,
    width = 5,
    height = 5,
    # 添加legend相关参数
     # 添加legend相关参数
    legend = TRUE,
    fontsize = 6,  # 控制整体字体大小
    fontsize_row = 10,  # 行名字体大小
    fontsize_col = 10,  # 列名字体大小
    main_fontsize = 10,  # 主标题字体大小
    legend_fontsize = 10  # legend标题和数字的字体大小
  )
}

pathways <- list(
  list(name = "Primary bile acid biosynthesis", file = "Bile_acid"),
  list(name = "Tryptophan metabolism", file = "Tryptophan_metabolism"),
  list(name = "Phenylalanine, tyrosine and tryptophan biosynthesis", file = "Tryptophan_metabolism_2"),
  list(name = "Arginine biosynthesis", file = "Arginine_biosynthesis"),
  list(name = "Arginine and proline metabolism", file = "Arginine_and_proline"),
  list(name = "Sphingolipid metabolism", file = "Sphingolipid_metabolism"),
  list(name = "Glycosphingolipid biosynthesis - ganglio series", file = "Glycosphingolipid_biosynthesis"),
  list(name = "Glycosphingolipid biosynthesis - globo and isoglobo series", file = "Glycosphingolipid_biosynthesis_2"),
  list(name = "Steroid biosynthesis", file = "Steroid_biosynthesis"),
  list(name = "Steroid hormone biosynthesis", file = "Steroid_hormone_biosynthesis")
)
for (p in pathways) {
  generate_heatmap_baseline(
    data = Pathway_activity_res,
    pathway_name = p$name,
    output_file = sprintf("figures/3_%s_baseline.pdf", p$file),
      pdf_file= sprintf("figures/3_%s_baseline.pdf", p$file),
    heatmap_title = p$name,
    scale_rows = FALSE
  )
}

library(gridExtra)
library(ggplot2)

# 修改函数调用，确保返回绘图对象而不是保存文件
heatmaps <- lapply(pathways, function(p) {
  generate_heatmap_baseline(
    data = Pathway_activity_res,
    pathway_name = p$name,
    output_file = sprintf("figures/3_%s_baseline.pdf", p$file),
    heatmap_title = p$name,
    scale_rows = FALSE
  )
})

# 每个热图保存到新页面
pdf("figures/4_combined_heatmaps_CVV.pdf", width = 6, height = 6, onefile = TRUE)
invisible(lapply(heatmaps, function(h) {
  if (!is.null(h)) {
    print(h)
    plot.new()  # 添加新页面
  }
}))
dev.off()



pathways <- list(
  list(name = "Primary bile acid biosynthesis", file = "Bile_acid"),
  list(name = "Tryptophan metabolism", file = "Tryptophan_metabolism"),
  list(name = "Phenylalanine, tyrosine and tryptophan biosynthesis", file = "Tryptophan_metabolism_2"),
  list(name = "Arginine biosynthesis", file = "Arginine_biosynthesis"),
  list(name = "Arginine and proline metabolism", file = "Arginine_and_proline"),
  list(name = "Sphingolipid metabolism", file = "Sphingolipid_metabolism"),
  list(name = "Glycosphingolipid biosynthesis - ganglio series", file = "Glycosphingolipid_biosynthesis"),
  list(name = "Glycosphingolipid biosynthesis - globo and isoglobo series", file = "Glycosphingolipid_biosynthesis_2"),
  list(name = "Steroid biosynthesis", file = "Steroid_biosynthesis"),
  list(name = "Steroid hormone biosynthesis", file = "Steroid_hormone_biosynthesis")
)

for (p in pathways) {
  generate_heatmap(
    data = Pathway_activity_res,
    pathway_name = p$name,
    output_file = sprintf("figures/2_%s_.pdf", p$file),
    heatmap_title = p$name,
    scale_rows = FALSE
  )
}







plot_Activity_pathway_heatmap <- function(pathway_ct3_merged, Pathway_Activity, margin_size) {
  # 加载必要的包
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  
  # 数据处理部分
  pathway_ct3_merged <- pathway_ct3_merged %>%
    mutate(
      Condition = sub(".*\\.", "", rownames(pathway_ct3_merged))  # 提取 '.' 后部分为 Condition
    )
  
  # 将宽表转换成长表
  pathway_ct3_long <- pathway_ct3_merged %>%
    pivot_longer(
      cols = -c(cell_type, Condition),  # 保留 CellType 和 Condition 列
      names_to = "Pathway",            # 将原列名转换为 Pathway
      values_to = "Value"              # 将值存储到 Value 列
    )
  
  pathway_ct3_long <- pathway_ct3_long %>%
    filter(Pathway %in% Pathway_Activity)
  
  # 自定义绘图函数，每个 Pathway 动态设置颜色映射范围
  ggplot_pathway <- function(df, pathway_name) {
    max_val <- max(df$Value, na.rm = TRUE)  # 获取 Pathway 最大值
    min_val <- min(df$Value, na.rm = TRUE)  # 获取 Pathway 最小值
    
    ggplot(df, aes(x = Condition, y = cell_type, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(
        low = "white", 
        high = "red",
        limits = c(min_val, max_val)  # 根据 Pathway 数据动态设置颜色范围
      ) +
      theme_minimal() +
      labs(
        title = paste(pathway_name),
        x = "Condition",
        y = "Cell Type",
        fill = "Activity Value"
      ) +
      theme(
        strip.text = element_text(size = 10),
        plot.margin = unit(margin_size, "char")
      )
  }
  
  # 按 Pathway 分组，生成各自的热图
  pathway_plots <- pathway_ct3_long %>%
    group_by(Pathway) %>%
    group_split() %>%
    lapply(function(df) {
      pathway_name <- unique(df$Pathway)
      ggplot_pathway(df, pathway_name)
    })
  
  # 返回所有生成的热图对象
  return(pathway_plots)
}

#Pathway_activity_res <- read.csv("data/1_All_ct_merge_merged_pathway_activity.csv",check.names = FALSE,row.names = 1)
colnames(Pathway_activity_res)
pathways <- colnames(Pathway_activity_res)[1:85]
length(pathways)
pathways

library(cowplot)  # 用于拼接多图
# 调用函数绘制Activity_pathway热图
pdf("figures/4_ALL_pathway_activity_heatmap.pdf", width = 40, height = 66)
pathway_heatmaps <-plot_Activity_pathway_heatmap(Pathway_activity_res,pathways,margin_size=c(9,3.5,9,3.5))
plot_grid(plotlist = pathway_heatmaps, ncol = 8) 
dev.off()








