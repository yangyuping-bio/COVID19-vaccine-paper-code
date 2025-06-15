#加载包
library(ggseqlogo)
library(ggplot2)
library(gridExtra)
library(cluster) #fanny
library(tidyverse)
library(reshape2)

getwd()

library(dplyr)
library(reshape2)

# 定义一个函数来计算比值和分类
process_sequences <- function(data, condition_col, cdr3_col, cluster_col, ct_level_col) {
  # 创建新的条件列
  data$Condition_2 <- data[[condition_col]]
  data$Condition_2[data[[condition_col]] %in% c("A0", "A1", "A2", "B0")] <- "before_B1"
  data$Condition_2[data[[condition_col]] %in% c("C0", "C1", "C2", "B1", "B2")] <- "after_B0"
  
  # 计算各个比值
  processed_data <- data %>%
     group_by(cdr3_b_aa) %>%
      mutate(cdr3_count_all = n()) %>%
      ungroup()%>% 
      group_by(Condition, cdr3_b_aa) %>%
      mutate(cdr3_count_condition = n()) %>%
      ungroup()%>% 
      group_by(cluster_id,Condition) %>%
      mutate(total_count_condition_cluster = n()) %>%
      ungroup()%>% 
      group_by(Condition) %>%
      mutate(total_cdr3_condition = n()) %>%
      ungroup()%>%
      mutate(ratio_condition = cdr3_count_condition/total_cdr3_condition)%>% 
      mutate(ratio_condition_cluster =total_count_condition_cluster/total_cdr3_condition)%>% 
      group_by(cluster_id,Condition_2) %>%
      mutate(condition_2_mean_ratio = mean(ratio_condition_cluster))%>% 
      ungroup()
  
  # 将数据转换为宽格式
  myFunction <- function(x) {
    if(any(!is.na(x))) {
      return(x[which(!is.na(x))[1]])
    } else {
      return(x[1])
    }
  }
  seqs_wide <- dcast(processed_data, cluster_id ~ Condition_2, 
                     value.var = 'condition_2_mean_ratio', 
                     fun.aggregate = myFunction)
  # 移除NA并计算expand_ratio
  seqs_wide <- na.omit(seqs_wide) %>%
    mutate(expand_ratio = after_B0 / before_B1)
  
  # 更新原始数据中的expand_ratio和expand_class
  processed_data$expand_ratio <- seqs_wide$expand_ratio[match(processed_data$cluster_id,seqs_wide$cluster_id)]
  processed_data$expand_class <- "Undefined"
  processed_data$expand_class[processed_data$expand_ratio < 1] <- "Negative"
  processed_data$expand_class[processed_data$expand_ratio >= 1 & processed_data$expand_ratio < 1.25] <- "Low"
  processed_data$expand_class[processed_data$expand_ratio >= 1.25] <- "High"
  processed_data$expand_class[is.na(processed_data$expand_ratio)] <- "Undefined"
  
  # 因子化expand_class和Condition
  processed_data$expand_class <- factor(
    processed_data$expand_class,
    levels = c("Negative", "Low", "High", "Undefined")
  )
  processed_data[[condition_col]] <- factor(
    processed_data[[condition_col]],
    levels = c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")
  )
  
  return(processed_data)
}


# 定义处理数据的函数
process_and_save_sequences <- function(input_file, output_file, condition_col, seq_col, cluster_col, ct_col) {
  # 读取数据
  seqs <- read.csv(input_file)
  
  # 打印行数和列名
  cat("Number of rows:", nrow(seqs), "\n")
  cat("Column names:", colnames(seqs), "\n")
  
  # 处理序列数据
  processed_data <- process_sequences(seqs, condition_col, seq_col, cluster_col, ct_col)
  
  # 查看 expand_class 的结果
  print(table(processed_data$expand_class))
  
  # 保存处理后的数据
  write.csv(processed_data, output_file, row.names = FALSE)
}



# 使用示例
process_and_save_sequences("data/2_CD8T_all_edge_35.csv", "data/3_CD8T_seqs2_all_35_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")

process_and_save_sequences("data/2_CD8T_all_edge_25.csv", "data/3_CD8T_seqs2_all_25_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")

process_and_save_sequences("data/2_CD4T_all_edge_25.csv", "data/3_CD4T_seqs2_all_25_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")

process_and_save_sequences("data/2_CD4T_all_edge_35.csv", "data/3_CD4T_seqs2_all_35_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")
# 使用示例
process_and_save_sequences("data/2_CD8T_all_edge_15.csv", "data/3_CD8T_seqs2_all_15_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")

process_and_save_sequences("data/2_CD4T_all_edge_15.csv", "data/3_CD4T_seqs2_all_15_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")


# 使用示例
process_and_save_sequences("data/2_CD8T_all_edge_50.csv", "data/3_CD8T_seqs2_all_50_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")
process_and_save_sequences("data/2_CD4T_all_edge_50.csv", "data/3_CD4T_seqs2_all_50_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")

process_and_save_sequences("data/2_CD8T_all_edge_60.csv", "data/3_CD8T_seqs2_all_60_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")
process_and_save_sequences("data/2_CD4T_all_edge_60.csv", "data/3_CD4T_seqs2_all_60_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")

process_and_save_sequences("data/2_CD8T_all_edge_80.csv", "data/3_CD8T_seqs2_all_80_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")
process_and_save_sequences("data/2_CD4T_all_edge_80.csv", "data/3_CD4T_seqs2_all_80_with_expand_ratio.csv", 
                           "Condition", "cdr3_b_aa", "cluster_id", "ct_level_4")

library(ggplot2)
library(dplyr)

# 定义函数来计算和绘制图形
plot_data_visualization <- function(data, bar_output, boxplot_output) {
  # 设置自定义颜色
  bar_colors <- c("#B22C2CFF","#51A3CCFF", "#f18800","#85B22CFF","#E57E7EFF", "#BFB2FFFF")
  boxplot_colors <- c("#3cb346", "#c93f00", "#B22C2CFF")
  
  # 计算每个 expand_class 和 ct_level_4 的频数
  counts <- data %>%
    group_by(expand_class, ct_level_4) %>%
    summarise(count = n(), .groups = 'drop')
  
  # 绘制堆叠条形图
  p1 <- ggplot(counts, aes(x = ct_level_4, y = count, fill = expand_class)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Frequency of expand_class by ct_level_4",
         x = "Cell Type (ct_level_4)",
         y = "Count",
         fill = "Expand Class") +
    scale_fill_manual(values = bar_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 保存堆叠条形图
  ggsave(bar_output, p1, width = 6, height = 3)
  
  # 合并 expand_class 中的某些类别
  data$expand_class[data$expand_class %in% c("Low", "Undefined", "Negative")] <- "Others"
  
  # 绘制箱线图
  p2 <- ggplot(data, aes(x = Condition, y = ratio_condition_cluster)) +
    geom_boxplot(aes(fill = expand_class), outlier.shape = NA) +
    scale_fill_manual(values = boxplot_colors) +  # 使用自定义颜色
    geom_point(aes(color = expand_ratio), size = 0.3) +
    geom_line(aes(color = expand_ratio, group = cluster_id), linewidth = 0.3, linetype = "longdash") +
    labs(y = "Expand ratio", x = "Condition", fill = "expand_class") +
    scale_color_gradient2(low = "black", mid = "lightgray", high = "red", midpoint = 0.99) +
    ylim(0, 0.02) +
    theme_minimal() +
    theme(
      legend.position = "bottom", 
      legend.direction = "horizontal",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.box = "vertical",
      legend.text = element_text(size = 4),
      plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    facet_wrap(~expand_class, scales = "free_y", ncol = 2)
  
  # 保存箱线图
  ggsave(boxplot_output, p2, height = 4, width = 8)
}


library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(ggseqlogo)

# 定义函数来处理数据并生成图形
generate_plots <- function(data, bar_output, boxplot_output, logo_output_prefix) {
    # 1. 根据条件筛选出符合要求的 top_clusters
    top_clusters <- data %>%
        group_by(cluster_id) %>%
        filter(n_distinct(Condition) > 8) %>%  # 筛选出 cluster_id 对应的 unique(data$Condition) 个数大于 5 的
        ungroup()%>% 
        arrange(desc(expand_ratio)) %>%
        distinct(cluster_id, .keep_all = TRUE) %>%
        slice_head(n = 10) %>%
        pull(cluster_id)
    
    # 2. 根据选出的 cluster_id 过滤数据
    seqs_choose <- data %>%
        filter(cluster_id %in% top_clusters)
  # # 1. 选出expand_ratio降序排列的不同cluster_id的前10个
  # top_clusters <- data %>%
  #   arrange(desc(expand_ratio)) %>%
  #   distinct(cluster_id, .keep_all = TRUE) %>%
  #   slice_head(n = 10) %>%
  #   pull(cluster_id)

  # # 2. 根据选出的cluster_id过滤数据
  # seqs_choose <- data %>%
  #   filter(cluster_id %in% top_clusters)
    
  # 打印过滤后的数据集行数
  print(nrow(seqs_choose))
  
  # 3. 生成序列logo图
  pdf(logo_output_prefix, height = 4.6, width = 3.2)
      plot_list <- list()
      position_freq_list <- list()
      unique_cluster_id <- unique(seqs_choose$cluster_id)
      
      for (id in unique_cluster_id) {
        seqs_aa <- seqs_choose %>% filter(cluster_id == id)
        max_char <- max(nchar(seqs_aa$cdr3_b_aa))
        
        # 初始化位置频率矩阵
        position_freq <- matrix(0, nrow = 26, ncol = max_char, dimnames = list(LETTERS, 1:max_char))
        
        # 统计每个位置的氨基酸频率
        for (seq in seqs_aa$cdr3_b_aa) {
          for (i in 1:nchar(seq)) {
            aa <- substr(seq, i, i)
            if (aa %in% LETTERS) {
              position_freq[match(aa, LETTERS), i] <- position_freq[match(aa, LETTERS), i] + 1
            }
          }
        }
        
        position_freq_list[[id]] <- position_freq
        
        # 生成序列logo图
        p1 <- ggseqlogo(position_freq, method = "bits") + ggtitle(paste("Cluster", id, "Bits"))
        p2 <- ggseqlogo(position_freq, method = "prob") + ggtitle(paste("Cluster", id, "Probability"))
        plot_save <- grid.arrange(p1, p2)
        plot_list[[id]] <- plot_save
        
      }
  dev.off()
  
  # 4. 绘制所有cluster的expand ratio图
  pdf(bar_output, height = 6, width = 12)
  p3 <- ggplot(seqs_choose, aes(x = Condition, y = ratio_condition_cluster,
                                group = cluster_id, color = expand_ratio)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.6, linetype = "longdash") +
    scale_color_gradient2(low = "black", mid = "lightgray", high = "red", midpoint = 0.99) +
    labs(color = "Expand_ratio") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.text = element_text(size = 12),
      legend.key.width = unit(2, "cm")
    )
  print(p3)
  dev.off()
  
  # 5. 绘制带有facet的expand ratio图
  pdf(boxplot_output, height = 6.5, width = 3)
  p4 <- ggplot(seqs_choose, aes(x = Condition, y = ratio_condition_cluster,
                                group = cluster_id, color = expand_ratio)) +
    geom_point(size = 1) +
    geom_line( linetype = "longdash") +#linewidth = 1,
    theme_classic() +
    scale_color_gradient2(low = "black", mid = "lightgray", high = "red", midpoint = 0.8) +
    labs(color = "expand_ratio") +
    theme(
      strip.background = element_blank(),
      strip.placement = 'outside',
      strip.switch.pad.grid = unit(1, "cm"),
      panel.spacing = unit(1, "lines"),
      panel.border = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5,color="black"),
        axis.text.y = element_text(size = 6,color="black"),
      legend.text = element_text(size = 6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.width = unit(1, "cm")
    ) +
    facet_wrap(~cluster_id, scales = "free", ncol = 2)
  print(p4)
  dev.off()
}


# 定义主函数
generate_visualizations <- function(file_path,dir_create) {
  # 读取数据
  processed_data <- read.csv(file_path)
  
  # 从文件路径提取信息并创建目录
  dir_create <- dir_create
  if (!dir.exists(dir_create)) {
    dir.create(dir_create, recursive = TRUE)
  }
  
  # 调用函数进行数据可视化
    # 画图展示 expand_class
    plot_data_visualization(processed_data, 
                            file.path(dir_create, "1_expand_class_count.pdf"), 
                            file.path(dir_create, "2_boxplot_cluster.pdf"))
    
    # 画图展示 top_10_seq
    generate_plots(processed_data, 
                   file.path(dir_create, "3_top_10_cluster.pdf"),
                   file.path(dir_create, "3_top_10_cluster_Condition.pdf"),  # 修正文件后缀
                   file.path(dir_create, "4_top_10_seq_Bits.pdf"))
}



generate_visualizations("data/3_CD4T_seqs2_all_35_with_expand_ratio.csv","figures/CD4T_35")
generate_visualizations("data/3_CD8T_seqs2_all_35_with_expand_ratio.csv","figures/CD8T_35")

generate_visualizations("data/3_CD4T_seqs2_all_35_with_expand_ratio.csv","figures/CD4T_35")
generate_visualizations("data/3_CD4T_seqs2_all_25_with_expand_ratio.csv","figures/CD4T_25")
generate_visualizations("data/3_CD4T_seqs2_all_15_with_expand_ratio.csv","figures/CD4T_15")

generate_visualizations("data/3_CD8T_seqs2_all_35_with_expand_ratio.csv","figures/CD8T_35")
generate_visualizations("data/3_CD8T_seqs2_all_25_with_expand_ratio.csv","figures/CD8T_25")
generate_visualizations("data/3_CD8T_seqs2_all_15_with_expand_ratio.csv","figures/CD8T_15")

generate_visualizations("data/3_CD8T_seqs2_all_50_with_expand_ratio.csv","figures/CD8T_50")
generate_visualizations("data/3_CD8T_seqs2_all_60_with_expand_ratio.csv","figures/CD8T_60")
generate_visualizations("data/3_CD8T_seqs2_all_80_with_expand_ratio.csv","figures/CD8T_80")

generate_visualizations("data/3_CD4T_seqs2_all_50_with_expand_ratio.csv","figures/CD4T_50")
generate_visualizations("data/3_CD4T_seqs2_all_60_with_expand_ratio.csv","figures/CD4T_60")
generate_visualizations("data/3_CD4T_seqs2_all_80_with_expand_ratio.csv","figures/CD4T_80")


