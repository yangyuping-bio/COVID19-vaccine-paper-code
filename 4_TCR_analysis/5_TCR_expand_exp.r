library("GSEABase")
library("tidyverse")
library(clusterProfiler)
library(ClusterGVis)
library(org.Hs.eg.db)
library(tidyverse)
library(Seurat)

CD8T <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD8T_TCR.rds")
CD4T <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD4T_TCR.rds")

library(Seurat)
library(EnhancedVolcano)
library(future)

# 定义函数
analyze_tcell_expansion <- function(input_file, tcr_data, output_diff_file, output_volcano_file,xlim=NULL,ylim=NULL) {
  # 读取数据
  CD8T_seqs <- read.csv(input_file)
  
  # 检查并重命名 expand_class
  CD8T_seqs$expand_class <- as.character(CD8T_seqs$expand_class)
  CD8T_seqs$expand_class[CD8T_seqs$expand_class %in% c("Low", "Undefined", "Negative")] <- "Others" 
  CD8T_seqs$expand_class <- factor(
    CD8T_seqs$expand_class,
    levels = c("High", "Others")
  )
  
  # 查看 expand_class 的分布
  print(table(CD8T_seqs$expand_class))
  
  # 匹配 expand_class 到 TCR 数据
  tcr_data$expand_class <- CD8T_seqs$expand_class[match(tcr_data$tcr_barcode, CD8T_seqs$barcode)]
  tcr_data$expand_class[is.na(tcr_data$expand_class)] <- "Others"
  
  # 查看更新后的 expand_class 分布
  print(table(tcr_data$expand_class))
  
  # 设定默认 Assay 并计算差异标记
  DefaultAssay(tcr_data) <- "RNA"
  Idents(tcr_data) <- "expand_class"
  
  # 设置并行计算
  plan("multicore", workers = 10)
  diff_markers <- FindMarkers(tcr_data, 
                               min.pct = 0.1, 
                               logfc.threshold = 0.05,
                               group.by = "expand_class",
                               ident.1 = "High",
                               ident.2 = "Others")    
  
  # 保存差异标记数据
  write.csv(diff_markers, output_diff_file)
  print(head(diff_markers))
  
  # 绘制火山图
  pdf(output_volcano_file, height = 6, width = 6.2)
  print(EnhancedVolcano(diff_markers,
                  lab = rownames(diff_markers),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  pCutoff = 0.05,
                  FCcutoff = 0.08,
                  pointSize = 2.0,
                  labSize = 6,
                  xlim = xlim,
                  ylim = ylim,
                  legendLabSize = 12,
                  legendIconSize = 6,
                  title = 'DEGs in High TCR Expansion')
        )
  dev.off()
}



unique(CD4T$ct_level_4)

Idents(CD4T) <- 'ct_level_4'
nrow(CD4T@meta.data)
CD4_Treg <- subset(CD4T,idents=c('CD4_Treg'))
nrow(CD4_Treg@meta.data)

# 使用示例
analyze_tcell_expansion("data/3_CD4T_seqs2_all_35_with_expand_ratio.csv", 
                         CD4_Treg, 
                         "data/6_diff_CD4_Treg_Markers_expand_class.csv", 
                         "figures/6_diff_CD4_Treg_expand_class.pdf",
                       xlim =c(-0.5, 1.5),ylim =c(0, 5))



# 使用示例
analyze_tcell_expansion("data/3_CD8T_seqs2_all_35_with_expand_ratio.csv", 
                         CD8T, 
                         "data/4_diff_CD8T_Markers_expand_class.csv", 
                         "figures/4_diff_CD8T_expand_class.pdf",
                       xlim =c(-1.5, 2),ylim =c(0, 145))

# 使用示例
analyze_tcell_expansion("data/3_CD4T_seqs2_all_35_with_expand_ratio.csv", 
                         CD4T, 
                         "data/4_diff_CD4T_Markers_expand_class.csv", 
                         "figures/4_diff_CD4T_expand_class.pdf",
                       xlim =c(-0.5, 1.5),ylim =c(0, 42))

source("/share/home/qlab/projects/project_cvv/yyp_results_3/yyp_function.R") #加载到当前的R环境中

diff_markers<- read.csv("data/4_diff_CD4T_Markers_expand_class.csv")
rownames(diff_markers) <- diff_markers$X
head(diff_markers)

# 绘制火山图
diff_markers_f <- diff_markers%>% filter(-log(p_val_adj) <150)
plot_volcano(diff_markers_f, ident_1 = "Expand_High", ident_2 = "Expand_Low",
             width = 5, height = 2.5,
            output_file = paste0("figures/2_CD4T_35_volcano.pdf"),
            FC_cutoff = 0.15, p_cutoff = 0.05,max.overlaps = 150,
             logFC_thresh = 0.25, pct_thresh = 0.05,
             pt_size = 0.5, geom_text_size = 4, 
             line_size = 0.5, geom_line_size = 0.05)


# GO富集分析
GO_result <- perform_GO_analysis(diff_markers, CD4T, top_n = 50,
                                 output_file="data/6_CD4T_35_GO_results_Top_50.csv",
                                 ident_1="Expand_High", ident_2="Expand_Low",
                                 FC_cutoff = 0.10, p_cutoff = 0.05) 
head(GO_result$Description,n=20)

GO_result <- read.csv("data/6_CD4T_35_GO_results_Top_50.csv")

# 绘制GO富集结果
 plot_GO_results(GO_result, top_n = 30, 
                 ident_1 = "Expand_High", ident_2 = "Expand_Low", Only_ident1 =F,
                 output_file="figures/2_CD4T_Expand_High_vs_Low_GO_results_ALL_30.pdf",
                 color_map = c("Expand_High" ="#ffd401", ident_2 = "#0f7b9f"),
                 width = 4, height = 6,ylim_range=c(-10, 15)
                )

# 绘制GO富集结果
 plot_GO_results(GO_result, top_n = 10, 
                 ident_1 = "Expand_High", ident_2 = NULL, Only_ident1 =T,
                 output_file="figures/2_CD4T_Expand_High_vs_Low_GO_results_10.pdf",
                 color_map = c("Expand_High" = "#ffd401", "Expand_Low" = "#0f7b9f"),
                 width = 4, height = 2,ylim_range=c(-10, 15)
                )





diff_markers<- read.csv("data/4_diff_CD8T_Markers_expand_class.csv")
rownames(diff_markers) <- diff_markers$X
head(diff_markers)

# 绘制火山图
#diff_markers_f <- diff_markers%>% filter(abs(avg_log2FC) <3)
plot_volcano(diff_markers, ident_1 = "Expand_High", ident_2 = "Expand_Low",
             width = 5, height = 2.5,
            output_file = paste0("figures/2_CD8T_35_volcano.pdf"),
            FC_cutoff = 0.5, p_cutoff = 0.05,max.overlaps = 30,
             logFC_thresh = 0.68, pct_thresh = 0.15,
             pt_size = 0.5, geom_text_size = 4, 
             line_size = 0.5, geom_line_size = 0.05)


# GO富集分析
GO_result <- perform_GO_analysis(diff_markers, CD8T, top_n = 50,
                                 output_file="data/6_CD8T_35_GO_results_Top_50.csv",
                                 ident_1="Expand_High", ident_2="Expand_Low",
                                 FC_cutoff = 0.5, p_cutoff = 0.05) 
head(GO_result$Description,n=20)

GO_result <- read.csv("data/6_CD8T_35_GO_results_Top_50.csv")

# 绘制GO富集结果
plot_GO_results(GO_result, top_n = 60, 
                 ident_1 = "Expand_High", ident_2 = "Expand_Low", Only_ident1 =F,
                 output_file="figures/2_CD8T_Expand_High_vs_Low_GO_results_ALL_30.pdf",
                 color_map = c("Expand_High" = "#ffd401", ident_2 = "#0f7b9f"),
                 width = 4, height = 6,ylim_range=c(-10, 10)
                )

# 绘制GO富集结果
 plot_GO_results(GO_result, top_n = 10, 
                 ident_1 = "Expand_High", ident_2 = NULL, Only_ident1 =T,
                 output_file="figures/2_CD8T_Expand_High_vs_Low_GO_results_10.pdf",
                 color_map = c("Expand_High" = "#ffd401", "Expand_Low" = "#0f7b9f"),
                 width = 4, height = 2,ylim_range=c(-10, 10)
                )



library(org.Hs.eg.db) ##加载人类
 require(clusterProfiler)

# perform GO
GO_a <- function(interest_gene, global_gene) {
  # Convert gene symbols to ENTREZ IDs for interest and global genes
  gene_entrez_id2GO <-
    clusterProfiler::bitr(interest_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)

  universe_gene_entrez_id2GO <-
    clusterProfiler::bitr(global_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)

  # Perform GO enrichment analysis and simplify results
  enr_res <- clusterProfiler::enrichGO(
    gene = gene_entrez_id2GO$ENTREZID,
    universe = universe_gene_entrez_id2GO$ENTREZID,
    ont = "BP",
    OrgDb = 'org.Hs.eg.db'
  )
  
  enr_res2 <- clusterProfiler::simplify(enr_res)
  
  # Return results and plots as a list, including the actual genes used
  list(
    #go_res = enr_res,
    go_res_simplify = enr_res2,
    actual_genes = gene_entrez_id2GO$SYMBOL # Store actual genes
   # p1 = goplot(enr_res),
    #p2 = barplot(enr_res2, showCategory = 20, color = "p.adjust")
  )
}

# Process and combine GO results
process_GO <- function(GO_data, celltype) {
  go_results <- GO_data$go_res_simplify@result %>%
    head(50) %>%
    mutate(
      celltype = celltype,
      log10q = -log10(abs(qvalue)),
      actual_genes = paste(GO_data$actual_genes, collapse = ", ")  # Add actual genes
    )
  
  return(go_results)
}


all_markers = read.csv("data/4_diff_CD4T_Markers_expand_class.csv")
all_markers$gene <- all_markers$X
head(all_markers)
Expand_low = all_markers %>% filter(p_val_adj < 0.05, avg_log2FC < 0)
Expand_high = all_markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5)

Expand_low_go = GO_a(Expand_low$gene,rownames(CD4T))
Expand_high_go = GO_a(Expand_high$gene,rownames(CD4T))

GO_result <- bind_rows(
  process_GO(Expand_high_go, "CD4T_high_Expand"),
  process_GO(Expand_low_go, "CD4T_low_Expand")
)
GO_result$log10q = ifelse(GO_result$celltype == "CD4T_high_Expand", 
                          GO_result$log10q, GO_result$log10q* -1)
write.csv(GO_result,"data/5_CD4T_35_GO_results_Top_50.csv")

ggplot(GO_result) +
  geom_bar(aes(x = reorder(Description, log10q),y = log10q,fill = celltype),stat = "identity",color = "white") +
  geom_text(aes(x = reorder(Description, log10q),y = 0,label = Description),size = 7 * 0.35,angle = 0) +
  scale_fill_manual(name = "",values = c("CD4T_high_Expand" = "#d83215", "CD4T_low_Expand" = "#0f7b9f") ) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    text = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "") #+
  #ylim(-60, 60)
ggsave("figures/5_CD4T_Expand_GO_results_Top_50.pdf",width = 5,height = 10) 


GO_result_f <- GO_result%>% filter(celltype == 'CD4T_high_Expand')%>% head(n=20)
ggplot(GO_result_f) +
  geom_bar(aes(x = reorder(Description, log10q),y = log10q,fill = celltype),stat = "identity",color = "white") +
  geom_text(aes(x = reorder(Description, log10q),y = 0,label = Description),size = 7 * 0.35,angle = 0) +
  scale_fill_manual(name = "",values = c("CD4T_high_Expand" = "#BFB2FFFF", "CD8T_low_Expand" = "#0f7b9f") ) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    text = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "")+
  ylim(-5, 5)
ggsave("figures/5_CD4T_Expand_GO_results_high.pdf",width = 5,height = 3.2) 

all_markers = read.csv("data/4_diff_CD8T_Markers_expand_class.csv")
all_markers$gene <- all_markers$X
head(all_markers)
Expand_low = all_markers %>% filter(p_val_adj < 0.05, avg_log2FC < 0)
Expand_high = all_markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5)

Expand_low_go = GO_a(Expand_low$gene,rownames(CD8T))
Expand_high_go = GO_a(Expand_high$gene,rownames(CD8T))

GO_result <- bind_rows(
  process_GO(Expand_high_go, "CD8T_high_Expand"),
  process_GO(Expand_low_go, "CD8T_low_Expand")
)
GO_result$log10q = ifelse(GO_result$celltype == "CD8T_high_Expand", 
                          GO_result$log10q, GO_result$log10q* -1)
write.csv(GO_result,"data/5_CD8T_35_GO_results_Top_50.csv")


ggplot(GO_result) +
  geom_bar(aes(x = reorder(Description, log10q),y = log10q,fill = celltype),stat = "identity",color = "white") +
  geom_text(aes(x = reorder(Description, log10q),y = 0,label = Description),size = 7 * 0.35,angle = 0) +
  scale_fill_manual(name = "",values = c("CD8T_high_Expand" = "#d83215", "CD8T_low_Expand" = "#0f7b9f") ) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    text = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "") 
  #ylim(-60, 60)
ggsave("figures/5_CD8T_Expand_GO_results_Top_50.pdf",width = 5,height = 10) 

GO_result_f <- GO_result%>% filter(celltype == 'CD8T_high_Expand')%>% head(n=20)
ggplot(GO_result_f) +
  geom_bar(aes(x = reorder(Description, log10q),y = log10q,fill = celltype),stat = "identity",color = "white") +
  geom_text(aes(x = reorder(Description, log10q),y = 0,label = Description),size = 7 * 0.35,angle = 0) +
  scale_fill_manual(name = "",values = c("CD8T_high_Expand" = "#fc9272", "CD8T_low_Expand" = "#0f7b9f") ) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    text = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "")+
  ylim(-10, 10)
ggsave("figures/5_CD8T_Expand_GO_results_high.pdf",width = 5,height = 3.2) 




