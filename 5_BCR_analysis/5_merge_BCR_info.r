library("tidyverse")
library("data.table")
library(ggplot2)
library("Matrix")
library(Seurat)
library(viridis)
library(RColorBrewer)

Bcell <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_BCR/1_BCR_info/data/1_Bcell_BCR.rds")

viral_info <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/6_BCR/2_BCR_viral/data/Bcell_data_res_distance_3.csv")

SHM_info <- read.csv("/share/home/qlab/projects/project_cvv/yyp_results_3/6_BCR/3_BCR_SHM/Bcell_BCR_SHM.csv")

Bcell@meta.data<- Bcell@meta.data[ ,colnames(Bcell@meta.data)[c(1:25,68:107)]]
colnames(Bcell@meta.data)

Bcell$cdr3_IGH_aa <- viral_info$cdr3_IGH_aa[match(Bcell$bcr_barcode,viral_info$bcr_barcode)]
Bcell$Virus_name <- viral_info$Virus_name[match(Bcell$bcr_barcode,viral_info$bcr_barcode)]
Bcell$Virus_distance <- viral_info$Chain_distance[match(Bcell$bcr_barcode,viral_info$bcr_barcode)]
table(Bcell$Virus_name)

Bcell$SHM_rate <- SHM_info$SHM_rate[match(Bcell$bcr_barcode,SHM_info$bcr_barcode)]
Bcell$shm_rate_ratio <- SHM_info$shm_rate_ratio[match(Bcell$bcr_barcode,SHM_info$bcr_barcode)]
summary(Bcell$SHM_rate)

saveRDS(Bcell, "data/1_Bcell_BCR_with_info.rds")

Bcell$SHM_levels <- ifelse(Bcell$SHM_rate < 4.688, "Low", "Median")
Bcell$SHM_levels <- ifelse(Bcell$SHM_rate > 5.523, "High", Bcell$SHM_levels)
Bcell$SHM_levels <- factor(Bcell$SHM_levels, levels = c("Low", "Median", "High"))
table(Bcell$SHM_levels)

sum(is.na(Bcell$bcr_cdr3_nt))

table(Bcell$SHM_levels,Bcell$ct_level_4)

plot_df = Bcell@meta.data %>%
  filter(! is.na("bcr_cdr3_nt"))%>% 
  group_by(ct_level_4) %>%
  mutate(median_SHM = median(SHM_rate, na.rm = TRUE)) %>% 
  mutate(group = ifelse(SHM_rate > median_SHM, "High","Low")) %>%
  as.data.frame()

head(plot_df)

sub_obj = subset(Bcell,bcr_barcode %in% plot_df$bcr_barcode)
sub_obj$SHM_group = plot_df$group
Idents(sub_obj) = sub_obj$SHM_group
table(sub_obj$SHM_group,sub_obj$ct_level_4)

# Perform differential expression analysis
bcr_COVID_diff <- FindMarkers(sub_obj, ident.1 = "High", ident.2 = "Low", 
                               logfc.threshold = 0.2, min.pct = 0.01)

# Set thresholds
FC_cutoff <- 0.3
p_cutoff <- 0.05

# Calculate percentage difference and classify differential expression
bcr_COVID_diff <- bcr_COVID_diff %>%
  dplyr::mutate(
    pct_value = abs(pct.1 - pct.2),
    DE = case_when(
      avg_log2FC > FC_cutoff & p_val < p_cutoff ~ "SHM-high group",
      avg_log2FC < -FC_cutoff & p_val < p_cutoff ~ "SHM-low group",
      TRUE ~ "Not significant"
    ),
    label = ifelse(
      abs(avg_log2FC) > 0.3 & p_val < 0.001 & pct_value > 0.1,
      rownames(bcr_COVID_diff),
      NA
    )
  )
write.csv("data/2.1_SHM_Diff_all.csv")

# Define colors for the plot
mycolors <- c("#0F7B9F", "#D83215", "grey")
names(mycolors) <- c("SHM-low group", "SHM-high group", "Not significant")

# Create the plot
ggplot(bcr_COVID_diff, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DE, label = label)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = mycolors, name = "Enriched in") +
  geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(size = 1.75, show.legend = FALSE, max.overlaps = 40) +
  cowplot::theme_cowplot() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.position = "right",
    plot.margin = margin(1, 1, 1, 1, "char"),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3)
  ) +
  labs(x = "Log2(fold change)", y = "-Log10(q-value)")

# Save the plot
ggsave("figures/2.1_SHM_Diff_ALL.pdf", width = 4, height = 2)

# perform DE

DE_sum = data.frame()

for (i in sort(unique(sub_obj$ct_level_4))) {
  tmp_obj = subset(sub_obj, ct_level_4 == i)
  DE = Seurat::FindMarkers(
    tmp_obj,
    ident.1 = "High",
    ident.2 = "Low",
    logfc.threshold = 0
  )
  DE$gene = rownames(DE)
  DE$ct_level_4 = i
  DE_sum = rbind(DE_sum, DE)
}

FC_cutoff <- 0.15
p_cutoff <- 0.05

plot_df <- DE_sum 
plot_df$pct_value = abs(plot_df$pct.1-plot_df$pct.2)

plot_df$DE <- "Not significant"
plot_df$DE[plot_df$avg_log2FC > FC_cutoff & plot_df$p_val < p_cutoff ] <- "SHM-high group"
plot_df$DE[plot_df$avg_log2FC < -FC_cutoff & plot_df$p_val < p_cutoff] <- "SHM-low group"

plot_df$label <-
  ifelse(
    abs(plot_df$avg_log2FC) > 0.15 &
      plot_df$p_val < 0.001 & plot_df$pct_value > 0.01,
    plot_df$gene,
    NA
  )
table(plot_df$DE)

mycolors <- c("#0F7B9F", "#D83215", "grey")
names(mycolors) <- c("SHM-low group", "SHM-high group", "Not significant")
plot_df$ct_level_4 = factor(plot_df$ct_level_4,levels = c('Bm_CD27+_IGHM-','Bm_CD27+_IGHM+',
                                                          'Bm_CD27+_IGHM+_SOX5','Plasma'))
write.csv("data/2.1_SHM_level_Diff_all.csv")
ggplot(data=plot_df, aes(x=avg_log2FC, y=-log10(p_val_adj), color = DE, label = label)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = mycolors, name = paste0("Enriched in")) +
  geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype="dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color = "black") +
  ggrepel::geom_text_repel(size = 5 * 0.35, show.legend = FALSE, max.overlaps = 40) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)) +
  xlab("Log2(fold change)") +
  ylab("-Log10(q-value)") + 
  facet_wrap(~ct_level_4, scales = "free")

ggsave("figures/2.1_SHM_Class_Diff_Celltype.pdf",
       width = 10,
       height = 4)


table(Bcell$Virus_name)

##--- Figure S2C; Differentially expressed genes between the SHM-high (red) and SHM-low (blue) groups
plot_df = Bcell@meta.data %>%
  filter(! is.na("bcr_cdr3_nt"))%>% 
  group_by(ct_level_4) %>% 
  mutate(COVID_BCR = ifelse(Virus_name == 'COVID', "YES","NO")) %>%
  as.data.frame()

head(plot_df)

sub_obj = subset(Bcell,bcr_barcode %in% plot_df$bcr_barcode)
sub_obj$COVID_BCR = plot_df$COVID_BCR
Idents(sub_obj) = sub_obj$COVID_BCR
table(sub_obj$COVID_BCR,sub_obj$ct_level_4)

# Perform differential expression analysis
bcr_COVID_diff <- FindMarkers(sub_obj, ident.1 = "YES", ident.2 = "NO", 
                               logfc.threshold = 0.2, min.pct = 0.01)

# Set thresholds
FC_cutoff <- 0.3
p_cutoff <- 0.05

# Calculate percentage difference and classify differential expression
bcr_COVID_diff <- bcr_COVID_diff %>%
  dplyr::mutate(
    pct_value = abs(pct.1 - pct.2),
    DE = case_when(
      avg_log2FC > FC_cutoff & p_val < p_cutoff ~ "COVID BCR group",
      avg_log2FC < -FC_cutoff & p_val < p_cutoff ~ "Not COVID group",
      TRUE ~ "Not significant"
    ),
    label = ifelse(
      abs(avg_log2FC) > 0.3 & p_val < 0.001 & pct_value > 0.1,
      rownames(bcr_COVID_diff),
      NA
    )
  )
write.csv("data/2.2_COVID_BCR_Diff_all.csv")

# Define colors for the plot
mycolors <- c("#0F7B9F", "#D83215", "grey")
names(mycolors) <- c("Not COVID group", "COVID BCR group", "Not significant")

# Create the plot
ggplot(bcr_COVID_diff, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DE, label = label)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = mycolors, name = "Enriched in") +
  geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(size = 1.75, show.legend = FALSE, max.overlaps = 40) +
  cowplot::theme_cowplot() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.position = "right",
    plot.margin = margin(1, 1, 1, 1, "char"),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3)
  ) +
  labs(x = "Log2(fold change)", y = "-Log10(q-value)")

# Save the plot
ggsave("figures/2.2_COVID_BCR_Diff_ALL.pdf", width = 4, height = 2)



# perform DE
DE_sum = data.frame()
for (i in sort(unique(sub_obj$ct_level_4))) {
  tmp_obj = subset(sub_obj, ct_level_4 == i)
  DE = Seurat::FindMarkers(
    tmp_obj,
    ident.1 = "YES",
    ident.2 = "NO",
    logfc.threshold = 0
  )
  DE$gene = rownames(DE)
  DE$ct_level_4 = i
  DE_sum = rbind(DE_sum, DE)
}

FC_cutoff <- 0.3
p_cutoff <- 0.05

plot_df <- DE_sum
plot_df$pct_value = abs(plot_df$pct.1-plot_df$pct.2)

plot_df$DE <- "Not significant"
plot_df$DE[plot_df$avg_log2FC > FC_cutoff & plot_df$p_val < p_cutoff ] <- "COVID BCR group"
plot_df$DE[plot_df$avg_log2FC < -FC_cutoff & plot_df$p_val < p_cutoff] <- "Not COVID BCR group"

plot_df$label <-
  ifelse(
    abs(plot_df$avg_log2FC) > 0.3 &
      plot_df$p_val < 0.001 & plot_df$pct_value > 0.15,
    plot_df$gene,
    NA
  )
table(plot_df$DE)
write.csv("data/2.2_COVID_BCR_Diff_celltype.csv")

# Define colors for the plot
mycolors <- c("#0F7B9F", "#D83215", "grey")
names(mycolors) <- c("Not COVID BCR group", "COVID BCR group", "Not significant")
plot_df$ct_level_4 = factor(plot_df$ct_level_4,levels = c("Bmn_IGHD",'Bm_CD27+_IGHM-',
                                                          'Bm_CD27+_IGHM+',
                                                          'Bm_CD27+_IGHM+_SOX5','Plasma'))

ggplot(data=plot_df, aes(x=avg_log2FC, y=-log10(p_val_adj), color = DE, label = label)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = mycolors, name = paste0("Enriched in")) +
  geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype="dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color = "black") +
  ggrepel::geom_text_repel(size = 5 * 0.35, show.legend = FALSE, max.overlaps = 40) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)) +
  xlab("Log2(fold change)") +
  ylab("-Log10(q-value)") + 
  facet_wrap(~ct_level_4, scales = "free")

ggsave("figures/2.2_COVID_BCR_Diff_celltype.pdf",
       width = 10,
       height = 4)




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
  cell_type_obj$Clonotype_num_ld[cell_type_obj$Clonotype_num > 1] <- "n > 1"
  
  # 将 Clonotype_num_ld 转换为因子，并设定顺序
  cell_type_obj$Clonotype_num_ld <- factor(
    cell_type_obj$Clonotype_num_ld,
    levels = c("Not detected", "n = 1", "n > 1")
  )
  
  return(cell_type_obj)
}

ct_seq <- c('Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM+','Plasma','Bmn_IGHD','Bm_CD27+_IGHM-')
clonotype_num <- as.data.frame(table(Bcell$bcr_clonotype_id_new))
Bcell$Clonotype_num <- clonotype_num$Freq[match(Bcell$bcr_clonotype_id_new, clonotype_num$Var1)]
table(Bcell$Clonotype_num)
Bcell <- process_clonotype_data(Bcell,ct_seq)
table(Bcell$ct_level_4,Bcell$Clonotype_num_ld)

# Perform differential expression analysis
bcr_expand_diff <- FindMarkers(sub_obj, ident.1 = "YES", ident.2 = "NO", 
                               logfc.threshold = 0.2, min.pct = 0.01)

# Set thresholds
FC_cutoff <- 0.3
p_cutoff <- 0.05

# Calculate percentage difference and classify differential expression
bcr_expand_diff <- bcr_expand_diff %>%
  dplyr::mutate(
    pct_value = abs(pct.1 - pct.2),
    DE = case_when(
      avg_log2FC > FC_cutoff & p_val < p_cutoff ~ "Expand group",
      avg_log2FC < -FC_cutoff & p_val < p_cutoff ~ "Not Expand group",
      TRUE ~ "Not significant"
    ),
    label = ifelse(
      abs(avg_log2FC) > 0.3 & p_val < 0.001 & pct_value > 0.12,
      rownames(bcr_expand_diff),
      NA
    )
  )

# Define colors for the plot
mycolors <- c("#0F7B9F", "#D83215", "grey")
names(mycolors) <- c("Not Expand group", "Expand group", "Not significant")

# Create the plot
ggplot(bcr_expand_diff, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DE, label = label)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = mycolors, name = "Enriched in") +
  geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(size = 1.75, show.legend = FALSE, max.overlaps = 40) +
  cowplot::theme_cowplot() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.position = "right",
    plot.margin = margin(1, 1, 1, 1, "char"),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3)
  ) +
  labs(x = "Log2(fold change)", y = "-Log10(q-value)")

# Save the plot
ggsave("figures/2.3_Expand_Class_Diff_ALL.pdf", width = 4, height = 2)

plot_df = Bcell@meta.data %>%
  filter(!Clonotype_num_ld == "Not detected")%>% 
  group_by(ct_level_4) %>%
  mutate(bcr_expand = ifelse(Clonotype_num_ld == "n = 1", "YES","NO")) %>%
  as.data.frame()

head(plot_df)

sub_obj = subset(Bcell,bcr_barcode %in% plot_df$bcr_barcode)
sub_obj$bcr_expand = plot_df$bcr_expand
Idents(sub_obj) = sub_obj$bcr_expand
table(sub_obj$bcr_expand,sub_obj$ct_level_4)

# perform DE
DE_sum = data.frame()
for (i in sort(unique(sub_obj$ct_level_4))) {
  tmp_obj = subset(sub_obj, ct_level_4 == i)
  DE = Seurat::FindMarkers(
    tmp_obj,
    ident.1 = "YES",
    ident.2 = "NO",
    logfc.threshold = 0
  )
  DE$gene = rownames(DE)
  DE$ct_level_4 = i
  DE_sum = rbind(DE_sum, DE)
}

FC_cutoff <- 0.3
p_cutoff <- 0.05

plot_df <- DE_sum
plot_df$pct_value = abs(plot_df$pct.1-plot_df$pct.2)

plot_df$DE <- "Not significant"
plot_df$DE[plot_df$avg_log2FC > FC_cutoff & plot_df$p_val < p_cutoff ] <- "Expand group"
plot_df$DE[plot_df$avg_log2FC < -FC_cutoff & plot_df$p_val < p_cutoff] <- "Not Expand group"

plot_df$label <-
  ifelse(
    abs(plot_df$avg_log2FC) > 0.3 &
      plot_df$p_val < 0.001 & plot_df$pct_value > 0.15,
    plot_df$gene,
    NA
  )
table(plot_df$DE)

write.csv("data/2.3_Expand_Class_Diff_celltype.csv")
# Define colors for the plot
mycolors <- c("#0F7B9F", "#D83215", "grey")
names(mycolors) <- c("Not Expand group", "Expand group", "Not significant")
plot_df$ct_level_4 = factor(plot_df$ct_level_4,levels = c("Bmn_IGHD",'Bm_CD27+_IGHM-',
                                                          'Bm_CD27+_IGHM+',
                                                          'Bm_CD27+_IGHM+_SOX5','Plasma'))

ggplot(data=plot_df, aes(x=avg_log2FC, y=-log10(p_val_adj), color = DE, label = label)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = mycolors, name = paste0("Enriched in")) +
  geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype="dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color = "black") +
  ggrepel::geom_text_repel(size = 5 * 0.35, show.legend = FALSE, max.overlaps = 40) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)) +
  xlab("Log2(fold change)") +
  ylab("-Log10(q-value)") + 
  facet_wrap(~ct_level_4, scales = "free")

ggsave("figures/2.3_Expand_Class_Diff_celltype.pdf",
       width = 10,
       height = 5)


Bcell$ct_level_4 = factor(Bcell$ct_level_4,levels = c("Bmn_IGHD",
                                                          'Bm_CD27+_IGHM+',
                                                          'Bm_CD27+_IGHM+_SOX5',
                                                      'Bm_CD27+_IGHM-','Plasma'))
obj <- Bcell
Idents(obj) = obj$ct_level_4
features <- grep(pattern = "^IG", x = rownames(obj), value = TRUE) 
percent.featureset <- colSums(x = GetAssayData(object = obj, slot = "counts")[features, , drop = FALSE])/ 
  colSums(x = obj@assays$RNA@counts) * 100 

obj$percent.ig = percent.featureset

ncols <- c('Bm_CD27+_IGHM-' ="#E5D2DD", 
           'Bm_CD27+_IGHM+' ="#F3B1A0", 
           'Bm_CD27+_IGHM+_SOX5'="#E59CC4",
           'Bmn_IGHD' ="#AB3282",
           'Plasma'= "#8C549C")

p = Seurat::VlnPlot(object = obj,
                features = "percent.ig",
                pt.size = 0) +
  scale_fill_manual(values = ncols) +
  theme(legend.position = "none") +
  labs(title = "The proportion of IG gene counts", x = "", y = "") + #immunoglobulin
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 6),
    plot.title = element_text(size = 7,face = "plain"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    plot.margin = unit(c(0, 1, 0, 0), "char"),
  )

p$layers[[1]]$aes_params$size = 0.1
p

ggsave("figures/3_IG_gene_count.pdf",width = 2,height = 2)

unique(Bcell$Condition)

Bcell$Condition = factor(Bcell$Condition,levels = c("A0","A1","A2",
                                                    "B0","B1","B2",
                                                    "C0","C1","C2"))
obj <- Bcell
Idents(obj) = obj$Condition
features <- grep(pattern = "^IG", x = rownames(obj), value = TRUE) 
percent.featureset <- colSums(x = GetAssayData(object = obj, slot = "counts")[features, , drop = FALSE])/ 
  colSums(x = obj@assays$RNA@counts) * 100 

obj$percent.ig = percent.featureset

ncols <- rep("#8C549C", 9)

p = Seurat::VlnPlot(object = obj,
                features = "percent.ig",
                pt.size = 0) +
  #scale_fill_manual(values = ncols) +
  theme(legend.position = "none") +
  labs(title = "The proportion of IG gene counts", x = "", y = "") + #immunoglobulin
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 6),
    plot.title = element_text(size = 7,face = "plain"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    plot.margin = unit(c(0, 1, 0, 0), "char"),
  )

p$layers[[1]]$aes_params$size = 0.1
p

ggsave("figures/3_IG_Condition_ALL_gene_count.pdf",width = 2,height = 2)

unique(Bcell$ct_level_4)


Bcell$Condition = factor(Bcell$Condition,levels = c("A0","A1","A2",
                                                    "B0","B1","B2",
                                                    "C0","C1","C2"))

# ---- subset cell
Idents(Bcell) <- "ct_level_4"
obj <- subset(Bcell, idents = c("Plasma"))

Idents(obj) = obj$Condition
features <- grep(pattern = "^IG", x = rownames(obj), value = TRUE) 
percent.featureset <- colSums(x = GetAssayData(object = obj, slot = "counts")[features, , drop = FALSE])/ 
  colSums(x = obj@assays$RNA@counts) * 100 

obj$percent.ig = percent.featureset

ncols <- rep("#8C549C", 9)

p = Seurat::VlnPlot(object = obj,
                features = "percent.ig",
                pt.size = 0) +
  scale_fill_manual(values = ncols) +
  theme(legend.position = "none") +
  labs(title = "The proportion of IG gene counts", x = "", y = "") + #immunoglobulin
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 6),
    plot.title = element_text(size = 7,face = "plain"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    plot.margin = unit(c(0, 1, 0, 0), "char"),
  )

p$layers[[1]]$aes_params$size = 0.1
p

ggsave("figures/3_IG_Condition_Plasma_gene_count.pdf",width = 2,height = 2)

summary(Bcell$SHM_rate)

##--- Figure S2B; Heatmap showing the median SHM rate for each TIB major lineage 
tmp = Bcell@meta.data %>% 
   filter(! is.na("bcr_cdr3_nt"))%>% 
  group_by(Condition) %>%
  mutate(sample_n = length(unique(Sample))) %>%
  mutate(Condition_new = paste0(Condition," (N = ", sample_n,")")) %>%
  group_by(Condition_new, ct_level_4) %>%
  dplyr::summarise(avg_SHM = median(SHM_rate, na.rm = TRUE), B_n = n())

colnames(tmp)[1] = "Condition"  
tmp = tidyr::spread(data = tmp[, 
                               c("Condition", "ct_level_4", "avg_SHM")],
                    key = ct_level_4, value = avg_SHM)

matrix = tmp %>% tibble::column_to_rownames("Condition")
matrix = matrix[,c("Bmn_IGHD",  'Bm_CD27+_IGHM+','Bm_CD27+_IGHM+_SOX5','Bm_CD27+_IGHM-','Plasma')]
col_fun <- circlize::colorRamp2(c(2.812 , 3.549, 4.688),
                                c("#0F7B9F", "white", "#D83215"))
head(tmp)
head(matrix)

library(ComplexHeatmap)
pdf("figures/4_SHM_heatmap_Condition.pdf", width = 3, height = 4)
ComplexHeatmap::Heatmap(
  matrix,
  col = col_fun,
  na_col = "grey",
  cluster_rows = F,
  cluster_columns = F,
  row_names_side = "left",
  column_names_gp = grid::gpar(fontsize = 7),
  row_names_gp = grid::gpar(fontsize = 7),
  rect_gp = gpar(col = "white", lwd = 1),
  width = ncol(matrix) * unit(0.2, "inch"),
  height = nrow(matrix) * unit(0.2, "inch"),
  name = "Median SHM rate",
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 8)
  ))
dev.off()

library(dplyr)
library(ComplexHeatmap)
library(Startrac)   # 假设这是 STARTRAC 分析的包名
library(circlize)   # 用于 colorRamp2 的包

# Run STARTRAC
in.dat <-  Bcell@meta.data [, c("bcr_barcode", 
                                "bcr_clonotype_id_new", 
                                "Donor", "ct_level_4", "Condition")]%>% 
             mutate(Donor = as.character(Donor),ct_level_4 = as.character(ct_level_4))%>% 
             filter(! is.na("bcr_cdr3_nt"))
colnames(in.dat) <- c("Cell_Name", 'clone.id', 'patient', 'majorCluster', 'loc')
out <- Startrac.run(in.dat, proj = "Bcell", verbose = F)


# plot
tmp = out@pIndex.tran %>% filter(aid == "Bcell")%>% as.data.frame()

head(tmp)
rownames(tmp) = tmp$majorCluster
tmp = subset(tmp, select = -c(majorCluster, aid,NCells))
m = as.matrix(tmp)

summary(m)
col_fun <- colorRamp2(c(0, 0.0065, 0.01, 0.05),
                      c("#fffde7", "#ffe0b2", "#ff9800", "#e65100"))
pdf("figures/4_STARTRAC_res.pdf", width = 5,height = 3)
Heatmap(
  m,
  col = col_fun,
  cluster_rows = T,
  cluster_columns = T,
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  width = ncol(m) * unit(0.2, "inch"),
  height = nrow(m) * unit(0.2, "inch"),
  name = "Number of overlaped clones",
)
dev.off()









