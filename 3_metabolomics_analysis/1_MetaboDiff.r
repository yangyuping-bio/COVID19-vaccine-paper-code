#library("devtools")
#install_github("andreasmock/MetaboDiff")
library(MetaboDiff)
library(readxl)

getwd()

colData <- read_excel("2.MetaPro-13-21_20211109-1.xlsx", sheet = "Sample Info")
#assay <- read_excel("2.MetaPro-13-21_20211109-1.xlsx", sheet = "Normalized data")
assay <-  read_excel("2.MetaPro-13-21_20211109-1.xlsx", sheet = "Peak Area Data")

head(assay)

## 去掉Food Component/Plant和Xenobiotics
nrow(assay)
assay <- as.data.frame(assay)%>% 
        filter(superPathway != 'Xenobiotics')%>% 
        filter(superPathway != 'Food Component/Plant')
nrow(assay)
length(unique(assay$superPathway))
length(unique(assay$subPathway))

rowData <- assay
rowData <-rowData[,1:11]
colnames(rowData) 
rowData <- as.data.frame(rowData)
rownames(rowData) <- rowData$name
write.csv(rowData, "data/1_Metabbolism_rowData.csv")
head(rowData, n=2)

colData <- as.data.frame(colData)
rownames(colData) <- colData$"SAMPLE ID"
colData$adaptive_group<- colData$"Group ID"
 colData$adaptive_group[colData$adaptive_group %in% c('A1','A2','A0','B0')] <- 'before_B1'
 colData$adaptive_group[colData$adaptive_group %in% c('B1','B2')] <- 'B1_B2'
# colData$adaptive_group[colData$adaptive_group %in% c('B0')] <- 'B0'
# colData$adaptive_group[colData$adaptive_group %in% c('B1')] <- 'B1'
#colData$adaptive_group[colData$adaptive_group %in% c('A1','A2','A0','B2')] <- NA
colData$innate_group <- colData$"Group ID"
colData$innate_group[colData$innate_group %in% c('A1')] <- 'A1'
colData$innate_group[colData$innate_group =='A0'] <- 'Baseline'
colData$innate_group[colData$innate_group %in% c('B1','B2','B0','A2')] <- NA
write.csv(colData, "data/1_Metabbolism_colData.csv")
colData

assay <- as.data.frame(assay)
colnames(assay)
colnames(assay)[c(2,12:44)]
assay <- assay[ , c(2,12:44)]
rownames(assay) <- assay$name
assay$name <- NULL
write.csv(assay, "data/1_Metabbolism_assay_Data.csv")
head(assay,n=2)

nrow(assay)
nrow(colData)
nrow(rowData)
#将三个数据集融合成一个以便于下游分析。
met <- create_mae(assay,rowData,colData)
met

db = read.csv("../0_data/All_smpdb_metabolites_2024_07_02.csv")
nrow(db)
head(db,n=1)

get_SMPDBanno_2 <- function(met, column_kegg_id, column_hmdb_id, column_cas_id) {
    rowData = rowData(met[["raw"]])
    db = read.csv("../0_data/All_smpdb_metabolites_2024_07_02.csv")
    
    # 初始化结果矩阵，包含SMPDB信息
    res = matrix(NA, nrow = nrow(rowData), ncol = ncol(db))
    colnames(res) = paste0("SMPDB|", colnames(db))

    # 精确匹配的函数，优先匹配最相似的项
    find_best_match <- function(row_value, db_column) {
        possible_matches = which(db_column == row_value)
        if (length(possible_matches) == 1) {
            return(db[possible_matches, , drop = FALSE])  # 使用 drop = FALSE 保持维度
        } else if (length(possible_matches) > 1) {
            # 当有多个匹配项时，选择第一个匹配项（或定义更复杂的选择逻辑）
            best_match = possible_matches[1]
            return(db[best_match, , drop = FALSE])
        } else {
            return(matrix(NA, nrow = 1, ncol = ncol(db)))
        }
    }

    # 逐步匹配 KEGG, HMDB, CAS
    if (!is.na(column_kegg_id)) {
        for (i in 1:nrow(rowData)) {
            match_result = find_best_match(rowData[i, column_kegg_id], db$KEGG.ID)
            if (!all(is.na(match_result))) {
                res[i, ] = as.matrix(match_result)  # 确保 match_result 是一个矩阵
            }
        }
    }
    
    if (!is.na(column_hmdb_id)) {
        for (i in 1:nrow(rowData)) {
            if (all(is.na(res[i, ]))) {  # 只在没有匹配到时再继续匹配
                match_result = find_best_match(rowData[i, column_hmdb_id], db$HMDB.ID)
                if (!all(is.na(match_result))) {
                    res[i, ] = as.matrix(match_result)  # 确保 match_result 是一个矩阵
                }
            }
        }
    }

    if (!is.na(column_cas_id)) {
        for (i in 1:nrow(rowData)) {
            if (all(is.na(res[i, ]))) {  # 只在没有匹配到时再继续匹配
                match_result = find_best_match(rowData[i, column_cas_id], db$CAS)
                if (!all(is.na(match_result))) {
                    res[i, ] = as.matrix(match_result)  # 确保 match_result 是一个矩阵
                }
            }
        }
    }

    # 更新 rowData，添加匹配到的 SMPDB 注释信息
    rowData(met[["raw"]]) = data.frame(rowData, res)
    
    return(met)
}


met <- get_SMPDBanno_2(met,
                           column_kegg_id=9,
                           column_hmdb_id=10,
                          column_cas_id= 8
                    )


SMPDB_res <-  rowData(met[["raw"]])
head(SMPDB_res)
write.csv(SMPDB_res, "data/2_SMPDB_res.csv")

pdf('figures/2_na_heatmap_group.pdf',10,10)
print(na_heatmap(met,
           group_factor="Group ID",
            label_colors=c("#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF", "#BFB2FFFF")))
dev.off()

#https://www.jianshu.com/p/80af83c2f630
#剔除缺失值，计算代谢物的相对丰度。
met = knn_impute(met,cutoff=0.25)

met

pdf('figures/2_outlier_heatmap.pdf',10,10)
print(outlier_heatmap(met,
                  group_factor="Group ID",
              label_colors=c("#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF", "#BFB2FFFF"),
                 k=5))
dev.off()

#如果有异常值
#met <- remove_cluster(met,cluster=1)

met <- normalize_met(met)

met

pdf('figures/2_quality_plot_adaptive_group.pdf',10,15)
quality_plot(met, group_factor="adaptive_group",
              label_colors=c("#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF", "#BFB2FFFF"))
dev.off()

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

label_colors <- c("#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF", "#BFB2FFFF")
pdf('figures/3_multiplot_pac_tsne_adaptive.pdf',12,5)
multiplot(
   pca_plot(met,
            group_factor="adaptive_group",
            label_colors=label_colors),
   tsne_plot(met,
             group_factor="adaptive_group",
             label_colors=label_colors),
   cols=2)
dev.off()

label_colors <- c("#B22C2CFF", "#85B22CFF", "#E5B17EFF", "#51A3CCFF", "#E57E7EFF", "#BFB2FFFF")
pdf('figures/3_multiplot_pac_tsne_innate.pdf',12,5)
multiplot(
   pca_plot(met,
            group_factor="innate_group",
            label_colors=label_colors),
   tsne_plot(met,
             group_factor="innate_group",
             label_colors=label_colors),
   cols=2)

dev.off()

#对单个代谢物进行差异分析，主要用T检验和ANOVA分析。
met = diff_test(met,
                 group_factors = c("adaptive_group","innate_group","Group ID")) #,"innate_group","Group ID"
 str(metadata(met), max.level=2)


volcano_plot <- function(met, group_factor, label_colors, dm_cutoff=0.5, p_adjust=TRUE, ...) {

    id = grep(group_factor, names(metadata(met)))[1]
    df = metadata(met)[[id]]
    name = names(metadata(met))[id]
    lv = levels(as.factor(colData(met)[[group_factor]]))
    xlabel = paste0("difference in means"," [",paste(lv, collapse = "-"),"]")

    if (p_adjust == TRUE) {
        plot(df$dm, -1*log10(df$adj_pval), bty="n", pch=20,
             xlab=xlabel, ylab="adjusted -log10 p-value", ...)
        abline(h=-log10(0.05), lty=2)
        abline(v=dm_cutoff, lty=2)
        abline(v=-dm_cutoff, lty=2)
        
        Tsig = df$dm < (-dm_cutoff) & df$adj_pval < 0.05
        points(x = df$dm[Tsig], y = -1*log10(df$adj_pval)[Tsig], col = label_colors[1], pch=20)
        text(x = df$dm[Tsig], y = -1*log10(df$adj_pval)[Tsig], labels = df$metabolite[Tsig], pos = 4, cex = 1.5)
        
        Nsig = df$dm > dm_cutoff & df$adj_pval < 0.05
        points(x = df$dm[Nsig], y = -1*log10(df$adj_pval)[Nsig], col = label_colors[2], pch=20)
        text(x = df$dm[Nsig], y = -1*log10(df$adj_pval)[Nsig], labels = df$metabolite[Nsig], pos = 4, cex = 1.5)
        
    } else {
        plot(df$dm, -1*log10(df$pval), bty="n", pch=20,
             xlab=xlabel, ylab="-log10 p-value", ...)
        abline(h=-log10(0.05), lty=2)
        abline(v=dm_cutoff, lty=2)
        abline(v=-dm_cutoff, lty=2)
        
        Tsig = df$dm < (-dm_cutoff) & df$pval < 0.05
        points(x = df$dm[Tsig], y = -1*log10(df$pval)[Tsig], col = label_colors[1], pch=20)
        text(x = df$dm[Tsig], y = -1*log10(df$pval)[Tsig], labels = df$metabolite[Tsig], pos = 4, cex = 1.5)
        
        Nsig = df$dm > dm_cutoff & df$pval < 0.05
        points(x = df$dm[Nsig], y = -1*log10(df$pval)[Nsig], col = label_colors[2], pch=20)
        text(x = df$dm[Nsig], y = -1*log10(df$pval)[Nsig], labels = df$metabolite[Nsig], pos = 4, cex = 1.5)
    }
}

# Point Plot of Cytokines ---------------------------

cyto.mhu.bind <- data.frame(row.names = colnames(df.cyto))
for (i in c(1:2)) {
  cyto.mhu <- read.xlsx("Serum.cytokine.MannWhitneyU.xlsx", # Source data
                        sheetIndex = i, startRow = 1, as.data.frame = T, header = T, check.names = F)
  cyto.mhu$Group <- group.name[i]
  cyto.mhu.bind <- rbind(cyto.mhu.bind, 
                         cyto.mhu)
}
cyto.mhu.bind$type <- ifelse(cyto.mhu.bind$BH < 0.05, "Significant", "Not Significant")
cyto.mhu.bind$type <- factor(cyto.mhu.bind$type, levels = c("Significant", "Not Significant"))

cyto_all <- data.frame(row.names = unique(cyto.mhu.bind$cytokine), 
                       MH_log2FC = cyto.mhu.bind$LOG2FC[match(unique(cyto.mhu.bind$cytokine), 
                                                              cyto.mhu.bind$cytokine[cyto.mhu.bind$Group == "MvsH"])], 
                       SH_log2FC = cyto.mhu.bind$LOG2FC[match(unique(cyto.mhu.bind$cytokine), 
                                                              cyto.mhu.bind$cytokine[cyto.mhu.bind$Group == "MvsH"])])
# sort the cytokines
cyto.ord <- hclust(dist(cyto_all, method = "euclidean"), method = "average")
cyto.ord <- rownames(cyto_all)[cyto.ord$order]
cyto.mhu.bind$cytokine <- factor(cyto.mhu.bind$cytokine, levels = rev(cyto.ord))

# point plot
ggplot(cyto.mhu.bind, aes(Group, cytokine)) +
  geom_point(aes(fill = LOG2FC, size = -log10(BH), color = type, shape = type)) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#3d67a3", "white", "#ce1020"))(99)[-c(47:49, 50, 51:53)], 
                       limits = c(-6,6), breaks = c(-6, -3, 0, 3, 6)) +
  scale_shape_manual(values = c(21, 23)) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_size_continuous(range = c(1, 3.5), breaks = c(1, 3, 5, 7)) +
  theme_dendro() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(fill = "Log2 (Fold change)", size = "-Log10 (P value)", shape = "", color = "") +
  theme(axis.text.x = element_text(size = 8, face = "plain", colour = "black", angle = 45), 
        axis.text.y = element_text(size = 8, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 8, face = "plain", colour = "black"), 
        legend.title = element_text(size = 8, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))   

pdf('figures/4_volcano_plot_adaptive.pdf',20,15)
#进行差异分析
volcano_plot(met, 
              group_factor="adaptive_group",
              label_colors=label_colors,
              p_adjust = FALSE)
dev.off()
pdf('figures/4_volcano_plot_innate.pdf',20,15)
volcano_plot(met, 
              group_factor="innate_group",
              label_colors=label_colors,
              p_adjust = FALSE)
dev.off()

write.csv(metadata(met)[[1]], "data/adaptive_diff_test.csv")
write.csv(metadata(met)[[2]], "data/innate_diff_test.csv")

rowData(met[["norm_imputed"]])$name[metadata(met)$`ttest_innate_group_Baseline_vs_A1`$pval<0.05]
# 获取代谢物名称
metabolites_of_interest <- rowData(met[["norm_imputed"]])$name[metadata(met)$`ttest_innate_group_Baseline_vs_A1`$pval<0.05]
length(metabolites_of_interest)
# 获取代谢物在 norm_imputed 数据中的行号
row_indices <- which(rowData(met[["norm_imputed"]])$name %in% metabolites_of_interest)
row_indices

met$innate_group

met$adaptive_group

pdf('figures/5_metabolic_barplot_innate_group_p_lower_than_0.05.pdf', 30, 15)
ids = met$innate_group %in% c("Baseline", "A1_A2")
group_factor = factor(met$innate_group[ids], levels=c("Baseline", "A1_A2"))  # 确保转换为因子

par(mfrow=c(4, length(row_indices)/4))  # 设置绘图窗口大小

for (i in seq_along(row_indices)) {
  # 获取代谢物的行号
  idx <- row_indices[i]
  
  # 获取数据
  y_data <- assay(met[["norm_imputed"]])[idx, ids]
  
  # 绘制图形
  plot(group_factor, y_data,
       main = rowData(met[["norm_imputed"]])$name[idx],  # 代谢物名称作为标题
       xlab = "", ylab = "Normalized values", 
       frame = FALSE, col = c("orange", "darkseagreen"))
  
  # 进行t检验
  t_test_result <- t.test(y_data ~ group_factor)
  p_value <- t_test_result$p.value
  
  # 选择统计显著性标记
  if (p_value < 0.001) {
    sign_marker <- "***"
  } else if (p_value < 0.01) {
    sign_marker <- "**"
  } else if (p_value < 0.05) {
    sign_marker <- "*"
  } else {
    sign_marker <- ""
  }
  
  # 确定标记位置
  y_max <- max(y_data)
  y_offset <- 0.001 * (y_max - min(y_data))  # 设定标记的偏移量
  
  # 添加统计显著性标记
  if (sign_marker != "") {
    text(x = 2, y = y_max - y_offset, sign_marker, cex = 2, col = "red")
  }
  
  # 打印p值
  print(p_value)
}
dev.off()


rowData(met[["norm_imputed"]])$name[metadata(met)$`ttest_adaptive_group_before_B1_vs_B1_B2`$pval<0.05]
# 获取代谢物名称
metabolites_of_interest <- rowData(met[["norm_imputed"]])$name[metadata(met)$`ttest_adaptive_group_before_B1_vs_B1_B2`$pval<0.05]
length(metabolites_of_interest)
# 获取代谢物在 norm_imputed 数据中的行号
row_indices <- which(rowData(met[["norm_imputed"]])$name %in% metabolites_of_interest)
row_indices

pdf('figures/5_metabolic_barplot_adaptive_group_p_lower_than_0.05.pdf', 20, 10)
ids = met$adaptive_group %in% c("before_B1", "B1_B2")
group_factor = factor(met$adaptive_group[ids], levels=c("before_B1", "B1_B2"))  # 确保转换为因子

par(mfrow=c(3, length(row_indices)/3))  # 设置绘图窗口大小

for (i in seq_along(row_indices)) {
  # 获取代谢物的行号
  idx <- row_indices[i]
  
  # 获取数据
  y_data <- assay(met[["norm_imputed"]])[idx, ids]
  
  # 绘制图形
  plot(group_factor, y_data,
       main = rowData(met[["norm_imputed"]])$name[idx],  # 代谢物名称作为标题
       xlab = "", ylab = "Normalized values", 
       frame = FALSE, col = c("orange", "darkseagreen"))
  
  # 进行t检验
  t_test_result <- t.test(y_data ~ group_factor)
  p_value <- t_test_result$p.value
  
  # 选择统计显著性标记
  if (p_value < 0.001) {
    sign_marker <- "***"
  } else if (p_value < 0.01) {
    sign_marker <- "**"
  } else if (p_value < 0.05) {
    sign_marker <- "*"
  } else {
    sign_marker <- ""
  }
  
  # 确定标记位置
  y_max <- max(y_data)
  y_offset <- 0.001 * (y_max - min(y_data))  # 设定标记的偏移量
  
  # 添加统计显著性标记
  if (sign_marker != "") {
    text(x = 2, y = y_max - y_offset, sign_marker, cex = 2, col = "red")
  }
  
  # 打印p值
  print(p_value)
}
dev.off()


names(metadata(met))

rowData(met[["norm_imputed"]])$name[metadata(met)$`anova_Group ID_A0_vs_A1_vs_A2_vs_B0_vs_B1_vs_B2`$pval<0.05]
# 获取代谢物名称
metabolites_of_interest <- rowData(met[["norm_imputed"]])$name[metadata(met)$`anova_Group ID_A0_vs_A1_vs_A2_vs_B0_vs_B1_vs_B2`$pval<0.05]
length(metabolites_of_interest)
# 获取代谢物在 norm_imputed 数据中的行号
row_indices <- which(rowData(met[["norm_imputed"]])$name %in% metabolites_of_interest)
row_indices

metabolites_of_interest <- c("arginine","choline","7-methylguanine",
                             "chenodeoxycholate","glucose","ethyl alpha-glucopyranoside",
                             "glycerol","Oxidative Phosphorylation","kynurenine")


##精氨酸
# 获取代谢物名称
# metabolites_of_interest <- c("2-oxoarginine*","argininate*","arginine",
#                              "dimethylarginine (SDMA + ADMA)",
#                              "homoarginine","N-acetylarginine")
metabolites_of_interest <- c("arginine","choline","7-methylguanine",
                             "chenodeoxycholate","glucose","ethyl alpha-glucopyranoside",
                             "glycerol","Oxidative Phosphorylation","kynurenine")

length(metabolites_of_interest)
# 获取代谢物在 norm_imputed 数据中的行号
row_indices <- which(rowData(met[["norm_imputed"]])$name %in% metabolites_of_interest)
row_indices

#pdf('figures/5_metabolic_barplot_arginine_0.05.pdf', 8, 6)
pdf('figures/5_metabolic_barplot_compass_cor_0.05.pdf', 10, 6)
ids = met$adaptive_group %in% c("before_B1", "B1_B2")
group_factor = factor(met$adaptive_group[ids], levels=c("before_B1", "B1_B2"))  # 确保转换为因子

par(mfrow=c(2, length(row_indices)/2))  # 设置绘图窗口大小

for (i in seq_along(row_indices)) {
  # 获取代谢物的行号
  idx <- row_indices[i]
  
  # 获取数据
  y_data <- assay(met[["norm_imputed"]])[idx, ids]
  
  # 绘制图形
  plot(group_factor, y_data,
       main = rowData(met[["norm_imputed"]])$name[idx],  # 代谢物名称作为标题
       xlab = "", ylab = "Normalized values", 
       frame = FALSE, col = c("orange", "darkseagreen"))
  
  # 进行t检验
  t_test_result <- t.test(y_data ~ group_factor)
  p_value <- t_test_result$p.value
  
  # 选择统计显著性标记
  if (p_value < 0.001) {
    sign_marker <- "***"
  } else if (p_value < 0.01) {
    sign_marker <- "**"
  } else if (p_value < 0.05) {
    sign_marker <- "*"
  } else {
    sign_marker <- ""
  }
  
  # 确定标记位置
  y_max <- max(y_data)
  y_offset <- 0.001 * (y_max - min(y_data))  # 设定标记的偏移量
  
  # 添加统计显著性标记
  if (sign_marker != "") {
    text(x = 2, y = y_max - y_offset, sign_marker, cex = 2, col = "red")
  }
  
  # 打印p值
  print(p_value)
}


ids = met$innate_group %in% c("Baseline", "A1_A2")
group_factor = factor(met$innate_group[ids], levels=c("Baseline", "A1_A2"))  # 确保转换为因子
par(mfrow=c(2, length(row_indices)/2))  # 设置绘图窗口大小
for (i in seq_along(row_indices)) {
  # 获取代谢物的行号
  idx <- row_indices[i]
  
  # 获取数据
  y_data <- assay(met[["norm_imputed"]])[idx, ids]
  
  # 绘制图形
  plot(group_factor, y_data,
       main = rowData(met[["norm_imputed"]])$name[idx],  # 代谢物名称作为标题
       xlab = "", ylab = "Normalized values", 
       frame = FALSE, col = c("orange", "darkseagreen"))
  
  # 进行t检验
  t_test_result <- t.test(y_data ~ group_factor)
  p_value <- t_test_result$p.value
  
  # 选择统计显著性标记
  if (p_value < 0.001) {
    sign_marker <- "***"
  } else if (p_value < 0.01) {
    sign_marker <- "**"
  } else if (p_value < 0.05) {
    sign_marker <- "*"
  } else {
    sign_marker <- ""
  }
  
  # 确定标记位置
  y_max <- max(y_data)
  y_offset <- 0.001 * (y_max - min(y_data))  # 设定标记的偏移量
  
  # 添加统计显著性标记
  if (sign_marker != "") {
    text(x = 2, y = y_max - y_offset, sign_marker, cex = 2, col = "red")
  }
  
  # 打印p值
  print(p_value)
}


# 确保组因子是有序的因子
ids = met$'Group ID' %in% c("A0", "A1", "A2", "B0", "B1", "B2")
group_factor = factor(met$'Group ID'[ids], levels = c("A0", "A1", "A2", "B0", "B1", "B2"))

# 设置绘图窗口大小
par(mfrow = c(2, ceiling(length(row_indices) / 2)))

for (i in seq_along(row_indices)) {
  # 获取代谢物的行号
  idx <- row_indices[i]
  
  # 获取数据
  y_data <- assay(met[["norm_imputed"]])[idx, ids]
  
  # 绘制图形
  plot(group_factor, y_data,
       main = rowData(met[["norm_imputed"]])$name[idx],  # 代谢物名称作为标题
       xlab = "", ylab = "Normalized values", 
       frame = FALSE, col = c("orange", "darkseagreen", "red", "blue", "purple", "green"))
  
  # 进行ANOVA
  anova_result <- aov(y_data ~ group_factor)
  p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]
  
  # 选择统计显著性标记
  if (p_value < 0.001) {
    sign_marker <- "***"
  } else if (p_value < 0.01) {
    sign_marker <- "**"
  } else if (p_value < 0.05) {
    sign_marker <- "*"
  } else {
    sign_marker <- ""
  }
  
  # 确定标记位置
  y_max <- max(y_data, na.rm = TRUE)
  y_offset <- 0.001 * (y_max - min(y_data, na.rm = TRUE))  # 设定标记的偏移量
  
  # 添加统计显著性标记
  if (sign_marker != "") {
    text(x = 3, y = y_max - y_offset, sign_marker, cex = 2, col = "red")
  }
  
  # 打印p值
  print(p_value)
}

dev.off()


pdf('figures/5_metabolic_barplot_arginine_outline_0.05.pdf', 8, 6)
# 确保组因子是有序的因子
ids = met$'Group ID' %in% c("A0", "A1", "A2", "B0", "B1", "B2")
group_factor = factor(met$'Group ID'[ids], levels = c("A0", "A1", "A2", "B0", "B1", "B2"))

# 设置绘图窗口大小
par(mfrow = c(2, ceiling(length(row_indices) / 2)))

for (i in seq_along(row_indices)) {
  # 获取代谢物的行号
  idx <- row_indices[i]
  
  # 获取数据
  y_data <- assay(met[["norm_imputed"]])[idx, ids]
  
  # 绘制箱线图，不显示异常值
  boxplot(y_data ~ group_factor, 
          main = rowData(met[["norm_imputed"]])$name[idx],  # 代谢物名称作为标题
          xlab = "", ylab = "Normalized values", 
          frame = FALSE, col = c("orange", "darkseagreen", "red", "blue", "purple", "green"), 
          outline = FALSE)
  
  # 计算新的中位数（不包括异常值）
  median_values <- tapply(y_data, group_factor, function(x) median(x, na.rm = TRUE))
  
  # 重新绘制中位线
  for (j in seq_along(median_values)) {
    segments(j - 0.4, median_values[j], j + 0.4, median_values[j], lwd = 2, col = "black")
  }
  
  # 进行ANOVA
  anova_result <- aov(y_data ~ group_factor)
  
  # 进行Tukey's HSD检验
  tukey_result <- TukeyHSD(anova_result)
  
  # 获取p值并选择统计显著性标记
  comparison_labels <- rownames(tukey_result$group_factor)
  significant_comparisons <- tukey_result$group_factor[,"p adj"] < 0.05
  
  # 在图中添加统计显著性标记
  y_max <- max(y_data[!is.na(y_data) & y_data != Inf], na.rm = TRUE)
  y_offset <- 0.5 * (y_max - min(y_data[!is.na(y_data) & y_data != -Inf], na.rm = TRUE))  # 设定标记的偏移量
  text_y_pos <- y_max - y_offset
  
  if (any(significant_comparisons)) {
    for (comp_idx in which(significant_comparisons)) {
      comp_label <- comparison_labels[comp_idx]
      comp_p_value <- tukey_result$group_factor[comp_idx, "p adj"]
      
      if (comp_p_value < 0.001) {
        sign_marker <- "***"
      } else if (comp_p_value < 0.01) {
        sign_marker <- "**"
      } else if (comp_p_value < 0.05) {
        sign_marker <- "*"
      } else {
        sign_marker <- ""
      }
      
      if (sign_marker != "") {
        text(x = 3, y = text_y_pos, labels = paste(comp_label, sign_marker), cex = 0.8, col = "red")
        text_y_pos <- text_y_pos + y_offset  # 更新文本位置以避免重叠
      }
    }
  }
  
  # 打印 Tukey's HSD p值
  print(tukey_result)
}

dev.off()

pdf('figures/5_matabolic_barplot.pdf',9,4)
ids = met$adaptive_group %in% c("before_B1", "B1_B2")
group_factor = factor(met$adaptive_group[ids],levels=c("before_B1", "B1_B2"))  # 确保转换为因子

par(mfrow=c(1,3))
#1
plot(assay(met[["norm_imputed"]])[1, ids] ~ group_factor,
     main = "beta-citrylglutamate", xlab = "", ylab = "Normalized values", 
     frame = FALSE, col = c("orange", "darkseagreen"))
text(x=2, y=17.55, "**", cex=2)

t.test(assay(met[["norm_imputed"]])[1, ids] ~ group_factor)[[3]]

#2
plot(assay(met[["norm_imputed"]])[2, ids] ~ group_factor,
     main = "5-oxoproline", xlab = "", ylab = "Normalized values", 
     frame = FALSE, col = c("orange", "darkseagreen"))
text(x=2, y=25.2, "*", cex=2)
t.test(assay(met[["norm_imputed"]])[2, ids] ~ group_factor)[[3]]

#3
plot(assay(met[["norm_imputed"]])[3, ids] ~ group_factor,
     main = "N-acetylthreonine", xlab = "", ylab = "Normalized values", 
     frame = FALSE, col = c("orange", "darkseagreen"))
text(x=2,y=22.62,"**",cex=2)
t.test(assay(met[["norm_imputed"]])[3, ids] ~ group_factor)[[3]]
dev.off()

#install.packages("VennDiagram")
library(VennDiagram)

rowData(met[["norm_imputed"]])$name[metadata(met)$`ttest_adaptive_group_before_B1_vs_B1_B2`$pval<0.05]

A = as.character(rowData(met[["norm_imputed"]])$name[metadata(met)$`ttest_adaptive_group_before_B1_vs_B1_B2`$pval<0.05])
B = as.character(rowData(met[["norm_imputed"]])$name[metadata(met)$`ttest_innate_group_Baseline_vs_A1_A2`$pval<0.05])
C = as.character(rowData(met[["norm_imputed"]])$name[metadata(met)$`anova_Group ID_A0_vs_A1_vs_A2_vs_B0_vs_B1_vs_B2`$pval<0.05])
pdf("figures/6_venn.pdf",width=5,height=5)
venn <- draw.triple.venn(area1 = length(A),
                         area2 = length(B),
                         area3 = length(C),
                         n12 = length(intersect(A,B)),
                         n13 = length(intersect(A,C)),
                         n23 = length(intersect(B,C)),
                         n123 = length(intersect(intersect(A,B),C)),
                         fill = c("dodgerblue","red3","yellow"),
                         alpha=c(0.1,0.1,0.1),#颜色的透明度
                         category = c("adaptive","innate","Condition"),
                         lwd = c(0.5,0.5,0.5),#维恩图的线条宽度
                         cex = 1.3,#维恩图中数字的大小
                         fontfamily = "sans",#字体样式
                         cat.cex = 1.3 #标签大小
                         )
grid.draw(venn)
dev.off()

met_example <- met %>% #met_example
   diss_matrix %>%    #构建相异矩阵
   identify_modules(min_module_size=6) %>%  #鉴定代谢相关模块
   name_modules(pathway_annotation="subPathway") %>%  #代谢相关模块命名
   calculate_MS(group_factors=c("adaptive_group","innate_group","Group ID")) #根据样本性状计算模块之间关联的显著性

#相关网络的属性
table(metadata(met_example)$modules)


#代谢相关模块可视化，分级聚类
pdf("figures/7_moudle_min_6.pdf",width=15,height=10)
 WGCNA::plotDendroAndColors(metadata(met_example)$tree, 
                            metadata(met_example)$module_color_vector, 
                            'Module colors', 
                            dendroLabels = FALSE, 
                           hang = 0.03,
                           addGuide = TRUE, 
                            guideHang = 0.05, main='')
dev.off()

#相关网络的属性
table(metadata(met_example)$modules)

MS_plot = function(met, group_factor,p_value_cutoff,p_adjust=FALSE){
    id = grep(group_factor,names(metadata(met)))[2]
    tree = ape::as.phylo(metadata(met)$METree)
    if(p_adjust==TRUE){
        x = -log10(metadata(met)[[id]]$av_adj_pval)
    } else {
        x = -log10(metadata(met)[[id]]$av_pval)
    }

    names(x) = tree$tip.label
    obj=phytools::contMap(tree=tree,
                      x=x,
                      res=400,
                      plot=FALSE,
                      lims=c(0,-log10(p_value_cutoff))
                      )
    obj = phytools::setMap(obj, colors=c(rep("white",9),"brown"))
    plot(obj)
}

pdf("figures/7_MS_plot_min_6.pdf",width=18,height=20)
MS_plot(met_example,
       group_factor="adaptive_group",
       p_value_cutoff=0.1,
       p_adjust=FALSE
)
MS_plot(met_example,
       group_factor="innate_group",
       p_value_cutoff=0.1,
       p_adjust=FALSE
)
MS_plot(met_example,
       group_factor="Group ID",
       p_value_cutoff=0.1,
       p_adjust=FALSE
)
dev.off()

# 2-oxoarginine*
# argininate*
# arginine
# dimethylarginine (SDMA + ADMA)
# homoarginine
# N-acetylarginine

# 设定要查询的模块编号（假设为 1），39是adptive
module_of_interest <- 86
modules_info <- metadata(met_example)$modules
metabolites_in_module <- which(modules_info == module_of_interest)
metabolite_names <- rowData(met[["norm_imputed"]])$name[metabolites_in_module]
print(metabolite_names)

# 设定要查询的模块编号（假设为 1），39是adptive
module_of_interest <- 39
modules_info <- metadata(met_example)$modules
metabolites_in_module <- which(modules_info == module_of_interest)
metabolite_names <- rowData(met[["norm_imputed"]])$name[metabolites_in_module]
print(metabolite_names)

# 设定要查询的模块编号（假设为 1）,24和48是innate
module_of_interest <- 24
modules_info <- metadata(met_example)$modules
metabolites_in_module <- which(modules_info == module_of_interest)
metabolite_names <- rowData(met[["norm_imputed"]])$name[metabolites_in_module]
print(metabolite_names)
module_of_interest <- 48
modules_info <- metadata(met_example)$modules
metabolites_in_module <- which(modules_info == module_of_interest)
metabolite_names <- rowData(met[["norm_imputed"]])$name[metabolites_in_module]
print(metabolite_names)

#代谢相关模块可视化，各模块直接的关系
par(mar=c(2,2,2,2))
 ape::plot.phylo(ape::as.phylo(metadata(met_example)$METree),
                type = 'fan',
                 show.tip.label = FALSE, 
                main='')
ape::tiplabels(frame = 'circle',
                col='black', 
                text=rep('',length(unique(metadata(met_example)$modules))), 
               bg = WGCNA::labels2colors(0:21))

# # 检查 tree$tip.label
# tree = ape::as.phylo(metadata(met_example)$METree)
# print(tree$tip.label)


MOI_plot_2 = function(met, group_factor, MOI, label_colors, p_adjust=TRUE, ...){
    id = grep(group_factor,names(metadata(met)))[1]
    mets = rowData(met[["norm_imputed"]])$name[metadata(met)$modules==MOI]
    if(p_adjust==TRUE){
        x = -log10(metadata(met)[[id]]$adj_pval[metadata(met)$modules==MOI])
        y =  metadata(met)$MM[,MOI+1][metadata(met)$modules==MOI]
        fc = metadata(met)[[id]]$dm[metadata(met)$modules==MOI]
        df = data.frame(mets=mets,x=x,y=y,fc=fc)
        ggplot(df, aes(x=x,y=y,color=fc)) +
            geom_point() +
            scale_color_gradient2(low=label_colors[1],mid="grey",high=label_colors[2],midpoint=0) +
            geom_text(aes(label=mets), vjust=1.5) +
            theme_classic() +
            xlab("adjusted p-value (-log10)") +
            ylab("module membership") +
            xlim(range(x) + c(-0.5,+1)) +
            geom_vline(aes(xintercept = 1.30103),
                       colour="darkorange3",
                       linetype="dashed")
    } else {
        x = -log10(metadata(met)[[id]]$pval[metadata(met)$modules==MOI])
        y =  metadata(met)$MM[,MOI+1][metadata(met)$modules==MOI]
        fc = metadata(met)[[id]]$dm[metadata(met)$modules==MOI]
        df = data.frame(mets=mets,x=x,y=y,fc=fc)
        ggplot(df, aes(x=x,y=y,color=fc)) +
            geom_point() +
            scale_color_gradient2(low=label_colors[1],mid="black",high=label_colors[2],midpoint=0) +
            geom_text(aes(label=mets), vjust=1.5) +
            theme_classic() +
            xlab("p-value (-log10)") +
            ylab("module membership") +
            xlim(range(x) + c(-0.5,+1)) +
            geom_vline(aes(xintercept = 1.30103),
                       colour="darkorange3",
                       linetype="dashed")
    }
}


#相关网络的属性
table(metadata(met_example)$modules)
num_MOI <- length(unique(metadata(met_example)$modules))
num_MOI

pdf("figures/7_MOI_plot_adaptive_group.pdf", width=8, height=8)
num_MOI <- length(unique(metadata(met_example)$modules))-1
# 循环遍历每个 MOI 并绘制图表
for (i in 0:num_MOI) {
    print(MOI_plot_2(met_example,
               group_factor="adaptive_group",
               MOI = i,
               label_colors=c("darkseagreen","dodgerblue"),
               p_adjust = FALSE))
}

dev.off()

pdf("figures/7_Hierarchy_of_metabolic_correlation_modules.pdf",width=12,height=18)
# plot phylogram with names
ape::plot.phylo(ape::as.phylo(metadata(met_example)$METree), cex=0.9)
dev.off()














