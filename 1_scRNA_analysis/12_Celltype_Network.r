library(Hmisc)
library(psych)
library(corrplot)
library(network)
library(ggnetwork)
library("tidyverse")




#All_Cell <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/5_new_level_1_2/3_All_Cell_CCA_integrated_by_Sample.rds")

# metadata.df <- All_Cell@meta.data
# write.csv(metadata.df ,'data/All_cell_metadata_df.csv')

metadata.df  <- read.csv('data/All_cell_metadata_df.csv')

colnames(metadata.df)

plot_cellprop_corr <- function(ratio.df, condition_seq, ct_level, 
                               output_pdf,title_name,width = 5.5, height  = 4,
                                color_seq = c("pos"= "red", 'neg' = "blue",
                                  Bcell="#AB3282",  # B细胞类型
                                  CD4T="#6778AE", # CD4 T细胞类型
                                  CD8T="#53A85F", # CD8 T细胞类型
                                  cDC="#9FA3A8",  # CDC细胞类型
                                  gdT="#23452F",  #gdT
                                  HPSC="#585658",#HPSC
                                  ILC="#F1BC62", #ILC
                                  MAIT="#3A6963",#MAIT
                                  Monocyte="#E95C39",# Mono细胞类型 
                                  NK="#E1A111", # NK细胞类型
                                  pDC="#938175", 
                                  Platelet="#BD956A",
                                  Tdn="#495608"# 其他细胞类型
                                                       )) {
    cold.ratio <- ratio.df %>% as.data.frame() %>% dplyr::filter(Condition %in% condition_seq )%>% .[,c('Sample', ct_level, 'ratio')] %>% 
                    pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>% column_to_rownames(var = 'Sample')
    cold.ratio
    cold.cor.pr = corr.test(cold.ratio, cold.ratio, method = 'pearson')
    cold.cor.pr
    
    # --- 02-plot HC CellProp cor network
    cor.df <- cold.cor.pr$r
    p.df <- cold.cor.pr$p
    n <- network(cor.df, directed = FALSE)
    n %v% ct_level <- colnames(p.df)
    
    #cor.gcmec
    e <- network.edgecount(n)
    cor.final <- NULL
    for(i in 1: (ncol(cor.df)-1)){
        cor.data <- cor.df[,i][-(1:i)]
        cor.final <- c(cor.final, cor.data)
    }
    #p.gcmec
    p.final <- NULL
    for(i in 1: (ncol(p.df)-1)){
        p.data <- p.df[,i][-(1:i)]
        p.final <- c(p.final, p.data)
    }
    set.edge.attribute(n, "cor", cor.final)
    set.edge.attribute(n, "p", p.final)
    set.edge.attribute(n, 'posneg', ifelse(cor.final >0, 'pos', 'neg'))
    all <- NULL
    for(i in 1:e){
        if(n$ mel[[i]]$ atl$ p < 0.05){
            all <- c(all,i)
        }
    }
    n$ mel <- n$ mel[all]
    
    set.seed(20)
    pdf(output_pdf,width = width , height  = height)
    print(ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(aes(linewidth = (abs(cor*cor)),color = posneg, alpha = abs(cor)), curvature = 0.2) +
      geom_nodes(aes(colour = !!sym(ct_level)),size = 5) +
      geom_nodelabel_repel(aes(label =  !!sym(ct_level)), size = 3, max.overlaps = 100, box.padding = unit(1, "lines"),colour = 'black') +
      scale_color_manual(values = color_seq) +
      #scale_size(range = c(12*min(abs(log2(wilcox_gc$com))), 12*max(abs(log2(wilcox_gc$com)))))+
      scale_linewidth(range = c(0.05,2))+
      scale_alpha(range = c(0.30046,0.75304))+
      theme_blank() +
       theme(
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = 'right',
        legend.spacing = unit(0.05, "lines"), # 调整两个 legend 之间的间距
        legend.key.size = unit(0.1, "inch"),
        plot.margin = unit(c(2, 2, 2, 2), "char")
      ) +
       guides(
        color = guide_legend(
          title = "Pos/Neg Correlation",
          override.aes = list(size = 1.5)  # 调小图例圆形标记大小
         )) +
      ggtitle(title_name)
          )
    dev.off()
}

# --- prepare the ratio data
classifier  <- "Condition"
cell_num.df  <- metadata.df %>%
  dplyr::group_by(Sample, !!sym(classifier)) %>% # !!sym(classifier) 这句可将字符串的双引号去除掉
  dplyr::summarise(total = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample,Condition,total)

head(cell_num.df ,n=2)
ratio.df  <- metadata.df %>%
  dplyr::group_by(Sample, celltypel2) %>%
  dplyr::summarise(sum = n()) %>%
  ungroup() %>%
  left_join(cell_num.df, by = "Sample") %>%
  rowwise() %>%
  dplyr::mutate(ratio = sum / total )
head(ratio.df)

color_seq = c("pos"= "red", 'neg' = "blue",
              Bcell="#AB3282",  # B细胞类型
              CD4T="#6778AE", # CD4 T细胞类型
              CD8T="#53A85F", # CD8 T细胞类型
              cDC="#9FA3A8",  # CDC细胞类型
              gdT="#23452F",  #gdT
              HPSC="#585658",#HPSC
              ILC="#F1BC62", #ILC
              MAIT="#3A6963",#MAIT
              Monocyte="#E95C39",# Mono细胞类型 
              NK="#E1A111", # NK细胞类型
              pDC="#938175", 
              Platelet="#BD956A",
              Tdn="#495608"# 其他细胞类型
                                   )

## 绘制A0时期的图
A0.ratio <- ratio.df %>% as.data.frame() %>% dplyr::filter(Condition %in% c('A0') )%>% .[,c('Sample', 'celltypel2', 'ratio')] %>% 
                pivot_wider(names_from = 'celltypel2', values_from = 'ratio', values_fill = 0) %>% column_to_rownames(var = 'Sample')
A0.ratio
A0.cor.pr = corr.test(A0.ratio, A0.ratio, method = 'pearson')
A0.cor.pr

# --- 02-plot HC CellProp cor network
cor.df <- A0.cor.pr$r
p.df <- A0.cor.pr$p
n <- network(cor.df, directed = FALSE)
n %v% "celltypel2" <- colnames(p.df)

#cor.gcmec
e <- network.edgecount(n)
cor.final <- NULL
for(i in 1: (ncol(cor.df)-1)){
    cor.data <- cor.df[,i][-(1:i)]
    cor.final <- c(cor.final, cor.data)
}
#p.gcmec
p.final <- NULL
for(i in 1: (ncol(p.df)-1)){
    p.data <- p.df[,i][-(1:i)]
    p.final <- c(p.final, p.data)
}
set.edge.attribute(n, "cor", cor.final)
set.edge.attribute(n, "p", p.final)
set.edge.attribute(n, 'posneg', ifelse(cor.final >0, 'pos', 'neg'))
all <- NULL
for(i in 1:e){
    if(n$ mel[[i]]$ atl$ p < 0.05){
        all <- c(all,i)
    }
}
n$ mel <- n$ mel[all]

set.seed(20)
pdf('figures/1_A0_CellProp_corr.pdf',width = 5.5, height  = 4)
ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linewidth = (abs(cor*cor)),color = posneg, alpha = abs(cor)), curvature = 0.2) +
  geom_nodes(aes(colour = celltypel2),size = 5) +
  geom_nodelabel_repel(aes(label =  celltypel2), size = 3, max.overlaps = 100, box.padding = unit(1, "lines"),colour = 'black') +
  scale_color_manual(values = c("pos"= "red", 'neg' = "blue",
                               Bcell="#AB3282",  # B细胞类型
                                  CD4T="#6778AE", # CD4 T细胞类型
                                  CD8T="#53A85F", # CD8 T细胞类型
                                  cDC="#9FA3A8",  # CDC细胞类型
                                  gdT="#23452F",  #gdT
                                  HPSC="#585658",#HPSC
                                  ILC="#F1BC62", #ILC
                                  MAIT="#3A6963",#MAIT
                                  Monocyte="#E95C39",# Mono细胞类型 
                                  NK="#E1A111", # NK细胞类型
                                  pDC="#938175", 
                                  Platelet="#BD956A",
                                  Tdn="#495608"# 其他细胞类型
                               )) +
  #scale_size(range = c(12*min(abs(log2(wilcox_gc$com))), 12*max(abs(log2(wilcox_gc$com)))))+
  scale_linewidth(range = c(0.05,2))+
  scale_alpha(range = c(0.30046,0.75304))+
  theme_blank() +
   theme(
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = 'right',
    legend.spacing = unit(0.05, "lines"), # 调整两个 legend 之间的间距
    legend.key.size = unit(0.1, "inch"),
    plot.margin = unit(c(2, 2, 2, 2), "char")
  ) +
   guides(
    color = guide_legend(
      title = "Pos/Neg Correlation",
      override.aes = list(size = 1.5)  # 调小图例圆形标记大小
     )) +
  ggtitle('A0: CellProp correlation')
dev.off()



plot_cellprop_corr(ratio.df, condition_seq = c('A0'), 
                   ct_level = 'celltypel2', 
                   output_pdf = 'figures/1_A0_CellProp_corr_2.pdf',
                  title_name='A0:CellProp correlation',
                  width = 5.5, height  = 4,color_seq=color_seq)




conditions <- unique(metadata.df$Condition)

# 循环每个 condition 调用 plot_cellprop_corr 函数
for (cond in conditions) {
  # 创建输出文件名，使用当前条件名称
  output_pdf <- paste0('figures/1_celltypel2_', cond, '_CellProp_corr.pdf')
  
  # 调用 plot_cellprop_corr 函数
  plot_cellprop_corr(
    ratio.df = ratio.df, 
    condition_seq = cond,  # 使用当前循环的 condition
    ct_level = 'celltypel2', 
    output_pdf = output_pdf,
    title_name = paste0(cond, ': CellProp correlation'),
    width = 5.5, 
    height = 4, 
    color_seq = color_seq
  )
}

## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A1','A2'), 
                   ct_level = 'celltypel2', 
                   output_pdf = 'figures/2_Celltyep2_A1A2_CellProp_corr.pdf',
                  title_name='A1 and A2: CellProp correlation',
                  width = 5.5, height  = 4,color_seq=color_seq)


## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('B1','B2'), 
                   ct_level = 'celltypel2', 
                   output_pdf = 'figures/2_Celltyep2_B1B2_CellProp_corr.pdf',
                  title_name='B1 and B2: CellProp correlation',
                  width = 5.5, height  = 4,color_seq=color_seq)
## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('C1','C2'), 
                   ct_level = 'celltypel2', 
                   output_pdf = 'figures/2_Celltyep2_C1C2_CellProp_corr.pdf',
                  title_name='C1 and C2: CellProp correlation',
                  width = 5.5, height  = 4,color_seq=color_seq)


## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A1','B1','C1'), 
                   ct_level = 'celltypel2', 
                   output_pdf = 'figures/2_Celltyep2_A1B1C1_CellProp_corr.pdf',
                  title_name='A1& B1& C1: CellProp correlation',
                  width = 5.5, height  = 4,color_seq=color_seq)

## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A2','B2','C2'), 
                   ct_level = 'celltypel2', 
                   output_pdf = 'figures/2_Celltyep2_A2B2C2_CellProp_corr.pdf',
                  title_name='A2& B2& C2: CellProp correlation',
                  width = 5.5, height  = 4,color_seq=color_seq)

## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A0','B0','C0','A1','B1','C1','A2','B2','C2'), 
                   ct_level = 'celltypel2', 
                   output_pdf = 'figures/2_Celltyep2_ALL_CellProp_corr.pdf',
                  title_name='ALL Condition: CellProp correlation',
                  width = 5.5, height  = 4,color_seq=color_seq)

metadata.df  <- read.csv('data/All_cell_metadata_df.csv')



# --- prepare the ratio data
classifier  <- "Condition"
cell_num.df  <- metadata.df %>%
  dplyr::group_by(Sample, !!sym(classifier)) %>% # !!sym(classifier) 这句可将字符串的双引号去除掉
  dplyr::summarise(total = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample,Condition,total)

head(cell_num.df ,n=2)
ratio.df  <- metadata.df %>%
  dplyr::group_by(Sample, celltypel3) %>%
  dplyr::summarise(sum = n()) %>%
  ungroup() %>%
  left_join(cell_num.df, by = "Sample") %>%
  rowwise() %>%
  dplyr::mutate(ratio = sum / total )
head(ratio.df,n=2)

color_seq = c("pos"= "red", 'neg' = "blue",
  'Bm_CD27+_IGHM-'="#E5D2DD",
  'Bm_CD27+_IGHM+'="#F3B1A0",
  'Bm_CD27+_IGHM+_SOX5'="#E59CC4",
  Bmn_IGHD="#AB3282",
  CD4_T_CCR6="#CCE0F5",
  CD4_Tcm="#C1E6F3",
  CD4_Tem="#91D0BE",
  CD4_Temra="#58A4C3",
  CD4_Th2="#57C3F3",
  CD4_Tn="#6778AE",
  CD4_Treg="#476D87",
  CD8_NELL2="#D6E7A9",
  CD8_Tem="#C5DEBA",
  CD8_Temra="#68A180",
  CD8_Tn="#53A85F",
  cDC1="#E0D4CA",
  cDC2="#9FA3A8",
  gdT_V1="#3A6963",
  gdT_V9="#23452F",
  HPSC="#585658",
  ILC1="#F1CC92",
  ILC2="#F1BC62",
  MAIT="#5F3D69",
  Mono_CD14="#E95C39",
  'Mono_CD14_ATG7+'="#A01319",
  Mono_CD14_CD16="#B53E2B",
  Mono_CD16="#712820",
  'Mono_CD16_ATG7+'="#A05401",
  NK_cyto_FCER1G="#E4C955",
  NK_cyto_KLRC2="#E1A111",
  NK_rest="#AA9A59",
  pDC="#968175",
  Plasma="#8C549C",
  Platelet="#BD956A",
  Tdn="#585658"
)

conditions <- unique(metadata.df$Condition)

# 循环每个 condition 调用 plot_cellprop_corr 函数
for (cond in conditions) {
  # 创建输出文件名，使用当前条件名称
  output_pdf <- paste0('figures/3_celltypel3_', cond, '_CellProp_corr.pdf')
  
  # 调用 plot_cellprop_corr 函数
  plot_cellprop_corr(
    ratio.df = ratio.df, 
    condition_seq = cond,  # 使用当前循环的 condition
    ct_level = 'celltypel3', 
    output_pdf = output_pdf,
    title_name = paste0(cond, ': CellProp correlation'),
    width = 10, 
    height = 6, 
    color_seq = color_seq
  )
}

## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A1','A2'), 
                   ct_level = 'celltypel3', 
                   output_pdf = 'figures/4_Celltyep3_A1A2_CellProp_corr.pdf',
                  title_name='A1 and A2: CellProp correlation',
                  width = 10, height  = 6,color_seq=color_seq)


## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('B1','B2'), 
                   ct_level = 'celltypel3', 
                   output_pdf = 'figures/4_Celltyep3_B1B2_CellProp_corr.pdf',
                  title_name='B1 and B2: CellProp correlation',
                  width = 10, height  = 6,color_seq=color_seq)
## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('C1','C2'), 
                   ct_level = 'celltypel3', 
                   output_pdf = 'figures/4_Celltyep3_C1C2_CellProp_corr.pdf',
                  title_name='C1 and C2: CellProp correlation',
                  width = 10, height  = 6,color_seq=color_seq)

## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A1','B1','C1'), 
                   ct_level = 'celltypel3', 
                   output_pdf = 'figures/4_Celltyep3_A1B1C1_CellProp_corr.pdf',
                  title_name='A1& B1& C1: CellProp correlation',
                  width = 10, height  = 6,color_seq=color_seq)

## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A2','B2','C2'), 
                   ct_level = 'celltypel3', 
                   output_pdf = 'figures/4_Celltyep3_A2B2C2_CellProp_corr.pdf',
                  title_name='A2& B2& C2: CellProp correlation',
                  width = 10, height  = 6,color_seq=color_seq)

c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2')

## 组合
plot_cellprop_corr(ratio.df, condition_seq = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
                   ct_level = 'celltypel3', 
                   output_pdf = 'figures/4_ALL_Celltyep3_CellProp_corr.pdf',
                  title_name='ALL: CellProp correlation',
                  width = 10, height  = 6,color_seq=color_seq)



metadata.df  <- read.csv('data/All_cell_metadata_df.csv')

plot_cellprop_corr_comparison <- function(ratio.df, condition_seq_all, condition_seq_a1, ct_level, 
                                          output_pdf, title_name, width = 5.5, height = 4,
                                          color_seq = c("pos"= "red", 'neg' = "blue",
                                                        Bcell="#AB3282",  
                                                        CD4T="#6778AE", 
                                                        CD8T="#53A85F", 
                                                        cDC="#9FA3A8",  
                                                        gdT="#23452F",  
                                                        HPSC="#585658",
                                                        ILC="#F1BC62", 
                                                        MAIT="#3A6963",
                                                        Monocyte="#E95C39",
                                                        NK="#E1A111", 
                                                        pDC="#938175", 
                                                        Platelet="#BD956A",
                                                        Tdn="#495608"),
                                          threshold = 0.3) {
    # --- Step 1: 计算全部时期相关性
    cold.ratio_all <- ratio.df %>%
        as.data.frame() %>%
        dplyr::filter(Condition %in% condition_seq_all) %>%
        .[, c('Sample', ct_level, 'ratio')] %>%
        pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>%
        column_to_rownames(var = 'Sample')

    cold.cor.pr <- corr.test(cold.ratio_all, cold.ratio_all, method = 'pearson')  # 获取相关性和显著性
    # --- Step 2: 构建网络
    cor.df <- cold.cor.pr$r
    p.df <- cold.cor.pr$p
    n <- network(cor.df, directed = FALSE)
    n %v% ct_level <- colnames(p.df)
    
    #cor.gcmec
    e <- network.edgecount(n)
    cor.final <- NULL
    for(i in 1: (ncol(cor.df)-1)){
        cor.data <- cor.df[,i][-(1:i)]
        cor.final <- c(cor.final, cor.data)
    }
    #p.gcmec
    p.final <- NULL
    for(i in 1: (ncol(p.df)-1)){
        p.data <- p.df[,i][-(1:i)]
        p.final <- c(p.final, p.data)
    }
    set.edge.attribute(n, "cor", cor.final)
    set.edge.attribute(n, "p", p.final)
    set.edge.attribute(n, 'posneg', ifelse(cor.final >0, 'pos', 'neg'))
    all <- NULL
    for(i in 1:e){
        if(n$ mel[[i]]$ atl$ p < 0.05){
            all <- c(all,i)
        }
    }
    n$ mel <- n$ mel[all]
    
    # --- Step 3: 计算A1时期相关性
    cold.ratio_a1 <- ratio.df %>%
        as.data.frame() %>%
        dplyr::filter(Condition %in% condition_seq_a1) %>%
        .[, c('Sample', ct_level, 'ratio')] %>%
        pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>%
        column_to_rownames(var = 'Sample')

    cor_a1 <- corr.test(cold.ratio_a1, cold.ratio_a1, method = 'pearson')$r

    # --- Step 4: 比较结果，计算节点属性（大小和颜色）
    cor_diff <- cor_a1 - cor.df  # 计算相关性差异
    node_diff <- rowMeans(cor_diff, na.rm = TRUE)  # 计算每个节点的平均变化

    # 设置节点大小和颜色
    node_size <- scales::rescale(abs(node_diff), to = c(3, 8))  # 节点大小范围
    node_color <- ifelse(node_diff > 0, "Down","UP")  # 正差异为红色，负差异为蓝色 

    # 添加节点属性
    network::set.vertex.attribute(n, "node_size", node_size)
    network::set.vertex.attribute(n, "node_color", node_color)

    

    # --- Step 5: 绘制网络图
    set.seed(20)
    pdf(output_pdf, width = width, height = height)
    print(ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
              # 边的属性保持不变
              geom_edges(aes(linewidth = (abs(cor * cor)), color = posneg, alpha = abs(cor)), curvature = 0.2) +
              # 节点属性根据比较结果变化
              geom_nodes(aes(size = node_size, color = node_color), alpha = 0.8) +
              geom_nodelabel_repel(aes(label = !!sym(ct_level)), size = 3, max.overlaps = 100, box.padding = unit(1, "lines"), colour = 'black') +
              scale_color_manual(values = c( "#68A180","blue", "red","#573333FF")) +  
              scale_size(range = c(3, 8)) +
              scale_linewidth(range = c(0.05, 2)) +
              scale_alpha(range = c(0.5, 0.75)) +
              theme_blank() +
              theme(
                  legend.title = element_text(size = 6),
                  legend.text = element_text(size = 6),
                  legend.position = 'right',
                  legend.spacing = unit(0.05, "lines"),
                  legend.key.size = unit(0.1, "inch"),
                  plot.margin = unit(c(2, 2, 2, 2), "char")
              ) +
              guides(
                  color = guide_legend(
                      title = "Node Changes",
                      override.aes = list(size = 1.5)
                  )
              ) +
              ggtitle(title_name))
    dev.off()
}

# --- prepare the ratio data
classifier  <- "Condition"
cell_num.df  <- metadata.df %>%
  dplyr::group_by(Sample, !!sym(classifier)) %>% # !!sym(classifier) 这句可将字符串的双引号去除掉
  dplyr::summarise(total = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample,Condition,total)

head(cell_num.df ,n=2)
ratio.df  <- metadata.df %>%
  dplyr::group_by(Sample, celltypel2) %>%
  dplyr::summarise(sum = n()) %>%
  ungroup() %>%
  left_join(cell_num.df, by = "Sample") %>%
  rowwise() %>%
  dplyr::mutate(ratio = sum / total )
head(ratio.df)

# 定义参数列表
comparison_params <- list(
    list(condition_seq_a1 = c('A1'), output_pdf = 'figures/5_Comparison_ALL_vs_A1_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A1 - CellProp Correlation'),
    list(condition_seq_a1 = c('B1'), output_pdf = 'figures/5_Comparison_ALL_vs_B1_CellProp_corr.pdf', title_name = 'Comparison: ALL vs B1 - CellProp Correlation'),
    list(condition_seq_a1 = c('A1', 'B1', 'C1'), output_pdf = 'figures/5_Comparison_ALL_vs_A1B1C1_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A1B1C1 - CellProp Correlation'),
    list(condition_seq_a1 = c('A1', 'A2'), output_pdf = 'figures/5_Comparison_ALL_vs_A1A2_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A1A2 - CellProp Correlation'),
    list(condition_seq_a1 = c('A2', 'B2', 'C2'), output_pdf = 'figures/5_Comparison_ALL_vs_A2B2C2_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A2B2C2 - CellProp Correlation'),
    list(condition_seq_a1 = c('B1', 'B2'), output_pdf = 'figures/5_Comparison_ALL_vs_B1B2_CellProp_corr.pdf', title_name = 'Comparison: ALL vs B1B2 - CellProp Correlation')
)

# 循环调用函数
for (params in comparison_params) {
    plot_cellprop_corr_comparison(
        ratio.df, 
        condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
        condition_seq_a1 = params$condition_seq_a1, 
        ct_level = 'celltypel2', 
        output_pdf = params$output_pdf,
        title_name = params$title_name,
        width = 5.5, 
        height = 4.5, 
        color_seq = color_seq, 
        threshold = 0
    )
}

plot_cellprop_corr_comparison(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('A1', 'B1', 'C1'), 
    ct_level = 'celltypel2', 
    output_pdf = 'figures/5_Comparison_ALL_vs_A1B1C1_CellProp_corr.pdf',
    title_name = 'Comparison: ALL vs A1B1C1 - CellProp Correlation',
    width = 5.5, height = 4, color_seq = color_seq,threshold = 0.4
)

# --- prepare the ratio data
classifier  <- "Condition"
cell_num.df  <- metadata.df %>%
  dplyr::group_by(Sample, !!sym(classifier)) %>% # !!sym(classifier) 这句可将字符串的双引号去除掉
  dplyr::summarise(total = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sample,Condition,total)

head(cell_num.df ,n=2)
ratio.df  <- metadata.df %>%
  dplyr::group_by(Sample, celltypel3) %>%
  dplyr::summarise(sum = n()) %>%
  ungroup() %>%
  left_join(cell_num.df, by = "Sample") %>%
  rowwise() %>%
  dplyr::mutate(ratio = sum / total )
head(ratio.df,n=2)

# 定义参数列表
comparison_params <- list(
    list(condition_seq_a1 = c('A1'), output_pdf = 'figures/6_L3_Comparison_ALL_vs_A1_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A1 - CellProp Correlation'),
    list(condition_seq_a1 = c('B1'), output_pdf = 'figures/6_L3_Comparison_ALL_vs_B1_CellProp_corr.pdf', title_name = 'Comparison: ALL vs B1 - CellProp Correlation'),
    list(condition_seq_a1 = c('A1', 'B1', 'C1'), output_pdf = 'figures/6_L3_Comparison_ALL_vs_A1B1C1_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A1B1C1 - CellProp Correlation'),
    list(condition_seq_a1 = c('A1', 'A2'), output_pdf = 'figures/6_L3_Comparison_ALL_vs_A1A2_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A1A2 - CellProp Correlation'),
    list(condition_seq_a1 = c('A2', 'B2', 'C2'), output_pdf = 'figures/6_L3_Comparison_ALL_vs_A2B2C2_CellProp_corr.pdf', title_name = 'Comparison: ALL vs A2B2C2 - CellProp Correlation'),
    list(condition_seq_a1 = c('B1', 'B2'), output_pdf = 'figures/6_L3_Comparison_ALL_vs_B1B2_CellProp_corr.pdf', title_name = 'Comparison: ALL vs B1B2 - CellProp Correlation')
)

# 循环调用函数
for (params in comparison_params) {
    plot_cellprop_corr_comparison(
        ratio.df, 
        condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
        condition_seq_a1 = params$condition_seq_a1, 
        ct_level = 'celltypel3', 
        output_pdf = params$output_pdf,
        title_name = params$title_name,
        width = 10, 
        height = 7, 
        color_seq = color_seq, 
        threshold = 0
    )
}

#把连接改成直线展示

plot_cellprop_corr_comparison <- function(ratio.df, condition_seq_all, condition_seq_a1, ct_level, 
                                          output_pdf, title_name, width = 5.5, height = 4,
                                          color_seq = c("pos"= "red", 'neg' = "blue",
                                                        Bcell="#AB3282",  
                                                        CD4T="#6778AE", 
                                                        CD8T="#53A85F", 
                                                        cDC="#9FA3A8",  
                                                        gdT="#23452F",  
                                                        HPSC="#585658",
                                                        ILC="#F1BC62", 
                                                        MAIT="#3A6963",
                                                        Monocyte="#E95C39",
                                                        NK="#E1A111", 
                                                        pDC="#938175", 
                                                        Platelet="#BD956A",
                                                        Tdn="#495608"),
                                          threshold = 0.3) {
    # --- Step 1: 计算全部时期相关性
    cold.ratio_all <- ratio.df %>%
        as.data.frame() %>%
        dplyr::filter(Condition %in% condition_seq_all) %>%
        .[, c('Sample', ct_level, 'ratio')] %>%
        pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>%
        column_to_rownames(var = 'Sample')

    cold.cor.pr <- corr.test(cold.ratio_all, cold.ratio_all, method = 'pearson')  # 获取相关性和显著性
    # --- Step 2: 构建网络
    cor.df <- cold.cor.pr$r
    p.df <- cold.cor.pr$p
    n <- network(cor.df, directed = FALSE)
    n %v% ct_level <- colnames(p.df)
    
    #cor.gcmec
    e <- network.edgecount(n)
    cor.final <- NULL
    for(i in 1: (ncol(cor.df)-1)){
        cor.data <- cor.df[,i][-(1:i)]
        cor.final <- c(cor.final, cor.data)
    }
    #p.gcmec
    p.final <- NULL
    for(i in 1: (ncol(p.df)-1)){
        p.data <- p.df[,i][-(1:i)]
        p.final <- c(p.final, p.data)
    }
    set.edge.attribute(n, "cor", cor.final)
    set.edge.attribute(n, "p", p.final)
    set.edge.attribute(n, 'posneg', ifelse(cor.final >0, 'pos', 'neg'))
    all <- NULL
    for(i in 1:e){
        if(n$ mel[[i]]$ atl$ p < 0.05){
            all <- c(all,i)
        }
    }
    n$ mel <- n$ mel[all]
    # --- Step 3: 计算A1时期相关性
    cold.ratio_a1 <- ratio.df %>%
        as.data.frame() %>%
        dplyr::filter(Condition %in% condition_seq_a1) %>%
        .[, c('Sample', ct_level, 'ratio')] %>%
        pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>%
        column_to_rownames(var = 'Sample')

    cor_a1 <- corr.test(cold.ratio_a1, cold.ratio_a1, method = 'pearson')$r

    # --- Step 4: 比较结果，计算节点属性（大小和颜色）
    cor_diff <- cor_a1 - cor.df  # 计算相关性差异
    node_diff <- rowMeans(cor_diff, na.rm = TRUE)  # 计算每个节点的平均变化

    # 设置节点大小和颜色
    node_size <- scales::rescale(abs(node_diff), to = c(3, 8))  # 节点大小范围
    node_color <- ifelse(node_diff > 0, "Down","UP")  # 正差异为红色，负差异为蓝色 

    # 添加节点属性
    network::set.vertex.attribute(n, "node_size", node_size)
    network::set.vertex.attribute(n, "node_color", node_color)

    # --- Step 5: 绘制网络图
    set.seed(20)
    pdf(output_pdf, width = width, height = height)
    print(ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
              # 边的属性保持不变
              geom_edges(aes(linewidth = (abs(cor * cor)), color = posneg, alpha = abs(cor)), curvature = 0) +
              # 节点属性根据比较结果变化，将 alpha 设置为 1（完全不透明）
              geom_nodes(aes(size = node_size, color = node_color), alpha = 1) +
              geom_nodelabel_repel(aes(label = !!sym(ct_level)), size = 3, max.overlaps = 100, box.padding = unit(1, "lines"), colour = 'black') +
              scale_color_manual(values = c( "#68A180","blue", "red","#573333FF")) +  
              scale_size(range = c(3, 8)) +
              scale_linewidth(range = c(0.05, 2)) +
              scale_alpha(range = c(0.5, 0.75)) +  # 边的透明度范围保持不变
              theme_blank() +
              theme(
                  legend.title = element_text(size = 6),
                  legend.text = element_text(size = 6),
                  legend.position = 'right',
                  legend.spacing = unit(0.05, "lines"),
                  legend.key.size = unit(0.1, "inch"),
                  plot.margin = unit(c(2, 2, 2, 2), "char")
              ) +
              guides(
                  color = guide_legend(
                      title = "Node Changes",
                      override.aes = list(size = 1.5)
                  )
              ) +
              ggtitle(title_name))
    dev.off()
}



plot_cellprop_corr_comparison(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('A1', 'B1'), 
    ct_level = 'celltypel3', 
    output_pdf = 'figures/7_L3_Comparison_ALL_vs_A1B1_CellProp_corr.pdf',
    title_name = 'Comparison: A1B1  vs ALL- CellProp Correlation',
    width = 10, height = 7, color_seq = color_seq,threshold = 0.4
)


plot_cellprop_corr_comparison(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('B2', 'C2'), 
    ct_level = 'celltypel3', 
    output_pdf = 'figures/7_L3_Comparison_ALL_vs_B2C2_CellProp_corr.pdf',
    title_name = 'Comparison: B2 C2  vs ALL- CellProp Correlation',
    width = 10, height = 7, color_seq = color_seq,threshold = 0.4
)

#把连接改成直线展示
plot_cellprop_corr_comparison_2 <- function(ratio.df, condition_seq_all, condition_seq_a1,
                                           condition_seq_a2,ct_level, 
                                          output_pdf, title_name, width = 5.5, height = 4,
                                          color_seq = c("pos"= "red", 'neg' = "blue",
                                                        Bcell="#AB3282",  
                                                        CD4T="#6778AE", 
                                                        CD8T="#53A85F", 
                                                        cDC="#9FA3A8",  
                                                        gdT="#23452F",  
                                                        HPSC="#585658",
                                                        ILC="#F1BC62", 
                                                        MAIT="#3A6963",
                                                        Monocyte="#E95C39",
                                                        NK="#E1A111", 
                                                        pDC="#938175", 
                                                        Platelet="#BD956A",
                                                        Tdn="#495608"),
                                          threshold = 0.3) {
    # --- Step 1: 计算全部时期相关性
    cold.ratio_all <- ratio.df %>%
        as.data.frame() %>%
        dplyr::filter(Condition %in% condition_seq_all) %>%
        .[, c('Sample', ct_level, 'ratio')] %>%
        pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>%
        column_to_rownames(var = 'Sample')

    cold.cor.pr <- corr.test(cold.ratio_all, cold.ratio_all, method = 'pearson')  # 获取相关性和显著性
    

    # --- Step 2: 构建网络
    cor.df <- cold.cor.pr$r
    p.df <- cold.cor.pr$p
    
    n <- network(cor.df, directed = FALSE)
    n %v% ct_level <- colnames(p.df)
    
    #cor.gcmec
    e <- network.edgecount(n)
    cor.final <- NULL
    for(i in 1: (ncol(cor.df)-1)){
        cor.data <- cor.df[,i][-(1:i)]
        cor.final <- c(cor.final, cor.data)
    }
    #p.gcmec
    p.final <- NULL
    for(i in 1: (ncol(p.df)-1)){
        p.data <- p.df[,i][-(1:i)]
        p.final <- c(p.final, p.data)
    }
    set.edge.attribute(n, "cor", cor.final)
    set.edge.attribute(n, "p", p.final)
    set.edge.attribute(n, 'posneg', ifelse(cor.final >0, 'pos', 'neg'))

    # Filter edges by p-value and threshold
      all <- NULL
      for (i in 1:e) {
        if (n$mel[[i]]$atl$p < 0.05 && abs(n$mel[[i]]$atl$cor) >= threshold) {
          all <- c(all, i)
        }
      }
      n$mel <- n$mel[all]

    # --- Step 3: 计算两个时期相关性
    cold.ratio_a1 <- ratio.df %>%
        as.data.frame() %>%
        dplyr::filter(Condition %in% condition_seq_a1) %>%
        .[, c('Sample', ct_level, 'ratio')] %>%
        pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>%
        column_to_rownames(var = 'Sample')

    cor_a1 <- corr.test(cold.ratio_a1, cold.ratio_a1, method = 'pearson')$r
    
    cold.ratio_a2 <- ratio.df %>%
        as.data.frame() %>%
        dplyr::filter(Condition %in% condition_seq_a2) %>%
        .[, c('Sample', ct_level, 'ratio')] %>%
        pivot_wider(names_from = ct_level, values_from = 'ratio', values_fill = 0) %>%
        column_to_rownames(var = 'Sample')

    cor_a2 <- corr.test(cold.ratio_a2, cold.ratio_a2, method = 'pearson')$r

    # --- Step 4: 比较结果，计算节点属性（大小和颜色）
    cor_diff <- cor_a1 - cor_a2  # 计算相关性差异
    node_diff <- rowMeans(cor_diff, na.rm = TRUE)  # 计算每个节点的平均变化

    # 设置节点大小和颜色
    node_size <- scales::rescale(abs(node_diff), to = c(3, 8))  # 节点大小范围
    node_color <- ifelse(node_diff > 0, "Down","UP")  # 正差异为红色，负差异为蓝色 

    # 添加节点属性
    network::set.vertex.attribute(n, "node_size", node_size)
    network::set.vertex.attribute(n, "node_color", node_color)

    # --- Step 5: 绘制网络图
    set.seed(20)
    pdf(output_pdf, width = width, height = height)
    print(ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
              # 边的属性保持不变
              geom_edges(aes(linewidth = (abs(cor * cor)), color = posneg, alpha = abs(cor)), curvature = 0) +
              # 节点属性根据比较结果变化，将 alpha 设置为 1（完全不透明）
              geom_nodes(aes(size = node_size, color = node_color), alpha = 1) +
              geom_nodelabel_repel(aes(label = !!sym(ct_level)), size = 3, max.overlaps = 100, box.padding = unit(1, "lines"), colour = 'black') +
              scale_color_manual(values = c( "#68A180","blue", "red","#573333FF")) +  
              scale_size(range = c(3, 8)) +
              scale_linewidth(range = c(0.05, 2)) +
              scale_alpha(range = c(0.5, 0.75)) +  # 边的透明度范围保持不变
              theme_blank() +
              theme(
                  legend.title = element_text(size = 6),
                  legend.text = element_text(size = 6),
                  legend.position = 'right',
                  legend.spacing = unit(0.05, "lines"),
                  legend.key.size = unit(0.1, "inch"),
                  plot.margin = unit(c(2, 2, 2, 2), "char")
              ) +
              guides(
                  color = guide_legend(
                      title = "Node Changes",
                      override.aes = list(size = 1.5)
                  )
              ) +
              ggtitle(title_name))
    dev.off()
}

plot_cellprop_corr_comparison_2(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('A1', 'B1'), 
    condition_seq_a2 = c('A0'), 
    ct_level = 'celltypel3', 
    output_pdf = 'figures/8_L3_Comparison_A0_vs_A1B1_CellProp_corr.pdf',
    title_name = 'Comparison: A1B1  vs A0 - CellProp Correlation',
    width = 10, height = 7, color_seq = color_seq,threshold =  0.366
)

plot_cellprop_corr_comparison_2(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('A1'), 
    condition_seq_a2 = c('A0'), 
    ct_level = 'celltypel3', 
    output_pdf = 'figures/8_L3_Comparison_A0_vs_A1_CellProp_corr.pdf',
    title_name = 'Comparison: A1  vs A0 - CellProp Correlation',
    width = 10, height = 7, color_seq = color_seq,threshold =  0.366
)

plot_cellprop_corr_comparison_2(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('A1', 'B1'), 
     condition_seq_a2 = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    ct_level = 'celltypel3', 
    output_pdf = 'figures/8_L3_Comparison_ALL_vs_A1B1_CellProp_corr.pdf',
    title_name = 'Comparison: A1B1  vs ALL- CellProp Correlation',
    width = 10, height = 7, color_seq = color_seq,threshold =  0.366
)


plot_cellprop_corr_comparison_2(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('A1', 'B1'), 
     condition_seq_a2 = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    ct_level = 'celltypel3', 
    output_pdf = 'figures/8_L3_Comparison_ALL_vs_A1B1_CellProp_corr_2.pdf',
    title_name = 'Comparison: A1B1  vs ALL- CellProp Correlation',
    width = 9, height = 6, color_seq = color_seq,threshold =  0.366
)


plot_cellprop_corr_comparison_2(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('B2', 'C2'), 
    condition_seq_a2 = c('A0'), 
    ct_level = 'celltypel3', 
    output_pdf = 'figures/8_L3_Comparison_A0_vs_B2C2_CellProp_corr.pdf',
    title_name = 'Comparison: B2 C2  vs A0- CellProp Correlation',
    width = 10, height = 7, color_seq = color_seq,threshold =  0.366
)

plot_cellprop_corr_comparison_2(
    ratio.df, 
    condition_seq_all = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'), 
    condition_seq_a1 = c('B2', 'C2'), 
    condition_seq_a2 = c('A0', 'B0', 'C0', 'A1', 'B1', 'C1', 'A2', 'B2', 'C2'),
    ct_level = 'celltypel3', 
    output_pdf = 'figures/8_L3_Comparison_ALL_vs_B2C2_CellProp_corr.pdf',
    title_name = 'Comparison: B2 C2  vs ALL- CellProp Correlation',
    width = 10, height = 7, color_seq = color_seq,threshold = 0.366
)








