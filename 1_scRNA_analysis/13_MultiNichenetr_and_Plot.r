library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)

organism = "human"

options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = 
    readRDS(
      "../files/lr_network_human_allInfo_30112033.rds"
      ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(
    "../files/ligand_target_matrix_nsga2r_final.rds"
    )
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(
    "../files/lr_network_mouse_allInfo_30112033.rds"
    ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(
    "../files/ligand_target_matrix_nsga2r_final_mouse.rds"
    )
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}

seurat_obj<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/5_new_level_1_2/5_merge_data_new_level_1.rds")

table(seurat_obj$ct_level_4,seurat_obj$Sample)

seurat_obj$ct_level_4[seurat_obj$ct_level_4 == "Bm_CD27+_IGHM-"] <- "Bm_CD27_IGHM_D"
seurat_obj$ct_level_4[seurat_obj$ct_level_4 == "Bm_CD27+_IGHM+"] <- "Bm_CD27_IGHM_H"
seurat_obj$ct_level_4[seurat_obj$ct_level_4 == "Bm_CD27+_IGHM+_SOX5"] <- "Bm_CD27_IGHM_SOX5"
seurat_obj$ct_level_4[seurat_obj$ct_level_4 == "Mono_CD14_ATG7+"] <- "Mono_CD14_ATG7"
seurat_obj$ct_level_4[seurat_obj$ct_level_4 == "Mono_CD16_ATG7+"] <- "Mono_CD16_ATG7"

table(seurat_obj$ct_level_4,seurat_obj$Sample)

sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")

sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

sample_id = "Sample"
group_id = "Condition"
celltype_id = "ct_level_4"

#无纠正的批处理效应或协变量
covariates = NA
batches = NA

#定义感兴趣的发送方和接收方单元格类型
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
          ]
# #只保留兴趣条件
# conditions_keep = c("M", "S", "A")
# sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% 
#             conditions_keep
#           ]
senders_oi
receivers_oi 

sce$ct_level_4 <- make.names(sce$ct_level_4)

#细胞类型过滤：确定哪些细胞类型充分存在
min_cells = 3
abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
  )

#细胞类型丰度信息的解释
pdf('4_abundance_info_abund_plot_sample.pdf',30,20)
abundance_info$abund_plot_sample
dev.off()

#基于细胞类型丰度信息的细胞类型过滤
sample_group_celltype_df = abundance_info$abundance_data %>% 
  filter(n > min_cells) %>% 
  ungroup() %>% 
  distinct(sample_id, group_id) %>% 
  cross_join(
    abundance_info$abundance_data %>% 
      ungroup() %>% 
      distinct(celltype_id)
    ) %>% 
  arrange(sample_id)

abundance_df = sample_group_celltype_df %>% left_join(
  abundance_info$abundance_data %>% ungroup()
  )

abundance_df$n[is.na(abundance_df$n)] = 0
abundance_df$keep[is.na(abundance_df$keep)] = FALSE
abundance_df_summarized = abundance_df %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition

celltypes_present_one_condition = abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

absent_celltypes = abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  count() %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)
  
print("condition-specific celltypes:")
## [1] "condition-specific celltypes:"
print(condition_specific_celltypes)
## character(0)
  
print("absent celltypes:")
## [1] "absent celltypes:"
print(absent_celltypes)
## character(0)

# 过滤掉特定条件的单元类型("HPSC" "cDC1")
analyse_condition_specific_celltypes = FALSE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
          ]

senders_oi
receivers_oi

#设置min_sample_prop = 0.50，这意味着如果样本数最低的组有4个样本，则基因应至少在2个样本中表达。
min_sample_prop = 0.30
#fraction_cutoff = 0.05，这意味着基因应在样本中至少5%的细胞中显示非零表达值。
fraction_cutoff = 0.005

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

#只保留至少一种细胞类型表达的基因：
genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

#每个基因/细胞类型/样本的标准化伪散量表达值可以通过以下方式检查：
abundance_expression_info$celltype_info$pb_df %>% head()

#每个条件/组的这些样本级表达式值的平均值可以通过以下方式检查：
abundance_expression_info$celltype_info$pb_df_group %>% head()

#检查配体-受体相互作用的这些值可以通过以下方式完成：
abundance_expression_info$sender_receiver_info$pb_df %>% head()

abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

# contrasts_oi = c("'A0','A1','A2','B0','B1','B2','C0','C1','C2'")
# contrast_tbl = tibble(contrast = 
#                        c("A0","A1","A2","B0","B1","B2","C0","C1","C2")
#                       group = c("A0","A1","A2","B0","B1","B2","C0","C1","C2"))

contrasts_oi = c("'A0-(A1+A2+B0+B1+B2+C0+C1+C2)/8','A1-(A0+A2+B0+B1+B2+C0+C1+C2)/8','A2-(A1+A0+B0+B1+B2+C0+C1+C2)/8','B0-(A1+A2+A0+B1+B2+C0+C1+C2)/8','B1-(A1+A2+B0+A0+B2+C0+C1+C2)/8','B2-(A1+A2+B0+B1+A0+C0+C1+C2)/8','C0-(A1+A2+B0+B1+B2+A0+C1+C2)/8','C1-(A1+A2+B0+B1+B2+C0+A0+C2)/8','C2-(A1+A2+B0+B1+B2+C0+C1+A0)/8'")

contrast_tbl = tibble(contrast = c('A0-(A1+A2+B0+B1+B2+C0+C1+C2)/8','A1-(A0+A2+B0+B1+B2+C0+C1+C2)/8','A2-(A1+A0+B0+B1+B2+C0+C1+C2)/8','B0-(A1+A2+A0+B1+B2+C0+C1+C2)/8','B1-(A1+A2+B0+A0+B2+C0+C1+C2)/8','B2-(A1+A2+B0+B1+A0+C0+C1+C2)/8','C0-(A1+A2+B0+B1+B2+A0+C1+C2)/8','C1-(A1+A2+B0+B1+B2+C0+A0+C2)/8','C2-(A1+A2+B0+B1+B2+C0+C1+A0)/8'), 
                      group = c("A0","A1","A2","B0","B1","B2","C0","C1","C2"))

head(frq_list$expressed_df)

min_cells = 3
DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)

#检查表格中的DE输出信息，其中包含每个基因-细胞类型-对比度的logFC和p值
DE_info$celltype_de$de_output_tidy %>% head()

pdf('7_DE_info_hist_pvals_2.pdf',40,35)
DE_info$hist_pvals
dev.off()

#估计经验p值,FALSE不计算
empirical_pval = FALSE

if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

#检查默认支架的geneset_oi-vs-background比率
logFC_threshold = 0.30
p_val_threshold = 0.05
p_val_adj = FALSE 

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
    ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval

#执行配体活动分析和配体目标推断
top_n_target = 250

verbose = TRUE
cores_system = 8
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

n.cores

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

ligand_activity_down = FALSE

sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    scenario = "regular", # all prioritization criteria will be weighted equally
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender,
    ligand_activity_down = ligand_activity_down
  ))

prioritization_tables$group_prioritization_tbl %>% head(20)

lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, #
  prioritization_tables = prioritization_tables, #
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
  )

head(p_val_threshold)





path = "./"
multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
multinichenet_output = make_lite_output(multinichenet_output)

save = FALSE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))

}

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))

}

head(prioritization_tables,n=1)

head(grouping_tbl)

head(lr_target_prior_cor)




# Plot Show ---------------------------------------------------------------

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
#devtools::install("/share/home/qlab/qlab_yyp/backup/packages/multinichenetr-main")

organism = "human"

options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = 
    readRDS(
      "../files/lr_network_human_allInfo_30112033.rds"
    ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(
    "../files/ligand_target_matrix_nsga2r_final.rds"
  )
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(
    "files/lr_network_mouse_allInfo_30112033.rds"
  ) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(
    "files/ligand_target_matrix_nsga2r_final_mouse.rds"
  )
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
}

multinichenet_output <- readRDS("multinichenet_output.rds")

source("new_plot_line.R") #make_sample_lr_prod_activity_plots_Omnipath_line

1

unique(multinichenet_output$celltype_info$avg_df$celltype)

# 创建并保存图像的函数
create_and_save_plot <- function(prioritized_tbl_oi_all, group_oi, top_n, height) {
  # 根据给定的组筛选数据
  prioritized_tbl_oi = prioritized_tbl_oi_all %>% 
    filter(group == group_oi)
  
  # 生成图像
  plot_oi = make_sample_lr_prod_activity_plots_Omnipath_line(
    filename=paste0("figures/1.0_group_top_10_20_50/data/1_", group_oi, "_top_", top_n),
    multinichenet_output$prioritization_tables, 
    prioritized_tbl_oi %>% inner_join(lr_network_all)
  )
  
  # 创建一个PDF文件并保存图像
  pdf(paste0("figures/1.0_group_top_10_20_50/1_top_", top_n, "_per_group_", group_oi, "_.pdf"), 35, height)
  print(plot_oi)
  dev.off()
}

# 循环展示所有组的图形的函数
create_and_save_plots_for_groups <- function(top_n, height) {
  # 获取优先级配对
  prioritized_tbl_oi_all = get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    top_n = top_n, 
    rank_per_group = TRUE #是否需要局部排序
  )
  
  # 定义组
  groups = c("A0","A1","A2","B0","B1","B2","C0","C1","C2")
  
  # 循环展示所有组的图形
  sapply(groups, function(group_oi) create_and_save_plot(prioritized_tbl_oi_all, group_oi, top_n, height))
}

# 使用函数
create_and_save_plots_for_groups(10, 5)
create_and_save_plots_for_groups(20, 8)
# create_and_save_plots_for_groups(30, 12)
# create_and_save_plots_for_groups(50, 18)

#对所有细胞类型进行展示
create_plots_for_celltype <- function(celltype, top_n_values, height_values) {
  for (i in seq_along(top_n_values)) {
    prioritized_tbl_oi_all = get_top_n_lr_pairs(
      multinichenet_output$prioritization_tables, 
      top_n = top_n_values[i], 
      senders_oi = celltype,
      #receivers_oi = celltype,
      rank_per_group = FALSE
    )
    
    plot_oi = make_sample_lr_prod_activity_plots_Omnipath_line(
      
      filename=paste0("figures/1.1_celltype_senders_top_10_20_50/data/1.1_", celltype, "_top_", top_n_values[i]),
      multinichenet_output$prioritization_tables, 
      prioritized_tbl_oi_all %>% inner_join(lr_network_all)
    )
    
    pdf(paste0("figures/1.1_celltype_senders_top_10_20_50/1.1_", celltype, "_top_", top_n_values[i], ".pdf"), 35, height_values[i])
    print(plot_oi)
    dev.off()
  }
}

# 对所有细胞类型进行处理
celltypes <- unique(multinichenet_output$celltype_info$avg_df$celltype)
top_n_values <- c(10, 20, 50)
height_values <- c(5, 8, 18)

for (celltype in celltypes) {
  create_plots_for_celltype(celltype, top_n_values, height_values)
}

#对所有细胞类型进行展示
create_plots_for_celltype <- function(celltype, top_n_values, height_values) {
  for (i in seq_along(top_n_values)) {
    prioritized_tbl_oi_all = get_top_n_lr_pairs(
      multinichenet_output$prioritization_tables, 
      top_n = top_n_values[i], 
      #senders_oi = celltype,
      receivers_oi = celltype,
      rank_per_group = FALSE
    )
    
    plot_oi = make_sample_lr_prod_activity_plots_Omnipath_line(
      
      filename=paste0("figures/1.2_celltype_receivers_top_10_20_50/data/1.2_", celltype, "_top_", top_n_values[i]),
      multinichenet_output$prioritization_tables, 
      prioritized_tbl_oi_all %>% inner_join(lr_network_all)
    )
    
    pdf(paste0("figures/1.2_celltype_receivers_top_10_20_50/1.2_", celltype, "_top_", top_n_values[i], ".pdf"), 35, height_values[i])
    print(plot_oi)
    dev.off()
  }
}

# 对所有细胞类型进行处理
celltypes <- unique(multinichenet_output$celltype_info$avg_df$celltype)
top_n_values <- c(10, 20, 50)
height_values <- c(5, 8, 18)

for (celltype in celltypes) {
  create_plots_for_celltype(celltype, top_n_values, height_values)
}

library(ggplot2)
library(RColorBrewer)

plot_func <- function(top_n,widths) {    
  conditions <- c("A0", "A1", "A2","B0","B1","B2","C0","C1","C2")
  result <- data.frame()
  for (condition in conditions) {
    condition_celltype <- read.csv(paste0("figures/1.0_group_top_10_20_50/data/1_", condition, "_top_",top_n,"_omnipath_data_p_omnipath.csv"))
    # 取前10行
    condition_data <- condition_celltype[1:top_n,]
    condition_data$Condition <- condition
    result <- rbind(result, condition_data)
  }
  result$sender = sapply(strsplit(result$sender_receiver, " --> "), "[", 1)
  CD4T <- c("CD4_Tn","CD4_Tcm","CD4_Tem","CD4_Treg","CD4_T_CCR6","CD4_Temra","CD4_Th2")
  CD8T <- c("CD8_Tn","CD8_Tem","CD8_Temra","CD8_NELL2")
  NK <- c("NK1","NK2","NK3")
  Monocyte <- c("Mono_CD14","Mono_CD16","Mono_CD14_CD16","Mono_CD14_ATG7","Mono_CD16_ATG7")
  Bcell <- c("Bmn_IGHD","Bm_CD27_IGHM_D","Bm_CD27_IGHM_H","Bm_CD27_IGHM_SOX5","Plasma")
  
  result$receiver_level_1 <- result$receiver
  result$receiver_level_1[result$receiver %in% CD4T] <- "CD4T"
  result$receiver_level_1[result$receiver %in% CD8T] <- "CD8T"
  result$receiver_level_1[result$receiver %in% NK] <- "NK"
  result$receiver_level_1[result$receiver %in% Bcell] <- "Bcell"
  result$receiver_level_1[result$receiver %in% Monocyte] <- "Monocyte"
  
  result$sender_level_1 <- result$sender
  result$sender_level_1[result$sender %in% CD4T] <- "CD4T"
  result$sender_level_1[result$sender %in% CD8T] <- "CD8T"
  result$sender_level_1[result$sender %in% NK] <- "NK"
  result$sender_level_1[result$sender %in% Bcell] <- "Bcell"
  result$sender_level_1[result$sender %in% Monocyte] <- "Monocyte"
  
  
  
  colpalette <- colorRampPalette(brewer.pal(11, "Spectral"))(36) 
  p1 <- ggplot(result, aes(x = Condition, fill = receiver)) + 
    geom_bar(position = "fill", width = 0.8) +  ## 通过调整 width 参数，控制条形图的宽度
    scale_fill_manual(values = colpalette) +
    theme_minimal(base_size = 16) +   ## 全局基础字体大小设定为16
    ylab("Number of receiver Cell Types")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.15, 'lines')  
    ) 
  
  p2 <- ggplot(result, aes(x = Condition, fill = sender)) + 
    geom_bar(position = "fill", width = 0.8) +  ## 通过调整 width 参数，控制条形图的宽度
    scale_fill_manual(values = colpalette) +
    theme_minimal(base_size = 16) +   ## 全局基础字体大小设定为16
    ylab("Number of sender Cell Types")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.15, 'lines')  
    ) 
  
  colpalette <- colorRampPalette(brewer.pal(11, "Spectral"))(16) 
  p3 <- ggplot(result, aes(x = Condition, fill = receiver_level_1)) + 
    geom_bar(position = "fill", width = 0.8) +  ## 通过调整 width 参数，控制条形图的宽度
    scale_fill_manual(values = colpalette) +
    theme_minimal(base_size = 16) +   ## 全局基础字体大小设定为16
    ylab("Number of receiver Cell Types")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.15, 'lines')  
    ) 
  p4 <- ggplot(result, aes(x = Condition, fill =sender_level_1)) + 
    geom_bar(position = "fill", width = 0.8) +  ## 通过调整 width 参数，控制条形图的宽度
    scale_fill_manual(values = colpalette) +
    theme_minimal(base_size = 16) +   ## 全局基础字体大小设定为16
    ylab("Number of sender Cell Types")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.15, 'lines')  
    ) 
  p5 <- ggplot(result, aes(x = Condition, fill = lr_interaction)) + 
    geom_bar(position = "fill", width = 0.8) +  ## 通过调整 width 参数，控制条形图的宽度
    theme_minimal(base_size = 16) +   ## 全局基础字体大小设定为16
    ylab("Number of lr_interaction")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.15, 'lines'),
      legend.key.size = unit(0.5, 'lines'),  # 调整图例键的大小
      legend.text = element_text(size = 5)  # 调整图例文本的大小  
    ) + guides(fill = guide_legend(ncol = 3)) 
  
  p6 <- ggplot(result, aes(x = Condition, fill = sender_receiver)) + 
    geom_bar(position = "fill", width = 0.8) +  ## 通过调整 width 参数，控制条形图的宽度
    theme_minimal(base_size = 16) +   ## 全局基础字体大小设定为16
    ylab("Number of sender_receiver")+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.15, 'lines') ,
      legend.key.size = unit(0.3, 'lines'),  # 调整图例键的大小
      legend.text = element_text(size = 3)  # 调整图例文本的大小
    )+ guides(fill = guide_legend(ncol = 3))  
  p = patchwork::wrap_plots(
    p1,p2,p3,p4,p5,p6,
    nrow = 6,
    #guides = "collect",
    widths = widths
  )
  return(p)
}


pdf("figures/1.3_Condition_CellType_receiver_and_sender_Top_10.pdf",10,30)
p=plot_func(10,3)
print(p)
dev.off()

pdf("figures/1.3_Condition_CellType_receiver_and_sender_Top_20.pdf",10,30)
p=plot_func(20,3)
print(p)
dev.off()

pdf("figures/1.3_Condition_CellType_receiver_and_sender_Top_30.pdf",10,30)
p=plot_func(30,3)
print(p)
dev.off()

pdf("figures/1.3_Condition_CellType_receiver_and_sender_Top_50.pdf",10,40)
p=plot_func(50,3)
print(p)
dev.off()

contrasts_oi = c("'A0-(A1+A2+B0+B1+B2+C0+C1+C2)/8','A1-(A0+A2+B0+B1+B2+C0+C1+C2)/8','A2-(A1+A0+B0+B1+B2+C0+C1+C2)/8','B0-(A1+A2+A0+B1+B2+C0+C1+C2)/8','B1-(A1+A2+B0+A0+B2+C0+C1+C2)/8','B2-(A1+A2+B0+B1+A0+C0+C1+C2)/8','C0-(A1+A2+B0+B1+B2+A0+C1+C2)/8','C1-(A1+A2+B0+B1+B2+C0+A0+C2)/8','C2-(A1+A2+B0+B1+B2+C0+C1+A0)/8'")
contrast_tbl = tibble(contrast = c('A0-(A1+A2+B0+B1+B2+C0+C1+C2)/8','A1-(A0+A2+B0+B1+B2+C0+C1+C2)/8','A2-(A1+A0+B0+B1+B2+C0+C1+C2)/8','B0-(A1+A2+A0+B1+B2+C0+C1+C2)/8','B1-(A1+A2+B0+A0+B2+C0+C1+C2)/8','B2-(A1+A2+B0+B1+A0+C0+C1+C2)/8','C0-(A1+A2+B0+B1+B2+A0+C1+C2)/8','C1-(A1+A2+B0+B1+B2+C0+A0+C2)/8','C2-(A1+A2+B0+B1+B2+C0+C1+A0)/8'), 
                      group = c("A0","A1","A2","B0","B1","B2","C0","C1","C2"))

prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 20, 
  senders_oi = "Mono_CD14_ATG7",
  groups_oi="A1",
  rank_per_group = FALSE
)

head(prioritized_tbl_oi_all)

nrow(prioritized_tbl_oi_all)
prioritized_tbl_oi_all

head(multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities,n=1)
head(contrast_tbl,n=1)

d <- multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
  distinct(ligand, target, direction_regulation, contrast) %>% inner_join(contrast_tbl) %>% ungroup()
nrow(d)
head(d,n=1)
lr_target_prior =lr_target_prior %>% inner_join(d)
nrow(lr_target_prior)
head(lr_target_prior,n=1)

lr_target_prior = prioritized_tbl_oi_all %>% inner_join(
  multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
    distinct(ligand, target, direction_regulation, contrast) %>% inner_join(contrast_tbl) %>% ungroup() 
) 
lr_target_df = lr_target_prior %>% distinct(group, sender, receiver, ligand, receptor, 
                                            id, target, direction_regulation,
                                            target,prioritization_rank,contrast)
nrow(lr_target_df)
head(lr_target_df,n=1)

table(lr_target_df$group)
table(lr_target_df$sender)
table(lr_target_df$receiver)
table(prioritized_tbl_oi_all$receiver)

lr_target_df_Temra <- lr_target_df%>% filter(receiver =="CD8_Temra")
head(lr_target_df_Temra,n=3) 
unique(lr_target_df_Temra $prioritization_rank)
unique(lr_target_df_Temra $target)

?infer_intercellular_regulatory_network

infer_intercellular_regulatory_network = function(lr_target_df, prioritized_tbl_oi){
  
  requireNamespace("dplyr")
  
  lr_target_df = lr_target_df %>% dplyr::inner_join(prioritized_tbl_oi, by = c("sender", "receiver", "ligand", "receptor", "id", "group"))
  
  source_df_lr = prioritized_tbl_oi %>% dplyr::mutate(
    celltype_ligand = paste(sender, ligand, sep = "_"), 
    celltype_receptor = paste(receiver, receptor, sep = "_")) %>% 
    dplyr::select(group, sender, receiver, celltype_ligand, celltype_receptor, ligand, receptor) 
  
  source_df_lrt = lr_target_df %>% dplyr::mutate(
    celltype_ligand = paste(sender, ligand, sep = "_"), 
    celltype_target = paste(receiver, target, sep = "_"), 
    celltype_receptor = paste(receiver, receptor, sep = "_")) %>% 
    dplyr::select(group, sender, receiver, celltype_ligand, celltype_receptor, celltype_target, ligand, target, receptor, direction_regulation) 
  
  lr_gr_network = dplyr::bind_rows(
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_ligand & !celltype_target %in% source_df_lr$celltype_receptor) %>% dplyr::mutate(type_target = "ligand"),
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_receptor  & !celltype_target %in% source_df_lr$celltype_ligand) %>% dplyr::mutate(type_target = "receptor")
  ) %>% dplyr::bind_rows(
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_ligand & celltype_target %in% source_df_lr$celltype_receptor) %>% dplyr::mutate(type_target = "ligand/receptor")
  )
  ligand_target_network = lr_gr_network %>% dplyr::select(celltype_ligand, celltype_target, direction_regulation, group) %>% dplyr::distinct() %>% dplyr::rename(sender_ligand = celltype_ligand, receiver_target = celltype_target) %>% dplyr::mutate(type = "Ligand-Target", weight = 1)
  
  nodes = 
    lr_gr_network %>% dplyr::select(celltype_ligand, sender, ligand) %>% dplyr::rename(celltype = sender, node = celltype_ligand, gene = ligand) %>% dplyr::mutate(type_gene = "ligand") %>% 
    dplyr::bind_rows(
      lr_gr_network %>% dplyr::select(celltype_receptor, receiver, receptor) %>% dplyr::rename(celltype = receiver, node = celltype_receptor, gene = receptor) %>% dplyr::mutate(type_gene = "receptor")
    ) %>% 
    dplyr::bind_rows(
      lr_gr_network %>% dplyr::select(celltype_target, receiver, target, type_target) %>% dplyr::rename(celltype = receiver, node = celltype_target, gene = target, type_gene = type_target)
    ) %>% dplyr::distinct() %>% 
    dplyr::filter(node %in% c(ligand_target_network$sender_ligand, ligand_target_network$receiver_target))
  
  double_nodes =  nodes %>% dplyr::group_by(node) %>% dplyr::count() %>% dplyr::filter(n > 1) %>% pull(node)
  nodes = dplyr::bind_rows(
    nodes %>% dplyr::filter(node %in% double_nodes) %>% dplyr::mutate(type_gene = "ligand/receptor") ,
    nodes %>% dplyr::filter(!node %in% double_nodes)
  ) %>% dplyr:: distinct()
  
  return(list(links = ligand_target_network, nodes = nodes, prioritized_lr_interactions = lr_gr_network %>% distinct(group, sender, receiver, ligand, receptor)))
}


network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)


network 

visualize_network_new = function(network, colors,layout = 'dh'){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggraph")
  
  nodes = network$nodes %>% data.frame() %>% magrittr::set_rownames(network$nodes$node)
  
  colors_regulation = NULL
  colors_regulation["up"] = "indianred1"
  colors_regulation["down"] = "steelblue2"
  
  
  # create the network object
  network_igraph = igraph::graph_from_data_frame(d=network$links %>% dplyr::filter(type == "Ligand-Target"), vertices = nodes, directed=T)
  network_tidygraph = tidygraph::as_tbl_graph(network_igraph) 
  set.seed(191)
  
  plot =  ggraph::ggraph(network_tidygraph, layout = layout) + 
    ggraph::geom_edge_fan(aes(color = direction_regulation), edge_width = 0.5, arrow = arrow(length = unit(1, 'mm')),
                          end_cap = ggraph::circle(1, 'mm'), 
                          start_cap = ggraph::circle(1, 'mm')) +    
    ggraph::geom_edge_loop(aes(color = direction_regulation), edge_width = 2, alpha = 0.70)  + 
    ggraph::geom_node_label(aes(label = gene, fill = celltype), fontface = "bold", 
                            size = 2, nudge_x = 0, nudge_y = 0, 
                            color = "whitesmoke") +
    ggraph::theme_graph(foreground = 'black', fg_text_colour = 'white', base_family = 'Helvetica') + 
    facet_grid(. ~group)  + 
    ggraph::scale_edge_color_manual(values = colors_regulation) + 
    scale_fill_manual(values = colors)
  
  
  return(list(plot = plot, network_igraph = network_igraph, network_tidygraph = network_tidygraph))
}

pdf('figures_point/2.1_Mono_CD14_ATG7_senders_top_20.pdf',10,4)
network_graph = visualize_network_new(network, colors_sender,'kk')
print(network_graph$plot)
dev.off()





senders_receivers = unique(multinichenet_output$celltype_info$avg_df$celltype)
color_pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(senders_receivers))
colors_sender = setNames(color_pal, senders_receivers)
colors_receiver = setNames(color_pal, senders_receivers)

generate_network_plot_sender <- function(sender_name, groups_name=NULL , top_n) {
  prioritized_tbl_oi_all = get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    top_n = top_n, 
    senders_oi = sender_name,
    groups_oi=groups_name,
    rank_per_group = FALSE
  )
  
  lr_target_prior = prioritized_tbl_oi_all %>% inner_join(
    multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
      distinct(ligand, target, direction_regulation, contrast) %>% inner_join(contrast_tbl) %>% ungroup() 
  ) 
  
  lr_target_df = lr_target_prior %>% distinct(group, sender, receiver, ligand, receptor, 
                                              id, target, direction_regulation,
                                              target,prioritization_rank,contrast) 
  network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
  network_graph = visualize_network(network, colors_sender)
  
  return(network_graph$plot)
}

# 定义组和文件名的组合
groups <- c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")

# 封装绘图函数
generate_senders_plots <- function(cell_type, output_file, top_n,width,hight) {
  pdf(output_file, width,hight)
  
  # 生成无分组的图
  tryCatch({
    p <- generate_network_plot_sender(cell_type, NULL, top_n)
    print(p)
  }, error = function(e) {
    message(paste("Error encountered with cell type", cell_type, "for ungrouped data: skipping"))
  })
  
  # 循环遍历每个组并生成相应的图
  for (group in groups) {
    tryCatch({
      p <- generate_network_plot_sender(cell_type, group, top_n)
      print(p)
    }, error = function(e) {
      message(paste("Error encountered with group", group, "for cell type", cell_type, ": skipping to next"))
    })
  }
  
  dev.off()
}

generate_senders_plots("Platelet", "figures/2.1_celltype_network_top/2.1_Platelet_senders_top_200.pdf",200,18,8)

generate_senders_plots("Platelet", "figures/2.1_celltype_network_top/2.1_Platelet_senders_top_100.pdf",100,16,6)
generate_senders_plots("CD8_Temra", "figures/2.1_celltype_network_top/2.1_CD8_Temra_senders_top_100.pdf",100,16,6)
generate_senders_plots("CD4_Treg", "figures/2.1_celltype_network_top/2.1_CD4_Treg_senders_top_100.pdf",100,16,6)
generate_senders_plots("CD8_NELL2", "figures/2.1_celltype_network_top/2.1_CD8_NELL2_senders_top_100.pdf",100,16,6)
generate_senders_plots("Mono_CD14_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD14_ATG7_senders_top_100.pdf",100,16,6)
generate_senders_plots("Mono_CD16_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD16_ATG7_senders_top_100.pdf",100,16,6)
generate_senders_plots("Mono_CD14", "figures/2.1_celltype_network_top/2.1_Mono_CD14_senders_top_100.pdf",100,16,6)
generate_senders_plots("Mono_CD16", "figures/2.1_celltype_network_top/2.1_Mono_CD16_senders_top_100.pdf",100,16,6)

generate_senders_plots("Platelet", "figures/2.1_celltype_network_top/2.1_Platelet_senders_top_50.pdf",50,14,6)
generate_senders_plots("CD8_Temra", "figures/2.1_celltype_network_top/2.1_CD8_Temra_senders_top_50.pdf",50,14,6)
generate_senders_plots("CD4_Treg", "figures/2.1_celltype_network_top/2.1_CD4_Treg_senders_top_50.pdf",50,14,6)
generate_senders_plots("CD8_NELL2", "figures/2.1_celltype_network_top/2.1_CD8_NELL2_senders_top_50.pdf",50,14,6)
generate_senders_plots("Mono_CD14_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD14_ATG7_senders_top_50.pdf",50,14,6)
generate_senders_plots("Mono_CD16_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD16_ATG7_senders_top_50.pdf",50,14,6)
generate_senders_plots("Mono_CD14", "figures/2.1_celltype_network_top/2.1_Mono_CD14_senders_top_50.pdf",50,14,6)
generate_senders_plots("Mono_CD16", "figures/2.1_celltype_network_top/2.1_Mono_CD16_senders_top_50.pdf",50,14,6)

# generate_senders_plots("Platelet", "figures/2.1_celltype_network_top/2.1_Platelet_senders_top_20.pdf")
# generate_senders_plots("CD8_Temra", "figures/2.1_celltype_network_top/2.1_CD8_Temra_senders_top_20.pdf")
# generate_senders_plots("CD4_Treg", "figures/2.1_celltype_network_top/2.1_CD4_Treg_senders_top_20.pdf")
# generate_senders_plots("CD8_NELL2", "figures/2.1_celltype_network_top/2.1_CD8_NELL2_senders_top_20.pdf")
# generate_senders_plots("Mono_CD14_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD14_ATG7_senders_top_20.pdf")
# generate_senders_plots("Mono_CD16_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD16_ATG7_senders_top_20.pdf")
generate_senders_plots("Mono_CD14", "figures/2.1_celltype_network_top/2.1_Mono_CD14_senders_top_100.pdf",100,14,6)
#generate_senders_plots("Mono_CD16", "figures/2.1_celltype_network_top/2.1_Mono_CD16_senders_top_20.pdf",20,12,4)

generate_senders_plots("Platelet", "figures/2.1_celltype_network_top/2.1_Platelet_senders_top_20.pdf",20,12,4)
generate_senders_plots("CD8_Temra", "figures/2.1_celltype_network_top/2.1_CD8_Temra_senders_top_20.pdf",20,12,4)
generate_senders_plots("CD4_Treg", "figures/2.1_celltype_network_top/2.1_CD4_Treg_senders_top_20.pdf",20,12,4)
generate_senders_plots("CD8_NELL2", "figures/2.1_celltype_network_top/2.1_CD8_NELL2_senders_top_20.pdf",20,12,4)
generate_senders_plots("Mono_CD14_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD14_ATG7_senders_top_20.pdf",20,12,4)
generate_senders_plots("Mono_CD16_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD16_ATG7_senders_top_20.pdf",20,12,4)
generate_senders_plots("Mono_CD14", "figures/2.1_celltype_network_top/2.1_Mono_CD14_senders_top_20.pdf",20,12,4)
generate_senders_plots("Mono_CD16", "figures/2.1_celltype_network_top/2.1_Mono_CD16_senders_top_20.pdf",20,12,4)

generate_network_plot <- function(receiver_name,groups_name, top_n) {
  prioritized_tbl_oi_all = get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    top_n = top_n, 
    groups_oi=groups_name,
    receivers_oi = receiver_name,
    rank_per_group = FALSE
  )
  
  lr_target_prior = prioritized_tbl_oi_all %>% inner_join(
    multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
      distinct(ligand, target, direction_regulation, contrast) %>% inner_join(contrast_tbl) %>% ungroup() 
  ) 
  
  lr_target_df = lr_target_prior %>% distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation,target,prioritization_rank,contrast) 
  network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
  network_graph = visualize_network(network, colors_sender)
  network_graph$plot
  return(network_graph$plot)
}

# 定义组和文件名的组合
groups <- c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2")

# 封装绘图函数
generate_receivers_plots <- function(cell_type, output_file) {
  pdf(output_file, 12, 4)
  
  # 生成无分组的图
  tryCatch({
    p <- generate_network_plot(cell_type, NULL, 20)
    print(p)
  }, error = function(e) {
    message(paste("Error encountered with cell type", cell_type, "for ungrouped data: skipping"))
  })
  
  # 循环遍历每个组并生成相应的图
  for (group in groups) {
    tryCatch({
      p <- generate_network_plot(cell_type, group, 20)
      print(p)
    }, error = function(e) {
      message(paste("Error encountered with group", group, "for cell type", cell_type, ": skipping to next"))
    })
  }
  
  dev.off()
}


# 调用函数生成不同的图
generate_receivers_plots("Platelet", "figures/2.1_celltype_network_top/2.1_Platelet_receivers_top_20.pdf")
generate_receivers_plots("CD8_Temra", "figures/2.1_celltype_network_top/2.1_CD8_Temra_receivers_top_20.pdf")
generate_receivers_plots("CD4_Treg", "figures/2.1_celltype_network_top/2.1_CD4_Treg_receivers_top_20.pdf")
generate_receivers_plots("CD8_NELL2", "figures/2.1_celltype_network_top/2.1_CD8_NELL2_receivers_top_20.pdf")
generate_receivers_plots("Mono_CD14_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD14_ATG7_receivers_top_20.pdf")
generate_receivers_plots("Mono_CD16_ATG7", "figures/2.1_celltype_network_top/2.1_Mono_CD16_ATG7_receivers_top_20.pdf")

prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 100, 
  # senders_oi = "Platelet",
  # receivers_oi = "CD8_Temra",
  rank_per_group = TRUE #是否需要局部排序
)

top_n_target = 20
lr_target_prior_cor_filtered = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% 
  lapply(function(group_oi){
    lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
      inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
          distinct(ligand, target, direction_regulation, contrast)
      ) %>% 
      inner_join(contrast_tbl) %>% filter(group == group_oi)
    
    lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "up") %>% 
      filter( (rank_of_target < top_n_target) & (pearson > 0.33))
    
    lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "down") %>% 
      filter( (rank_of_target < top_n_target) & (pearson < -0.33))
    lr_target_prior_cor_filtered = bind_rows(
      lr_target_prior_cor_filtered_up, 
      lr_target_prior_cor_filtered_down
    )
  }) %>% bind_rows()

lr_target_df = lr_target_prior_cor_filtered %>% 
  distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 

network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
network$links %>% head()

unique(lr_target_df$group)
head(lr_target_df)

network$nodes %>% head()

# 加载RColorBrewer包
library(RColorBrewer)

# 显示所有调色板
display.brewer.all()

colors_sender

#colors_sender["Platelet_JAM3"] = "pink" # the  original yellow background with white font is not very readable
network_graph = visualize_network(network, colors_sender)
network_graph$plot
pdf('figures/2.2_network_graph_plot_all_top_100.pdf',60,4)
network_graph$plot
dev.off()

network$prioritized_lr_interactions

prioritized_tbl_oi_network = prioritized_tbl_oi_all %>% inner_join(
  network$prioritized_lr_interactions)
prioritized_tbl_oi_network

#celltype_to_celltype
create_plots_for_celltype_sr <- function(celltype_s,celltype_r, 
                                         celltype_s_name,celltype_r_name,
                                         top_n_values, height_values) {
  for (i in seq_along(top_n_values)) {
    prioritized_tbl_oi_all = get_top_n_lr_pairs(
      multinichenet_output$prioritization_tables, 
      top_n = top_n_values[i], 
      senders_oi = celltype_s,
      receivers_oi = celltype_r,
      rank_per_group = FALSE
    )
    
    plot_oi = make_sample_lr_prod_activity_plots_Omnipath_line(
      
      filename=paste0("figures/2.3_celltype_to_celltype_top/data/2.3_", 
                      celltype_s_name, "_to_", celltype_r_name,"_top_",top_n_values[i]),
      multinichenet_output$prioritization_tables, 
      prioritized_tbl_oi_all %>% inner_join(lr_network_all)
    )
    
    pdf(paste0("figures/2.3_celltype_to_celltype_top/2.3_", 
               celltype_s_name, "_to_", celltype_r_name,"_top_",
               top_n_values[i],".pdf"), 35, height_values[i])
    print(plot_oi)
    dev.off()
  }
}


#celltype_to_celltype_in_condition
create_plots_for_celltype_sr_and_group <- function(celltype_s,celltype_r,
                                                   celltype_s_name,celltype_r_name,
                                                   group_oi, top_n_values, height_values) {
  for (i in seq_along(top_n_values)) {
    prioritized_tbl_oi_all = get_top_n_lr_pairs(
      multinichenet_output$prioritization_tables, 
      top_n = top_n_values[i], 
      senders_oi = celltype_s,
      receivers_oi = celltype_r,
      groups_oi = group_oi,
      rank_per_group = FALSE
    )
    
    plot_oi = make_sample_lr_prod_activity_plots_Omnipath_line(
      
      filename=paste0("figures/2.3_celltype_to_celltype_top/data/2.3_", 
                      celltype_s_name, "_to_", celltype_r_name,"_in_",group_oi,"_top_",top_n_values[i]),
      multinichenet_output$prioritization_tables, 
      prioritized_tbl_oi_all %>% inner_join(lr_network_all)
    )
    
    pdf(paste0("figures/2.3_celltype_to_celltype_top/2.3_", 
               celltype_s_name, "_to_", celltype_r_name,"_in_",
               group_oi,"_top_",top_n_values[i],".pdf"), 35, height_values[i])
    print(plot_oi)
    dev.off()
  }
}

# create_plots_for_celltype_and_group(celltype_s,celltype_r,group_oi, top_n_values, height_values) 
create_plots_for_celltype_sr("CD8_NELL2","ILC2","CD8_NELL2","ILC2", 20, 7) 
create_plots_for_celltype_sr_and_group("CD8_NELL2","ILC2","CD8_NELL2","ILC2","A1", 20, 7) 
create_plots_for_celltype_sr("CD8_NELL2","pDC","CD8_NELL2","pDC", 20, 7) 
create_plots_for_celltype_sr_and_group("CD8_NELL2","pDC","CD8_NELL2","pDC","A1", 20, 7) 

# create_plots_for_celltype_and_group(celltype_s,celltype_r,group_oi, top_n_values, height_values) 
create_plots_for_celltype_sr("Platelet","CD8_Temra","Platelet","CD8_Temra", 20, 7) 
create_plots_for_celltype_sr_and_group("Platelet","CD8_Temra","Platelet","CD8_Temra", "B0", 20, 7) 
create_plots_for_celltype_sr_and_group("Platelet","CD8_Temra","Platelet","CD8_Temra", "B1", 20, 7) 
create_plots_for_celltype_sr("Platelet","pDC","Platelet","pDC", 20, 7) 
create_plots_for_celltype_sr_and_group("Platelet","pDC","Platelet","pDC", "B1", 20, 7) 
create_plots_for_celltype_sr("Mono_CD14_ATG7","CD8_Temra","Mono_CD14_ATG7","CD8_Temra", 20, 7) 
create_plots_for_celltype_sr_and_group("Mono_CD14_ATG7","CD8_Temra", 
                                       "Mono_CD14_ATG7","CD8_Temra",
                                       "A1", 20, 7) 

create_plots_for_celltype_sr_and_group("Mono_CD14_ATG7",CD4T_CD8T, 
                                       "Mono_CD14_ATG7","CD4T_CD8T","C1", 20, 8) 

CD4T_CD8T <-c("CD4_Tn","CD4_Tcm","CD4_Tem","CD4_Treg","CD4_T_CCR6","CD4_Temra","CD4_Th2",
              "CD8_Tn","CD8_Tem","CD8_Temra","CD8_NELL2")
create_plots_for_celltype_sr("Mono_CD14_ATG7",CD4T_CD8T, "Mono_CD14_ATG7","CD4T_CD8T",20, 8) 
create_plots_for_celltype_sr_and_group("Mono_CD14_ATG7",CD4T_CD8T, 
                                       "Mono_CD14_ATG7","CD4T_CD8T","A1", 20, 8) 

create_plots_for_celltype_sr("Mono_CD14",CD4T_CD8T, 
                             "Mono_CD14","CD4T_CD8T",20, 8) 
create_plots_for_celltype_sr_and_group("Mono_CD14",CD4T_CD8T, 
                                       "Mono_CD14","CD4T_CD8T","A1", 20, 8) 


create_plots_for_celltype_sr("Mono_CD16",CD4T_CD8T, 
                             "Mono_CD16","CD4T_CD8T",20, 8) 
create_plots_for_celltype_sr_and_group("Mono_CD16",CD4T_CD8T, 
                                       "Mono_CD16","CD4T_CD8T","A1", 20, 8) 
create_plots_for_celltype_sr("Mono_CD16_ATG7",CD4T_CD8T, 
                             "Mono_CD16_ATG7","CD4T_CD8T",20, 8) 
create_plots_for_celltype_sr_and_group("Mono_CD16_ATG7",CD4T_CD8T,  #A1算不出来，在A1时期无信号对
                                       "Mono_CD16_ATG7","CD4T_CD8T","C1", 20, 8) 



ligands_oi = multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% 
  inner_join(contrast_tbl) %>% 
  group_by(group, receiver) %>% 
  distinct(ligand, receiver, group, activity) %>% #确保在每个组-受体组合中，每个配体的活性值都是唯一的
  top_n(1, activity) %>% #对每个组-受体组合选择活性值最高的配体
  pull(ligand) %>% unique()
plot_oi = make_ligand_activity_plots(
  multinichenet_output$prioritization_tables, 
  ligands_oi, 
  contrast_tbl,
  widths = NULL)
pdf("figures/3_ligand_gene/3_all_sender_top1_plot_oi.pdf",30,20)
plot_oi
dev.off()

ligands_oi_B1 = multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == 'B1')%>% 
  group_by(group, receiver) %>% 
  distinct(ligand, receiver, group, activity) %>% #确保在每个组-受体组合中，每个配体的活性值都是唯一的
  top_n(20, activity) %>% #对每个组-受体组合选择活性值最高的配体
  pull(ligand) %>% unique()
plot_oi = make_ligand_activity_plots(
  multinichenet_output$prioritization_tables, 
  ligands_oi_B1, 
  contrast_tbl,
  widths = NULL)
pdf("figures/3_ligand_gene/3_B1_receiver_top20.pdf",30,6)
plot_oi
dev.off()

ligands_oi_TGFB1 = multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% 
  inner_join(contrast_tbl) %>% 
  filter(ligand == 'TGFB1')%>% 
  group_by(group, receiver) %>% 
  distinct(ligand, receiver, group, activity) %>% #确保在每个组-受体组合中，每个配体的活性值都是唯一的
  top_n(20, activity) %>% #对每个组-受体组合选择活性值最高的配体
  pull(ligand) %>% unique()
plot_oi = make_ligand_activity_plots(
  multinichenet_output$prioritization_tables, 
  ligands_oi_TGFB1, 
  contrast_tbl,
  widths = NULL)
pdf("figures/3_ligand_gene/3_TGFB1_receiver_top20.pdf",30,2.3)

plot_oi
dev.off()

prioritized_tbl_oi_M_20 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  20, 
  groups_oi = "B0", 
  receivers_oi = "CD8_Temra",
  senders_oi = "Platelet")

combined_plot = make_ligand_activity_target_plot(
  group_oi =group_oi, 
  receiver_oi=c("CD8_Temra"), 
  prioritized_tbl_oi_M_20,
  multinichenet_output$prioritization_tables, 
  multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  ligand_target_matrix, 
  heights=c(3,1),
  widths=c(1,1,3),
  plot_legend = FALSE)
pdf("figures/4.1_ligand_target_gene/4_B0_Platelet_to_CD8_Temra.pdf",17,5.5)
combined_plot
dev.off()

group_oi = "A1"
senders_oi = "CD8_NELL2"
prioritized_tbl_oi_M_10 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  20, 
  groups_oi = group_oi, 
  senders_oi = senders_oi)

combined_plot = make_ligand_activity_target_plot(
  group_oi, 
  "ILC2", 
  prioritized_tbl_oi_M_10,
  multinichenet_output$prioritization_tables, 
  multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  ligand_target_matrix, 
  heights=c(3,1),
  widths=c(1,1,6),
  plot_legend = FALSE)
pdf("figures/4.1_ligand_target_gene/4_A1_CD8_NELL2_to_ILC2.pdf",25,5)
combined_plot
dev.off()

group_oi = "A1"
senders_oi = "Mono_CD14_ATG7"
prioritized_tbl_oi_M_10 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  10, 
  groups_oi = group_oi, 
  senders_oi = senders_oi)

combined_plot = make_ligand_activity_target_plot(
  group_oi, 
  "CD8_Temra", 
  prioritized_tbl_oi_M_10,
  multinichenet_output$prioritization_tables, 
  multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  ligand_target_matrix, 
  heights=c(3,1),
  widths=c(1,1,4),
  plot_legend = FALSE)
pdf("figures/4.1_ligand_target_gene/4_A1_Mono_CD14_ATG7_to_CD8_Temra.pdf",22,5.5)
combined_plot
dev.off()

make_lr_target_scatter_plot = function(prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, receiver_info, grouping_tbl, lr_target_prior_cor_filtered){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  #library(viridis)
  library(ggsci)
  
  grouping_tbl = grouping_tbl %>% dplyr::inner_join(prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, keep_receiver, keep_sender))
  
  ligand_receptor_pb_prod_df = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(keep_receiver == 1 & keep_sender == 1) %>% dplyr::filter(ligand == ligand_oi, receptor == receptor_oi, sender == sender_oi, receiver == receiver_oi) %>% dplyr::distinct(sample, ligand_receptor_pb_prod, id)
  
  targets_oi = ligand_receptor_pb_prod_df %>% dplyr::inner_join(lr_target_prior_cor_filtered, by = "id") %>% dplyr::pull(target) %>% unique()
  
  target_pb_df = receiver_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = "sample") %>% dplyr::filter(keep_receiver == 1 & keep_sender == 1) %>% dplyr::filter(gene %in% targets_oi & celltype == receiver_oi) %>% dplyr::distinct(gene, sample, pb_sample) %>% dplyr::rename(target_pb = pb_sample)
  
  ligand_receptor_target_pb_df = grouping_tbl %>% dplyr::inner_join(ligand_receptor_pb_prod_df, by = "sample") %>% dplyr::inner_join(target_pb_df, by = "sample")
  
  p = ligand_receptor_target_pb_df %>% ggplot(aes(ligand_receptor_pb_prod, target_pb)) + 
    geom_point(aes(color = group), size = 2) + 
    geom_smooth(method = "lm",alpha = 0.10, color = "grey50") + 
    facet_grid(.~gene) +
    theme_bw() + 
    #scale_color_brewer(palette = "Set2") + 
    #scale_color_viridis(discrete = TRUE, option = "D") +
    scale_color_manual(values = ggsci::pal_jco("default")(10))+
    ggtitle(paste(ligand_oi,"-",receptor_oi," expression vs target expression between ", sender_oi," and ", receiver_oi, " cells.", sep = "")) +
    xlab("LR pair pseudobulk expression product") + ylab("Target gene pseudobulk expression")
  return(p)
} 

group_oi = "A1"
receiver_oi = "CD8_NELL2"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
  inner_join(
    multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% 
      distinct(ligand, target, direction_regulation, contrast)
  ) %>% 
  inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)

lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "up") %>% 
  filter( (rank_of_target < top_n_target) & (pearson > 0.33)) # replace pearson by spearman if you want to filter on the spearman correlation
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "down") %>% 
  filter( (rank_of_target < top_n_target) & (pearson < -0.33)) # downregulation -- negative correlation - # replace pearson by spearman if you want to filter on the spearman correlation
lr_target_prior_cor_filtered = bind_rows(
  lr_target_prior_cor_filtered_up, 
  lr_target_prior_cor_filtered_down)

prioritized_tbl_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  20, 
  groups_oi = group_oi, 
  receivers_oi = receiver_oi)

lr_target_correlation_plot = make_lr_target_correlation_plot(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi,  
  lr_target_prior_cor_filtered , 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  receiver_oi,
  plot_legend = FALSE)
pdf("figures/4.2_A1_CD8_NELL2_receiver_top_20_correlation_plot.pdf",50,30)
lr_target_correlation_plot$combined_plot
dev.off()

ligand_oi = "TNFSF10"
receptor_oi = "TNFRSF10D"
sender_oi = "Mono_CD16"
receiver_oi = "CD8_NELL2"

lr_target_scatter_plot = make_lr_target_scatter_plot(
  multinichenet_output$prioritization_tables, 
  ligand_oi, receptor_oi, sender_oi, receiver_oi, 
  multinichenet_output$celltype_info, 
  multinichenet_output$grouping_tbl, 
  lr_target_prior_cor_filtered)

pdf("figures/4.2_A1_CD8_NELL2_TNFSF10_lr_target_cor.pdf",30,4)
lr_target_scatter_plot
dev.off()

group_oi = "B1"
senders_oi = "Platelet"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
  inner_join(
    multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% 
      distinct(ligand, target, direction_regulation, contrast)
  ) %>% 
  inner_join(contrast_tbl) %>% filter(group == group_oi, sender == senders_oi)

lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "up") %>% 
  filter( (rank_of_target < top_n_target) & (pearson > 0.33)) # replace pearson by spearman if you want to filter on the spearman correlation
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
  filter(direction_regulation == "down") %>% 
  filter( (rank_of_target < top_n_target) & (pearson < -0.33)) # downregulation -- negative correlation - # replace pearson by spearman if you want to filter on the spearman correlation
lr_target_prior_cor_filtered = bind_rows(
  lr_target_prior_cor_filtered_up, 
  lr_target_prior_cor_filtered_down)

prioritized_tbl_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  200, 
  groups_oi = group_oi, 
  senders_oi = senders_oi)

lr_target_correlation_plot = make_lr_target_correlation_plot(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi,  
  lr_target_prior_cor_filtered , 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  senders_oi,
  plot_legend = FALSE)
pdf("figures/4.2_B1_Platelet_senders_top_200_correlation_plot.pdf",50,35)
lr_target_correlation_plot$combined_plot
dev.off()

ligand_oi = "TGFB1"
receptor_oi = "ITGB1"
sender_oi = "Platelet"
receiver_oi = "Platelet"

lr_target_scatter_plot = make_lr_target_scatter_plot(
  multinichenet_output$prioritization_tables, 
  ligand_oi, receptor_oi, sender_oi, receiver_oi, 
  multinichenet_output$celltype_info, 
  multinichenet_output$grouping_tbl, 
  lr_target_prior_cor_filtered)

pdf("figures/4.2_B1_Platelet_TGFB1_lr_target_scatter_plot.pdf",30,4)
lr_target_scatter_plot
dev.off()

if(organism == "human"){
  sig_network = readRDS("../files/signaling_network_human_21122021.rds") %>% 
    mutate(from = make.names(from), to = make.names(to))
  
  gr_network = readRDS("../files/gr_network_human_21122021.rds") %>% 
    mutate(from = make.names(from), to = make.names(to))
  
  ligand_tf_matrix = readRDS("../files/ligand_tf_matrix_nsga2r_final.rds")
  colnames(ligand_tf_matrix) = colnames(ligand_tf_matrix) %>% make.names()
  rownames(ligand_tf_matrix) = rownames(ligand_tf_matrix) %>% make.names()
  
  weighted_networks = readRDS("../files/weighted_networks_nsga2r_final.rds")
  weighted_networks$lr_sig = weighted_networks$lr_sig %>% mutate(from = make.names(from), to = make.names(to))
  weighted_networks$gr = weighted_networks$gr %>% mutate(from = make.names(from), to = make.names(to))
  
} else if(organism == "mouse"){
  sig_network = readRDS("files/signaling_network_mouse_21122021.rds") %>% 
    mutate(from = make.names(from), to = make.names(to))
  
  gr_network = readRDS("files/gr_network_mouse_21122021.rds") %>% 
    mutate(from = make.names(from), to = make.names(to))
  
  ligand_tf_matrix = readRDS("files/ligand_tf_matrix_nsga2r_final_mouse.rds")
  colnames(ligand_tf_matrix) = colnames(ligand_tf_matrix) %>% make.names()
  rownames(ligand_tf_matrix) = rownames(ligand_tf_matrix) %>% make.names()
  
  weighted_networks = readRDS("files/weighted_networks_nsga2r_final_mouse.rds")
  weighted_networks$lr_sig = weighted_networks$lr_sig %>% mutate(from = make.names(from), to = make.names(to))
  weighted_networks$gr = weighted_networks$gr %>% mutate(from = make.names(from), to = make.names(to))
}

prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 200, 
  rank_per_group = FALSE #是否需要局部排序
)

top_n_target = 200
lr_target_prior_cor_filtered = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% 
  lapply(function(group_oi){
    lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
      inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
          distinct(ligand, target, direction_regulation, contrast)
      ) %>% 
      inner_join(contrast_tbl) %>% filter(group == group_oi)
    
    lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "up") %>% 
      filter( (rank_of_target < top_n_target) & (pearson > 0.33))
    
    lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "down") %>% 
      filter( (rank_of_target < top_n_target) & (pearson < -0.33))
    lr_target_prior_cor_filtered = bind_rows(
      lr_target_prior_cor_filtered_up, 
      lr_target_prior_cor_filtered_down
    )
  }) %>% bind_rows()

lr_target_df = lr_target_prior_cor_filtered %>% 
  distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 

network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
network$links %>% head()

network$links %>% filter(sender_ligand == "Platelet_TGFB1" & direction_regulation == "up" & group == "B1")

network$links %>% filter(sender_ligand == "Platelet_TGFB1" & direction_regulation == "up" & group == "B0")

network$links %>% filter(sender_ligand == "Platelet_PF4" )

network$links %>% filter(sender_ligand == "CD8_NELL2_TNFRSF10"  )

ligand_oi = c("PF4")
receptor_oi = c("CXCR3")
targets_all = c("CXCR3") 

active_signaling_network = nichenetr::get_ligand_signaling_path_with_receptor(
  ligand_tf_matrix = ligand_tf_matrix, 
  ligands_all = ligand_oi, 
  receptors_all = receptor_oi, 
  targets_all = targets_all, 
  weighted_networks = weighted_networks, 
  top_n_regulators = 3
)

data_source_network = nichenetr::infer_supporting_datasources(
  signaling_graph_list = active_signaling_network,
  lr_network = lr_network %>% dplyr::rename(from = ligand, to = receptor), 
  sig_network = sig_network, 
  gr_network = gr_network
)

active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
colors = c("ligand" = "purple", "receptor" = "orange", "target" = "royalblue", "mediator" = "grey60")
#ggraph_signaling_path = suppressWarnings(make_ggraph_signaling_path(active_signaling_network_min_max, colors, ligand_oi, receptor_oi, targets_all))
ggraph_signaling_path = make_ggraph_signaling_path(
  active_signaling_network_min_max, 
  colors, 
  ligand_oi, 
  receptor_oi, 
  targets_all)
#ggraph_signaling_path$plot
pdf("figures/4.3_Platelet_PF4_to_CXCR3.pdf",10,6)
ggraph_signaling_path$plot
dev.off()

ligand_oi = c("ICAM2")
receptor_oi = c("ITGB2")
targets_all = c("ITGB2") #,"ITGB1"

active_signaling_network = nichenetr::get_ligand_signaling_path_with_receptor(
  ligand_tf_matrix = ligand_tf_matrix, 
  ligands_all = ligand_oi, 
  receptors_all = receptor_oi, 
  targets_all = targets_all, 
  weighted_networks = weighted_networks, 
  top_n_regulators = 3
)

data_source_network = nichenetr::infer_supporting_datasources(
  signaling_graph_list = active_signaling_network,
  lr_network = lr_network %>% dplyr::rename(from = ligand, to = receptor), 
  sig_network = sig_network, 
  gr_network = gr_network
)

active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
colors = c("ligand" = "purple", "receptor" = "orange", "target" = "royalblue", "mediator" = "grey60")
#ggraph_signaling_path = suppressWarnings(make_ggraph_signaling_path(active_signaling_network_min_max, colors, ligand_oi, receptor_oi, targets_all))
ggraph_signaling_path = make_ggraph_signaling_path(
  active_signaling_network_min_max, 
  colors, 
  ligand_oi, 
  receptor_oi, 
  targets_all)
#ggraph_signaling_path$plot
pdf("figures/4.3_Platelet_ICAM2_to_ITGB2_CD8_Temra.pdf",10,6)
ggraph_signaling_path$plot
dev.off()

# ligand_oi = c("TGFB1")
# receptor_oi = c("TGFBR3","ITGB1")
# targets_all = c("TGFBR3","ITGB1") 

ligand_oi = c("TGFB1")
receptor_oi = c("TGFBR3")
targets_all = c("TGFBR3") #,"ITGB1"

active_signaling_network = nichenetr::get_ligand_signaling_path_with_receptor(
  ligand_tf_matrix = ligand_tf_matrix, 
  ligands_all = ligand_oi, 
  receptors_all = receptor_oi, 
  targets_all = targets_all, 
  weighted_networks = weighted_networks, 
  top_n_regulators = 3
)

data_source_network = nichenetr::infer_supporting_datasources(
  signaling_graph_list = active_signaling_network,
  lr_network = lr_network %>% dplyr::rename(from = ligand, to = receptor), 
  sig_network = sig_network, 
  gr_network = gr_network
)

active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
colors = c("ligand" = "purple", "receptor" = "orange", "target" = "royalblue", "mediator" = "grey60")
#ggraph_signaling_path = suppressWarnings(make_ggraph_signaling_path(active_signaling_network_min_max, colors, ligand_oi, receptor_oi, targets_all))
ggraph_signaling_path = make_ggraph_signaling_path(
  active_signaling_network_min_max, 
  colors, 
  ligand_oi, 
  receptor_oi, 
  targets_all)
#ggraph_signaling_path$plot
pdf("figures/4.3_Platelet_TGFB1_to_TGFBR3.pdf",10,6)
ggraph_signaling_path$plot
dev.off()

seurat_obj<- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/5_new_level_1_2/5_merge_data_new_level_1.rds")

sce = Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")

sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

sce$ct_level_4 <- make.names(sce$ct_level_4)

# 解决Warning message in RColorBrewer::brewer.pal(n, pal):
# “n too large, allowed maximum for palette Set2 is 8
# Returning the palette you asked for with that many colors
#让C2时期也被画上颜色

##./R/plotting.R

make_ligand_receptor_violin_plot = function(sce, ligand_oi, receptor_oi, sender_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id, batch_oi = NA, background_groups = NULL){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  library(viridis)
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }
  
  # ligand plot
  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% sender_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]
  
  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  coldata_df = SummarizedExperiment::colData(sce_subset)
  if(! "cell" %in% colnames(coldata_df)){
    coldata_df = coldata_df %>% data.frame() %>% tibble::rownames_to_column("cell")
  }
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[ligand_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[ligand_oi,])) %>% dplyr::inner_join(coldata_df %>% tibble::as_tibble())
  
  if(is.na(batch_oi)){
    
    p_sender = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(.~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + 
      #scale_color_brewer(palette = "Set2") + 
      scale_color_viridis(discrete = TRUE, option = "D") +
      ggtitle(paste("Expression of the ligand ",ligand_oi, " in sender cell type ", sender_oi, sep = ""))
    
  } else {
    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)
    
    p_sender = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(batch_oi ~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + 
      #scale_color_brewer(palette = "Set2") + 
      scale_color_viridis(discrete = TRUE, option = "D") +
      ggtitle(paste("Expression of the ligand ",ligand_oi, " in sender cell type ", sender_oi, sep = ""))
    
  }
  
  # receptor plot
  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% receiver_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]
  
  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  coldata_df = SummarizedExperiment::colData(sce_subset)
  if(! "cell" %in% colnames(coldata_df)){
    coldata_df = coldata_df %>% data.frame() %>% tibble::rownames_to_column("cell")
  }
  
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[receptor_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[receptor_oi,])) %>% dplyr::inner_join(coldata_df %>% tibble::as_tibble())
  
  if(is.na(batch_oi)){
    p_receiver = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(.~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) +  
      #scale_color_brewer(palette = "Set2") + 
      scale_color_viridis(discrete = TRUE, option = "D") +
      ggtitle(paste("Expression of the receptor ",receptor_oi, " in receiver cell type ", receiver_oi, sep = ""))
    
  } else {
    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)
    
    p_receiver = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(batch_oi~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + 
      #scale_color_brewer(palette = "Set2") + 
      scale_color_viridis(discrete = TRUE, option = "D") +
      ggtitle(paste("Expression of the receptor ",receptor_oi, " in receiver cell type ", receiver_oi, sep = ""))
    
  }
  
  return(patchwork::wrap_plots(p_sender, p_receiver, nrow = 2))
  
}

make_target_violin_plot = function(sce, target_oi, receiver_oi, group_oi, group_id, sample_id, celltype_id, batch_oi = NA, background_groups = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  library(viridis)
  
  if(is.null(background_groups)){
    background_groups = SummarizedExperiment::colData(sce)[,group_id] %>% unique() %>% generics::setdiff(group_oi)
  }
  
  sce_subset =  sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% receiver_oi]
  sce_subset =  sce_subset[, SummarizedExperiment::colData(sce_subset)[,group_id] %in% c(group_oi,background_groups)]
  
  sce_subset$id = sce_subset[[sample_id]]
  sce_subset = muscat::prepSCE(sce_subset,
                               kid = celltype_id, # subpopulation assignments
                               gid = group_id,  # group IDs (ctrl/stim)
                               sid = "id",   # sample IDs (ctrl/stim.1234)
                               drop = FALSE)  #
  
  coldata_df = SummarizedExperiment::colData(sce_subset)
  if(! "cell" %in% colnames(coldata_df)){
    coldata_df = coldata_df %>% data.frame() %>% tibble::rownames_to_column("cell")
  }
  
  exprs_df = tibble::tibble(expression = SingleCellExperiment::logcounts(sce_subset)[target_oi,], cell = names(SingleCellExperiment::logcounts(sce_subset)[target_oi,])) %>% dplyr::inner_join(coldata_df %>% tibble::as_tibble())
  
  if(is.na(batch_oi)){
    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(.~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) +
      #scale_color_brewer(palette = "Set2") + 
      scale_color_viridis(discrete = TRUE, option = "D") +
      ggtitle(paste("Expression of the target ",target_oi, " in receiver cell type ", receiver_oi, sep = ""))
  } else {
    
    extra_metadata = SummarizedExperiment::colData(sce_subset) %>% tibble::as_tibble() %>% dplyr::select(all_of(sample_id), all_of(batch_oi)) %>% dplyr::distinct() %>% dplyr::mutate_all(factor)
    colnames(extra_metadata) = c("sample_id","batch_oi")
    exprs_df = exprs_df %>% dplyr::inner_join(extra_metadata)
    
    p_violin = exprs_df %>% ggplot(aes(sample_id, expression, group = sample_id, color = group_id)) + geom_violin(color = "grey10") + ggbeeswarm::geom_quasirandom(bandwidth = 1.25, varwidth = TRUE)  + 
      facet_grid(batch_oi~ group_id, scales = "free", space = "free") +
      scale_x_discrete(position = "bottom") +
      theme_light() +
      theme(
        axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
        strip.text.x.top = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.5, "lines"),
        panel.spacing.y = unit(0.25, "lines"),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0),
        strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
      ) + 
      #scale_color_brewer(palette = "Set2") +
      scale_color_viridis(discrete = TRUE, option = "D") +
      ggtitle(paste("Expression of the target ",target_oi, " in receiver cell type ", receiver_oi, sep = ""))
  }
  
  return(p_violin)
  
}

ligand_oi = "EFNB2"#"CD47"
receptor_oi = "EPHA4"
sender_oi = "CD8_NELL2"
receiver_oi = "ILC2"
group_oi = "A1"
sample_id = "Sample"
group_id = "Condition"
celltype_id = "ct_level_4"
targets_all=c("EPHA4")


p_violin = make_ligand_receptor_violin_plot(
  sce = sce, 
  ligand_oi = ligand_oi,
  receptor_oi = receptor_oi, 
  group_oi = group_oi, 
  group_id = group_id, 
  sender_oi = sender_oi, 
  receiver_oi = receiver_oi, 
  sample_id = sample_id, 
  celltype_id = celltype_id)
pdf("figures/4.4_scRNA_gene/4.4_CD8_NELL2_to_ILC2.pdf",14,6)
p_violin
dev.off()


ligand_oi = "EFNB2"#"CD47"
receptor_oi = "PECAM1"
sender_oi = "CD8_NELL2"
receiver_oi = "pDC"
group_oi = "A1"
sample_id = "Sample"
group_id = "Condition"
celltype_id = "ct_level_4"
targets_all=c("PECAM1","EPHB1","C19orf38")

p_violin = make_ligand_receptor_violin_plot(
  sce = sce, 
  ligand_oi = ligand_oi,
  receptor_oi = receptor_oi, 
  group_oi = group_oi, 
  group_id = group_id, 
  sender_oi = sender_oi, 
  receiver_oi = receiver_oi, 
  sample_id = sample_id, 
  celltype_id = celltype_id)
pdf("figures/4.4_scRNA_gene/4.4_CD8_NELL2_to_pDC.pdf",15,5)
p_violin
dev.off()

#放大特定的配体-目标相互作用：目标基因单细胞表达小提琴图
list_target_plots = lapply(targets_all, function(target_oi) {
  p = make_target_violin_plot(sce = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id = celltype_id)
})

pdf("figures/4.4_scRNA_gene/4.4_CD8_NELL2_to_pDC_target.pdf",15,2.5)
list_target_plots
dev.off()

ligand_oi = "ADAM15"
receptor_oi = "ITGA5"
sender_oi = "Mono_CD14_ATG7"
receiver_oi = "CD8_Temra"
group_oi = "A1"
sample_id = "Sample"
group_id = "Condition"
celltype_id = "ct_level_4"
#targets_all=c("ITGA5","TSPAN17","TSPAN14")

p_violin = make_ligand_receptor_violin_plot(
  sce = sce, 
  ligand_oi = ligand_oi,
  receptor_oi = receptor_oi, 
  group_oi = group_oi, 
  group_id = group_id, 
  sender_oi = sender_oi, 
  receiver_oi = receiver_oi, 
  sample_id = sample_id, 
  celltype_id = celltype_id)
pdf("figures/4.4_scRNA_gene/4.4_Mono_CD14_ATG7_to_CD8_Temra.pdf",15,5)
p_violin
dev.off()

# #放大特定的配体-目标相互作用：目标基因单细胞表达小提琴图
# list_target_plots = lapply(targets_all, function(target_oi) {
#   p = make_target_violin_plot(sce = sce, 
#                               target_oi = target_oi, receiver_oi = receiver_oi,
#                               group_oi = group_oi, 
#                               group_id = group_id, sample_id, 
#                               celltype_id = celltype_id)
# })

# pdf("figures/4.4_scRNA_gene/4.4_Mono_CD14_ATG7_to_CD8_Temra_target.pdf",15,2.5)
# list_target_plots
# dev.off()

ligand_oi = "TGFB1"
receptor_oi = "TGFBR3"
sender_oi = "Platelet"
receiver_oi = "CD8_Temra"
group_oi = "B0"
sample_id = "Sample"
group_id = "Condition"
celltype_id = "ct_level_4"
targets_all=c("TGFBR3","TGFBR1","ITGAV","ENG","APP")

p_violin = make_ligand_receptor_violin_plot(
  sce = sce, 
  ligand_oi = ligand_oi,
  receptor_oi = receptor_oi, 
  group_oi = group_oi, 
  group_id = group_id, 
  sender_oi = sender_oi, 
  receiver_oi = receiver_oi, 
  sample_id = sample_id, 
  celltype_id = celltype_id)
pdf("figures/4.4_scRNA_gene/4.4_Platelet_to_CD8_Temra_B0.pdf",15,5)
p_violin
dev.off()

#放大特定的配体-目标相互作用：目标基因单细胞表达小提琴图
list_target_plots = lapply(targets_all, function(target_oi) {
  p = make_target_violin_plot(sce = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id = celltype_id)
})

pdf("figures/4.4_scRNA_gene/4.4_Platelet_to_CD8_Temra_target_B0.pdf",15,2.5)
list_target_plots
dev.off()

ligand_oi = "TGFB1"
receptor_oi = "TGFBR3"
sender_oi = "Platelet"
receiver_oi = "CD8_Temra"
group_oi = "B1"
sample_id = "Sample"
group_id = "Condition"
celltype_id = "ct_level_4"
targets_all=c("TGFBR3","TGFBR1","ITGAV","ENG","APP")

p_violin = make_ligand_receptor_violin_plot(
  sce = sce, 
  ligand_oi = ligand_oi,
  receptor_oi = receptor_oi, 
  group_oi = group_oi, 
  group_id = group_id, 
  sender_oi = sender_oi, 
  receiver_oi = receiver_oi, 
  sample_id = sample_id, 
  celltype_id = celltype_id)
pdf("figures/4.4_scRNA_gene/4.4_Platelet_to_CD8_Temra_B1.pdf",15,5)
p_violin
dev.off()

#放大特定的配体-目标相互作用：目标基因单细胞表达小提琴图
list_target_plots = lapply(targets_all, function(target_oi) {
  p = make_target_violin_plot(sce = sce, target_oi = target_oi, receiver_oi = receiver_oi, group_oi = group_oi, group_id = group_id, sample_id, celltype_id = celltype_id)
})

pdf("figures/4.4_scRNA_gene/4.4_Platelet_to_CD8_Temra_target_B1.pdf",15,2.5)
list_target_plots
dev.off()

group_oi = "A1"
receiver_oi = "CD8_NELL2"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 2 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()

p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)
pdf("figures/4.5_CD8_NELL2_A1_target_pseudobulk_plot.pdf",20,12)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

group_oi = "B0"
receiver_oi = "Platelet"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 1 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()
DE_genes = DE_genes %>% intersect(lr_network$ligand)
p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)

pdf("figures/4.5_Platelet_B0_p_target_pseudobulk_plot.pdf",18,10)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

group_oi = "B1"
receiver_oi = "Platelet"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 1 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()
DE_genes = DE_genes %>% intersect(lr_network$ligand)
p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)

pdf("figures/4.5_Platelet_B1_p_target_pseudobulk_plot.pdf",18,10)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

#受体
group_oi = "A1"
receiver_oi = "Mono_CD14"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 1 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()
DE_genes = DE_genes %>% intersect(lr_network$receptor)
p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)
pdf("figures/4.5_Mono_CD14_A1_p_target_pseudobulk_plot.pdf",20,12)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

#受体
group_oi = "A1"
receiver_oi = "Mono_CD14_ATG7"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 1 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()
DE_genes = DE_genes %>% intersect(lr_network$receptor)
p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)
pdf("figures/4.5_Mono_CD14_ATG7_A1_p_target_pseudobulk_plot.pdf",20,12)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

#受体
group_oi = "A1"
receiver_oi = "Mono_CD16_ATG7"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 1 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()
DE_genes = DE_genes %>% intersect(lr_network$receptor)
p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)
pdf("figures/4.5_Mono_CD16_ATG7_A1_p_target_pseudobulk_plot.pdf",20,12)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

#受体
group_oi = "B1"
receiver_oi = "CD8_Temra"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 1 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()
DE_genes = DE_genes %>% intersect(lr_network$receptor)
p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)
pdf("figures/4.5_CD8_Temra_B1_p_target_pseudobulk_plot.pdf",20,12)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

#受体
group_oi = "C1"
receiver_oi = "CD8_Temra"
DE_genes = multinichenet_output$ligand_activities_targets_DEgenes$de_genes_df %>% 
  inner_join(contrast_tbl) %>% 
  filter(group == group_oi) %>% 
  arrange(p_val) %>% 
  filter(
    receiver == receiver_oi & 
      logFC > 1 & 
      p_val <= 0.05 &
      contrast == contrast_tbl %>% filter(group == group_oi) %>% pull(contrast)) %>% 
  pull(gene) %>% unique()
DE_genes = DE_genes %>% intersect(lr_network$receptor)
p_target = make_DEgene_dotplot_pseudobulk(
  genes_oi = DE_genes, 
  celltype_info = multinichenet_output$celltype_info, 
  prioritization_tables = multinichenet_output$prioritization_tables, 
  celltype_oi = receiver_oi, 
  multinichenet_output$grouping_tbl)
pdf("figures/4.5_CD8_Temra_C1_p_target_pseudobulk_plot.pdf",20,12)
p_target$pseudobulk_plot + ggtitle("DE genes (pseudobulk expression)")
p_target$singlecell_plot + ggtitle("DE ligands (single-cell expression)")
dev.off()

library(RColorBrewer)
# display.brewer.all()
library("circlize")

# 使用 colorRampPalette 生成足够的颜色
# senders_receivers <- union(prioritized_tbl_oi_all$sender %>% unique(), prioritized_tbl_oi_all$receiver %>% unique()) %>%
#   sort()
# colors_sender <- colorRampPalette(brewer.pal(11, "Spectral"))(length(senders_receivers)) %>%
#   magrittr::set_names(senders_receivers)

# colors_receiver <- colorRampPalette(brewer.pal(11, "Spectral"))(length(senders_receivers)) %>%
#   magrittr::set_names(senders_receivers)
# colors_sender
# colors_receiver 

#celltypel3
ncols <- c(
  'Bm_CD27_IGHM_D'="#E5D2DD",
  'Bm_CD27_IGHM_H'="#F3B1A0",
  'Bm_CD27_IGHM_SOX5'="#E59CC4",
  Bmn_IGHD="#AB3282",
  CD4_T_CCR6="#CCE0F5",
  CD4_Tcm="#C1E6F3",
  CD4_Tem="#91D0BE",
  CD4_Temra="#5F3D69",
  CD4_Th2="#58A4C3",
  CD4_Tn="#6778AE",
  CD4_Treg="#476D87",
  CD8_NELL2="#23452F",
  CD8_Tem="#C5DEBA",
  CD8_Temra="#3A6963",
  CD8_Tn="#D6E7A9",
  cDC1="#E0D4CA",
  cDC2="#9FA3A8",
  gdT="#53A85F",
  gdT_V9="#57C3F3",
  HPSC="#585658",
  ILC="#F1CC92",
  ILC2="#F1BC62",
  MAIT="#85B22CFF",
  Mono_CD14="#A01319",
  'Mono_CD14_ATG7'="#712820",
  Mono_CD14_CD16="#B53E2B",
  Mono_CD16="#E95C39",
  'Mono_CD16_ATG7'="#A05401",
  NK1="#E4C955",
  NK2="#E1A111",
  NK3="#AA9A59",
  pDC="#968175",
  Platelet="#942d8d",
  Plasma="#8C549C",
  Tdn="#585658"
)
# 使用 colorRampPalette 生成足够的颜色

colors_sender <- ncols 

colors_receiver <- ncols 
colors_sender
colors_receiver 

# 定义一个函数封装重复逻辑
process_circos <- function(group_name, sender, immune_cell, top_n, colors_sender, colors_receiver) {
  prioritized_tbl_oi_all <- get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    groups_oi = group_name,
    senders_oi = sender,
    receivers_oi = setdiff(immune_cell, sender),  # 过滤掉自身
    top_n = top_n,
    rank_per_group = FALSE
  )
  print(nrow(prioritized_tbl_oi_all))  # 打印行数
  make_circos_group_comparison_2(prioritized_tbl_oi_all, colors_sender, colors_receiver)
}


make_circos_group_comparison_2 = function(prioritized_tbl_oi, colors_sender, colors_receiver){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() # if grouped: things will be messed up downstream
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  # Make the plot for each condition
  groups_oi = prioritized_tbl_oi$group %>% unique()
  all_plots = groups_oi %>% lapply(function(group_oi){
    
    # Make the plot for condition of interest - title of the plot
    title = group_oi
    circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
    
    # deal with duplicated sector names
    # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
    # only do it with duplicated ones!
    # this is in iterative process because changing one ligand to try to solve a mistake here, can actually induce another problem
    circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
    
    df = circos_links
    
    
    ligand.uni = unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
      df.i = df[df$ligand == ligand.uni[i], ]
      sender.uni = unique(df.i$sender)
      for (j in 1:length(sender.uni)) {
        df.i.j = df.i[df.i$sender == sender.uni[j], ]
        df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
        df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
      }
    }
    receptor.uni = unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
      df.i = df[df$receptor == receptor.uni[i], ]
      receiver.uni = unique(df.i$receiver)
      for (j in 1:length(receiver.uni)) {
        df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
        df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
        df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
      }
    }
    
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    
    while(length(intersecting_ligands_receptors) > 0){
      df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
      df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
      df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
      df = dplyr::bind_rows(df_unique, df_duplicated)
      intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
    }
    
    circos_links = df
    
    # Link ligands/Receptors to the colors of senders/receivers
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
    links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
    grid_col =c(grid_ligand_color,grid_receptor_color)
    
    # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
    transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    
    # Define order of the ligands and receptors and the gaps
    ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
      ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    }) %>% unlist()
    
    receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
      receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
    }) %>% unlist()
    
    order = c(ligand_order,receptor_order)
    
    width_same_cell_same_ligand_type = 0.275
    width_different_cell = 3
    width_ligand_receptor = 9
    width_same_cell_same_receptor_type = 0.275
    
    # sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    #   sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # sender_gaps = sender_gaps[-length(sender_gaps)]
    
    sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% 
      lapply(function(sender_oi){
        # 计算当前 sender 对应的配体数量
        num_ligands = circos_links %>% 
          dplyr::filter(sender == sender_oi) %>% 
          dplyr::distinct(ligand) %>% 
          nrow()
        # 如果配体数量大于 1，则设置间隙宽度
        if (num_ligands > 1) {
          sector = rep(width_same_cell_same_ligand_type, times = num_ligands - 1)
        } else {
          sector = c()  # 没有足够的配体时，不设置间隙
        }
        # 每组 sender 之后添加大间隙
        gap = width_different_cell
        return(c(sector, gap))
      }) %>% unlist()
    # 移除最后一个大间隙
    sender_gaps = sender_gaps[-length(sender_gaps)]
    
    # receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    #   sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    #   gap = width_different_cell
    #   return(c(sector,gap))
    # }) %>% unlist()
    # receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    
    receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% 
      lapply(function(receiver_oi) {
        # 计算当前 receiver 对应的受体数量
        num_receptors = circos_links %>% 
          dplyr::filter(receiver == receiver_oi) %>% 
          dplyr::distinct(receptor) %>% 
          nrow()
        # 如果受体数量大于 1，则设置间隙宽度
        if (num_receptors > 1) {
          sector = rep(width_same_cell_same_receptor_type, times = num_receptors - 1)
        } else {
          sector = c()  # 没有足够的受体时，不设置间隙
        }
        # 每组 receiver 之后添加大间隙
        gap = width_different_cell
        return(c(sector, gap))
      }) %>% unlist()
    # 移除最后一个大间隙
    receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    
    #  receiver_gaps = prioritized_tbl_oi$receiver %>% 
    #      unique() %>% 
    #      sort() %>% 
    #      lapply(function(receiver_oi) {
    #        # 强制捕获 receiver_oi
    #        force(receiver_oi)
    #        num_receptors = circos_links %>% 
    #          dplyr::filter(receiver == receiver_oi) %>% 
    #          dplyr::distinct(receptor) %>% 
    #          nrow()
    #        if (num_receptors > 1) {
    #          sector = rep(width_same_cell_same_receptor_type, times = num_receptors - 1)
    #        } else {
    #          sector = c()
    #        }
    #        gap = width_different_cell
    #        return(c(sector, gap))
    #      }) %>% 
    #      unlist()
    # receiver_gaps = receiver_gaps[-length(receiver_gaps)]
    
    print(group_oi)    
    print(length(sender_gaps))
    print(length(receiver_gaps))
    print(length(union(circos_links$ligand, circos_links$receptor)))
    
    gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
    
    if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
      warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
    }
    
    links_circle$weight[links_circle$weight == 0] = 0.01
    circos.clear()
    circos.par(gap.degree = gaps)
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = TRUE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight > 0.01,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.175),
                 grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 reduce = 0,
                 scale = TRUE)
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA) #
    
    title(title)
    p_circos = recordPlot()
    return(p_circos)
    
  })
  names(all_plots) = groups_oi
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  all_plots$legend = p_legend
  
  return(all_plots)
}

unique(multinichenet_output$prioritization_tables$group_prioritization_tbl$receiver)

# 获取所有细胞类型
all_cell_types <- unique(multinichenet_output$prioritization_tables$group_prioritization_tbl$receiver)
immune_cell <- setdiff(all_cell_types, c('HPSC'))
#setdiff(all_cell_types, c('HPSC','pDC','Tdn','cDC1','cDC2'))
# 去掉自身 'Mono_CD14'
setdiff(immune_cell, 'Mono_CD14')

pdf("figures/5_circos/2_Mono_CD14_ATG7_Condition_A1_top20_without_itself.pdf",8,6.5)
prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  groups_oi='A1',
  senders_oi='Mono_CD14',
  receivers_oi = setdiff(immune_cell, 'Mono_CD14'),  # 设置过滤
  top_n = 20, 
  rank_per_group = FALSE
)
nrow(prioritized_tbl_oi_all)
circos_list <- make_circos_group_comparison_2(prioritized_tbl_oi_all, colors_sender, colors_receiver)

prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  groups_oi='C1',
  senders_oi='Mono_CD14',
  receivers_oi = setdiff(immune_cell, 'Mono_CD14'),  # 设置过滤
  top_n = 20, 
  rank_per_group = FALSE
)
nrow(prioritized_tbl_oi_all)
circos_list <- make_circos_group_comparison_2(prioritized_tbl_oi_all, colors_sender, colors_receiver)



prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  groups_oi='A1',
  senders_oi='Mono_CD14_ATG7',
  receivers_oi = setdiff(immune_cell, 'Mono_CD14_ATG7'),  # 设置过滤
  top_n = 20, 
  rank_per_group = FALSE
)
nrow(prioritized_tbl_oi_all)
circos_list <- make_circos_group_comparison_2(prioritized_tbl_oi_all, colors_sender, colors_receiver)

prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  groups_oi='C1',
  senders_oi='Mono_CD14_ATG7',
  receivers_oi = setdiff(immune_cell, 'Mono_CD14_ATG7'),  # 设置过滤
  top_n = 20, 
  rank_per_group = FALSE
)
nrow(prioritized_tbl_oi_all)
circos_list <- make_circos_group_comparison_2(prioritized_tbl_oi_all, colors_sender, colors_receiver)



dev.off()

pdf("figures/5_circos/2_Mono_CD16_ATG7_Condition_C1_top20_without_itself.pdf",8,6.5)

groups <- c('B2','C0', 'C1')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'Mono_CD16_ATG7',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}

groups <- c('B2','C0', 'C1')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'Mono_CD16',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}

dev.off()

pdf("figures/5_circos/2_Platelet_Condition_B1_B0_top20_legend.pdf",8,6)
groups <- c('B0', 'B1')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'Platelet',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}
dev.off()

pdf("figures/5_circos/2_CD8_NELL2_Condition_A1_top20.pdf",8,7.2)
prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  groups_oi='A1',
  senders_oi='CD8_NELL2',
  receivers_oi = setdiff(immune_cell, 'CD8_NELL2'),  # 设置过滤
  top_n = 20, 
  rank_per_group = FALSE
)
nrow(prioritized_tbl_oi_all)
circos_list <- make_circos_group_comparison_2(prioritized_tbl_oi_all, colors_sender, colors_receiver)

dev.off()

pdf("figures/5_circos/2_CD4_Treg_Condition_top20.pdf",8,6.5)
groups <- c('A1', 'B2', 'C0', 'C1')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'CD4_Treg',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}
dev.off()

pdf("figures/5_circos/2_CD4_Temra_Condition_top20.pdf",8,6.5)
groups <- c('A1', 'A2', 'C0', 'C1', 'C2')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'CD4_Temra',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}
dev.off()

pdf("figures/5_circos/2_CD8_Temra_Condition_top20.pdf", 8, 6.5)
# 定义 groups_oi 参数
groups <- c('A1', 'A2', 'C0', 'C1', 'C2')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'CD8_Temra',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}
dev.off()

pdf("figures/5_circos/2_Plasma_Condition_top20.pdf", 8, 6.5)
# 定义 groups_oi 参数
groups <- c('A1', 'B1', 'C1')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'Plasma',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}
dev.off()

pdf("figures/5_circos/2_Bm_CD27_IGHM_SOX5_Condition_top20.pdf", 8, 6.5)
# 定义 groups_oi 参数
groups <- c('A1', 'B1', 'C1')
# 循环调用函数
for (group in groups) {
  circos_list <- process_circos(
    group_name = group,
    sender = 'Bm_CD27_IGHM_SOX5',
    immune_cell = immune_cell,
    top_n = 20,
    colors_sender = colors_sender,
    colors_receiver = colors_receiver
  )
}
dev.off()







# 定义所有的时期
periods <- c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')

# 定义输出 PDF 文件
pdf("figures/5_circos/1_All_Conditions_Top20_with_Legend.pdf", width = 8, height = 6.8)

# 循环遍历每个时期
for (period in periods) {
  
  # 提取当前时期的前 30 个配体-受体对
  prioritized_tbl_oi_all <- get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    groups_oi = period,
    top_n = 20, 
    rank_per_group = FALSE
  )
  
  # 如果有数据才绘图
  if (nrow(prioritized_tbl_oi_all) > 0) {
    
    # 输出当前时期的信息
    message(paste("Generating circos plot for period:", period))
    
    # 生成弦图
    circos_list <- make_circos_group_comparison_2(
      prioritized_tbl_oi_all, 
      colors_sender, 
      colors_receiver
    )
  } else {
    # 如果没有数据，提示并跳过
    message(paste("No data available for period:", period, "- skipping"))
  }
}

# 关闭 PDF 设备
dev.off()

# 打印完成信息
message("All circos plots have been generated and saved to a single PDF file.")

# 定义所有的时期
periods <- c('A0', 'A1', 'A2', 'B0', 'B1', 'B2', 'C0', 'C1', 'C2')

# 定义输出 PDF 文件
pdf("figures/5_circos/1_All_Conditions_Top30_with_Legend.pdf", width = 8, height = 8)

# 循环遍历每个时期
for (period in periods) {
  
  # 提取当前时期的前 30 个配体-受体对
  prioritized_tbl_oi_all <- get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    groups_oi = period,
    top_n = 30, 
    rank_per_group = FALSE
  )
  
  # 如果有数据才绘图
  if (nrow(prioritized_tbl_oi_all) > 0) {
    
    # 输出当前时期的信息
    message(paste("Generating circos plot for period:", period))
    
    # 生成弦图
    circos_list <- make_circos_group_comparison_2(
      prioritized_tbl_oi_all, 
      colors_sender, 
      colors_receiver
    )
  } else {
    # 如果没有数据，提示并跳过
    message(paste("No data available for period:", period, "- skipping"))
  }
}

# 关闭 PDF 设备
dev.off()

# 打印完成信息
message("All circos plots have been generated and saved to a single PDF file.")







# 弦图可视化每个组中，感兴趣的 sender 与 receiver 最优先级的受配体对

# 首先提取数据，获取排序前 top n 的受配体对
prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 30,
  rank_per_group = TRUE
)

# 提取感兴趣的受配体对并合并优先级分数
prioritized_tbl_oi <- multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(prioritized_tbl_oi_all)

# 处理缺失值，如果存在 NA 值将其替换为 0
# prioritized_tbl_oi <- na.omit(prioritized_tbl_oi)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] <- 0

# 获取 senders 和 receivers 的集合
senders_receivers <- union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>%
  sort()

# 设置 senders 和 receivers 的颜色
colors_sender <- RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>%
  magrittr::set_names(senders_receivers)
colors_receiver <- RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>%
  magrittr::set_names(senders_receivers)

# 设置绘图布局
par(mfrow = c(2, 2), xpd = TRUE)

# 调用 make_circos_group_comparison 函数生成每个组的 Circos 图
circos_list <- make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)

library(circlize)
library(ComplexHeatmap)
library(dplyr)

# 获取前 30 个配体-受体对
prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  groups_oi='B1',
  senders_oi='Platelet',
  top_n = 30, 
  rank_per_group = FALSE
)
nrow(prioritized_tbl_oi_all)

# 获取前 30 个配体-受体对
prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  #groups_oi='B1',
  senders_oi='Platelet',
  top_n = 30, 
  rank_per_group = TRUE
)
nrow(prioritized_tbl_oi_all)

plot_circos <- function(prioritized_tbl_oi_all, sender_colors, receiver_colors) {
  
  # 确保数据格式正确
  if (!all(c("group", "sender", "receiver", "ligand", "receptor", "id", "prioritization_score", "prioritization_rank") %in% colnames(prioritized_tbl_oi_all))) {
    stop("Input data must contain columns: 'group', 'sender', 'receiver', 'ligand', 'receptor', 'id', 'prioritization_score', 'prioritization_rank'")
  }
  
  # 去掉数据分组（如果有）
  data <- prioritized_tbl_oi_all %>% ungroup()
  
  # 提取唯一组别
  groups <- unique(data$group)
  
  # 绘图循环
  plots <- lapply(groups, function(group) {
    # 提取特定组的数据
    group_data <- data %>% filter(group == !!group)
    
    # 合并 ligand 和 receptor 的唯一值，用于定义扇区
    sectors <- unique(c(group_data$ligand, group_data$receptor))
    
    # 设置颜色：sender 和 receiver
    ligand_colors <- setNames(sender_colors[group_data$sender], group_data$ligand)
    receptor_colors <- setNames(receiver_colors[group_data$receiver], group_data$receptor)
    grid_colors <- c(ligand_colors, receptor_colors)
    
    # 链接数据准备
    links <- group_data %>% 
      select(ligand, receptor, prioritization_score) %>% 
      rename(from = ligand, to = receptor, value = prioritization_score)
    
    # 归一化分数（用于透明度）
    links$value <- (links$value - min(links$value)) / (max(links$value) - min(links$value) + 1e-6)
    transparency <- 1 - links$value
    
    # 绘制 Circos 图
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
    chordDiagram(
      links,
      order = sectors,
      grid.col = grid_colors,
      transparency = transparency,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.1),
      directional = 1,
      link.arr.length = 0.05,
      link.arr.type = "big.arrow"
    )
    
    # 添加扇区标签
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    }, bg.border = NA)
    
    # 添加标题
    title(group, cex.main = 1.2)
    recordPlot()
  })
  
  names(plots) <- groups
  return(plots)
}

library(RColorBrewer)
prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%  
  filter(id %in% prioritized_tbl_oi_all$id) %>%  
  distinct(id, sender, receiver, ligand, receptor, group) %>%   
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi <- na.omit(prioritized_tbl_oi)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()
# 使用 colorRampPalette 生成足够的颜色
colors_sender <- colorRampPalette(brewer.pal(11, "Spectral"))(length(senders_receivers)) %>%
  magrittr::set_names(senders_receivers)

colors_receiver <- colorRampPalette(brewer.pal(11, "Spectral"))(length(senders_receivers)) %>%
  magrittr::set_names(senders_receivers)
colors_sender
colors_receiver 

plot_circos_legend <- function(prioritized_tbl_oi_all, 
                               sender_colors, receiver_colors, 
                               legend_pdf_path = "legend.pdf",
                               width = 8, height = 6) {
  
  # 确保数据格式正确
  if (!all(c("group", "sender", "receiver", "ligand", "receptor", "id", "prioritization_score", "prioritization_rank") %in% colnames(prioritized_tbl_oi_all))) {
    stop("Input data must contain columns: 'group', 'sender', 'receiver', 'ligand', 'receptor', 'id', 'prioritization_score', 'prioritization_rank'")
  }
  
  # 创建图例
  sender_legend <- Legend(
    at = names(sender_colors),
    legend_gp = gpar(fill = sender_colors),
    title = "Sender Cell Types",
    ncol = 2  # 设置为两列
  )
  
  receiver_legend <- Legend(
    at = names(receiver_colors),
    legend_gp = gpar(fill = receiver_colors),
    title = "Receiver Cell Types",
    ncol = 2  # 设置为两列
  )
  
  # 合并图例
  legend_combined <- packLegend(sender_legend, receiver_legend, direction = "vertical")
  
  # 保存图例为 PDF
  pdf(legend_pdf_path, width = width, height = height)  # 设置PDF文件大小
  ComplexHeatmap::draw(legend_combined, just = "bottom-left")  # 绘制图例
  dev.off()  # 关闭PDF设备
  
  message("Legend saved as ", legend_pdf_path)  # 打印保存路径
}

plot_circos_legend(prioritized_tbl_oi, colors_sender, colors_receiver ,
                   legend_pdf_path ="figures/5_circos/1_Platelet_Condition_top30_legend.pdf",
                   width = 8, height = 6)

pdf("figures/5_circos/1_Platelet_Condition_top30.pdf",10,10)
#plots <- plot_circos(prioritized_tbl_oi, colors_sender, colors_receiver)
plots <- plot_circos_with_arrow(prioritized_tbl_oi, colors_sender, colors_receiver)
dev.off()



plot_circos_with_arrow <- function(prioritized_tbl_oi_all, sender_colors, receiver_colors) {
  # 确保数据格式正确
  if (!all(c("group", "sender", "receiver", "ligand", "receptor", "id", "prioritization_score", "prioritization_rank") %in% colnames(prioritized_tbl_oi_all))) {
    stop("Input data must contain columns: 'group', 'sender', 'receiver', 'ligand', 'receptor', 'id', 'prioritization_score', 'prioritization_rank'")
  } 
  
  data <- prioritized_tbl_oi_all %>% ungroup()  # 去掉数据分组（如果有）
  groups <- unique(data$group) ## 提取唯一组别
  
  # 绘图循环
  plots <- lapply(groups, function(group) {
    # 提取特定组的数据
    group_data <- data %>% filter(group == !!group)
    
    # 合并 ligand 和 receptor 的唯一值，用于定义扇区
    sectors <- unique(c(group_data$ligand, group_data$receptor))
    
    # 设置颜色：sender 和 receiver
    ligand_colors <- setNames(sender_colors[group_data$sender], group_data$ligand)
    receptor_colors <- setNames(receiver_colors[group_data$receiver], group_data$receptor)
    grid_colors <- c(ligand_colors, receptor_colors)
    
    # 链接数据准备
    links <- group_data %>% 
      select(ligand, receptor, prioritization_score) %>% 
      rename(from = ligand, to = receptor, value = prioritization_score)
    
    # 归一化分数（用于透明度）
    links$value <- (links$value - min(links$value)) / (max(links$value) - min(links$value) + 1e-6)
    transparency <- 0.5  # 设置透明度为0.5，确保箭头可见
    
    # 绘制 Circos 图
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
    chordDiagram(
      links,
      order = sectors,
      grid.col = grid_colors,
      transparency = transparency,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.1),
      directional = 1,  # 开启箭头方向
      link.arr.length = 0.1,  # 设置箭头长度，适中的长度
      link.arr.width = 1,     # 设置箭头宽度，适中的宽度
      link.arr.type = "big.arrow"  # 使用大的箭头
    )    
    
    # 添加扇区标签
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    }, bg.border = NA)
    
    # 添加标题
    title(group, cex.main = 1.2)
    p_circos <- recordPlot()
    
    # 返回 Circos 图
    return(p_circos)
  })
  return(plots)
}







plot_circos_with_arrow <- function(prioritized_tbl_oi_all, sender_colors, receiver_colors) {
  # 确保数据格式正确
  if (!all(c("group", "sender", "receiver", "ligand", "receptor", "id", "prioritization_score", "prioritization_rank") %in% colnames(prioritized_tbl_oi_all))) {
    stop("Input data must contain columns: 'group', 'sender', 'receiver', 'ligand', 'receptor', 'id', 'prioritization_score', 'prioritization_rank'")
  } 
  
  data <- prioritized_tbl_oi_all %>% ungroup()  # 去掉数据分组（如果有）
  groups <- unique(data$group)## 提取唯一组别
  
  # 绘图循环
  plots <- lapply(groups, function(group) {
    # 提取特定组的数据
    group_data <- data %>% filter(group == !!group)
    
    # 合并 ligand 和 receptor 的唯一值，用于定义扇区
    sectors <- unique(c(group_data$ligand, group_data$receptor))
    
    # 设置颜色：sender 和 receiver
    ligand_colors <- setNames(sender_colors[group_data$sender], group_data$ligand)
    receptor_colors <- setNames(receiver_colors[group_data$receiver], group_data$receptor)
    grid_colors <- c(ligand_colors, receptor_colors)
    
    # 链接数据准备
    links <- group_data %>% 
      select(ligand, receptor, prioritization_score) %>% 
      rename(from = ligand, to = receptor, value = prioritization_score)
    
    # 归一化分数（用于透明度）
    links$value <- (links$value - min(links$value)) / (max(links$value) - min(links$value) + 1e-6)
    transparency <- 1 - links$value  # 使得值越大透明度越低
    
    # 绘制 Circos 图
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
    chordDiagram(
      links,
      order = sectors,
      grid.col = grid_colors,
      transparency = transparency,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.1),
      directional = 1,  # 开启箭头方向
      link.arr.length = 0.5,  # 增加箭头长度
      link.arr.width = 25,     # 调整箭头宽度
      #link.lwd = 1.25,
      #link.lty = 1,
      # link.border = "gray35",
      link.arr.type = "big.arrow"  # 使用大的箭头
    )    
    # circos.clear()
    # circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
    # chordDiagram(links,
    #              directional = 1,  # 开启箭头
    #              order = sectors,
    #              link.sort = TRUE,
    #              link.decreasing = TRUE,
    #              grid.col = grid_colors,
    #              transparency = transparency,
    #              diffHeight = 0.0075,
    #              direction.type = c("diffHeight", "arrows"),  # 设置箭头显示
    #              annotationTrack = "grid",
    #              preAllocateTracks = list(track.height = 0.175),
    #              grid.border = "gray35",
    #              link.arr.length = 0.1,  # 适当增加箭头长度
    #              link.arr.type = "big.arrow",  # 使用大箭头
    #              link.lwd = 1.25,
    #              link.lty = 1,
    #              link.border = "gray35",
    #              reduce = 0,
    #              scale = TRUE)
    # 添加扇区标签
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    }, bg.border = NA)
    
    # 添加标题
    title(group, cex.main = 1.2)
    p_circos <- recordPlot()
    
    # 返回 Circos 图
    return(p_circos)
  })
  return(plots)
}






