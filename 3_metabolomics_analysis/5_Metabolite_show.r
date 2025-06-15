library(dplyr)
library(tidyr)
library(MetaboDiff)
library(readxl)

metabolites_df <- read.csv("data/0_All_smpdb_metabolites_2024_07_02.csv")
colnames(metabolites_df)
nrow(metabolites_df)
head(metabolites_df,n=1)

colData<- read_excel("data/0_MetaPro-13-21_20211109.xlsx", sheet = "Sample Info")
assay <-  read_excel("data/0_MetaPro-13-21_20211109.xlsx", sheet = "Peak Area Data")  
colnames(assay)

## 去掉Food Component/Plant和Xenobiotics
nrow(assay)
assay <- assay%>% 
        filter(superPathway != 'Xenobiotics')%>% 
        filter(superPathway != 'Food Component/Plant')
nrow(assay)
length(unique(assay$superPathway))
length(unique(assay$subPathway))

# 0. 将 'NA' 字符替换为真正的 NA 值
columns_to_check <- c("Match_name", "PubChem", "kegg", "hmdb")  # 指定列
assay[columns_to_check] <- lapply(assay[columns_to_check], function(x) {
  x[x == "NA"] <- NA  # 替换为 NA
  return(x)
})

# 1. 填充 NA 值
assay$Match_name[is.na(assay$Match_name)] <- assay$name[is.na(assay$Match_name)]
assay$keggId[is.na(assay$keggId)] <- assay$kegg[is.na(assay$keggId)]
assay$hmdbId[is.na(assay$hmdbId)] <- assay$hmdb[is.na(assay$hmdbId)]

# 3. 删除不需要的列
assay <- assay[, !colnames(assay) %in% c("kegg", "hmdb")]

# 检查处理后的数据
head(assay)
colnames(assay)

rowData <-assay[,1:13]
rowData <- as.data.frame(rowData)
rownames(rowData) <- rowData$name
write.csv(rowData, "data/0_Metabbolism_rowData.csv")
head(rowData, n=1)
colData <- as.data.frame(colData)
rownames(colData) <- colData$"SAMPLE ID"
colData$adaptive_group<- colData$"Group ID"
colData$adaptive_group[colData$adaptive_group %in% c('A1','A2','A0','B0')] <- 'before_B1'
colData$adaptive_group[colData$adaptive_group %in% c('B1','B2')] <- 'B1_B2'
# colData$adaptive_group[colData$adaptive_group %in% c('B0')] <- 'B0'
# colData$adaptive_group[colData$adaptive_group %in% c('B1')] <- 'B1'
# colData$adaptive_group[colData$adaptive_group %in% c('A1','A2','A0','B2')] <- NA
colData$innate_group <- colData$"Group ID"
colData$innate_group[colData$innate_group %in% c('A1')] <- 'A1'
colData$innate_group[colData$innate_group =='A0'] <- 'Baseline'
colData$innate_group[colData$innate_group %in% c('B1','B2','B0','A2')] <- NA
write.csv(colData, "data/0_Metabbolism_colData.csv")
head(colData,n=1)

assay <- as.data.frame(assay)
assay <- assay[ , c(2,14:46)]
rownames(assay) <- assay$name
assay$name <- NULL
write.csv(assay, "data/0_Metabbolism_assay_Data.csv")
head(assay,n=1)

nrow(assay)
nrow(colData)
nrow(rowData)

#将三个数据集融合成一个以便于下游分析。
met <- create_mae(assay,rowData,colData)
met

get_SMPDBanno_2 <- function(met,metabolites_df, column_kegg_id, column_hmdb_id, column_cas_id) {
    rowData = rowData(met[["raw"]])    
    # 初始化结果矩阵，包含SMPDB信息
    res = matrix(NA, nrow = nrow(rowData), ncol = ncol(metabolites_df))
    colnames(res) = paste0("SMPDB_", colnames(metabolites_df))

    # 精确匹配的函数，优先匹配最相似的项
    find_best_match <- function(row_value, db_column) {
        possible_matches = which(db_column == row_value)
        if (length(possible_matches) == 1) {
            return(metabolites_df[possible_matches, , drop = FALSE])  # 使用 drop = FALSE 保持维度
        } else if (length(possible_matches) > 1) {
            # 当有多个匹配项时，选择第一个匹配项（或定义更复杂的选择逻辑）
            best_match = possible_matches[1]
            return(metabolites_df[best_match, , drop = FALSE])
        } else {
            return(matrix(NA, nrow = 1, ncol = ncol(metabolites_df)))
        }
    }

    # 逐步匹配 KEGG, HMDB, CAS
    if (!is.na(column_kegg_id)) {
        for (i in 1:nrow(rowData)) {
            match_result = find_best_match(rowData[i, column_kegg_id], metabolites_df$KEGG.ID)
            if (!all(is.na(match_result))) {
                res[i, ] = as.matrix(match_result)  # 确保 match_result 是一个矩阵
            }
        }
    }
    
    if (!is.na(column_hmdb_id)) {
        for (i in 1:nrow(rowData)) {
            if (all(is.na(res[i, ]))) {  # 只在没有匹配到时再继续匹配
                match_result = find_best_match(rowData[i, column_hmdb_id], metabolites_df$HMDB.ID)
                if (!all(is.na(match_result))) {
                    res[i, ] = as.matrix(match_result)  # 确保 match_result 是一个矩阵
                }
            }
        }
    }

    if (!is.na(column_cas_id)) {
        for (i in 1:nrow(rowData)) {
            if (all(is.na(res[i, ]))) {  # 只在没有匹配到时再继续匹配
                match_result = find_best_match(rowData[i, column_cas_id], metabolites_df$CAS)
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


rowData(met[["raw"]])[1, 11]
rowData(met[["raw"]])[1, 12]
rowData(met[["raw"]])[1, 10]

## 精准匹配：逐步匹配 KEGG, HMDB, CAS
met <- get_SMPDBanno_2(met,metabolites_df,
                           column_kegg_id=11,
                           column_hmdb_id=12,
                          column_cas_id= 10 )


SMPDB_res <-  rowData(met[["raw"]])
nrow(SMPDB_res)
write.csv(SMPDB_res, "data/1_SMPDB_res_KEGG_HMDB_CAS.csv")

head(SMPDB_res$SMPDB_SMPDB.ID)
# 计算非空值的数量
sum(!is.na(SMPDB_res$SMPDB_SMPDB.ID))
colnames(SMPDB_res)

length(unique(metabolites_df$SMPDB.ID))

# 数据验证
nrow(metabolites_df)
sum(!is.na(metabolites_df$CAS))
sum(!is.na(metabolites_df$HMDB.ID))
length(unique(metabolites_df$CAS))
length(unique(metabolites_df$HMDB.ID))
length(unique(metabolites_df$Pathway.Name))
nrow(rowData)
length(unique(rowData$casId))
length(unique(rowData$hmdbId))
sum(unique(rowData$casId) %in% unique(metabolites_df$CAS))
sum(unique(rowData$keggId) %in% unique(metabolites_df$KEGG.ID))
sum(unique(rowData$hmdbId) %in% unique(metabolites_df$HMDB.ID))
head(metabolites_df,n=1)
colnames(metabolites_df)

library(dplyr)
library(stringdist)

# 将匹配的 Metabolite.Name 和对应的 distance 添加到 matched_data 中，并合并 metabolites_df 中的其他列
match_metabolites <- function(rowData, metabolites_df, method,distance_threshold = 1) {
  
  # 创建一个空的结果数据框
  matched_data <- rowData
  
  # 遍历 rowData 中的每一行
  for (i in 1:nrow(rowData)) {
    #row_name <- rowData$name[i]
    row_name <- rowData$Match_name[i]
      
    # 计算当前 row_name 和 metabolites_df 中每个 Metabolite.Name 的距离
    distances <- stringdist::stringdist(row_name, metabolites_df$Metabolite.Name, method = method)
    
    # 找到最小的距离并确定是否满足匹配条件
    min_distance <- min(distances)
    
    if (min_distance <= distance_threshold) {
      # 找到最匹配的 Metabolite.Name
      best_match_idx <- which.min(distances)
      best_match_name <- metabolites_df$Metabolite.Name[best_match_idx]
      
      # 将匹配的 Metabolite.Name 和对应的 distance 添加到 matched_data 中
      matched_data$Matched_Metabolite[i] <- best_match_name
      matched_data$Distance[i] <- min_distance
      
      # 合并 metabolites_df 中与 Metabolite.Name 对应的其他列
      matched_data[i, paste0("SMPDB_", colnames(metabolites_df))] <- metabolites_df[best_match_idx, ]
    } else {
      # 如果没有匹配，保持 NA
      matched_data$Matched_Metabolite[i] <- NA
      matched_data$Distance[i] <- NA
      # 如果没有匹配，metabolites_df 列也填充 NA
      matched_data[i, paste0("SMPDB_", colnames(metabolites_df))] <- NA
    }
  }
  
  return(matched_data)
}

# 1.	“jw” (Jaro-Winkler Distance)
# 	•	该方法适用于代谢物名称中有一定的拼写差异、相似的前缀部分以及格式差异。特别是在代谢物名称中常见的拼写变体和前缀匹配时，jw 方法可以处理得很好。
# 	2.	“osa” (Optimal String Alignment Distance)
# 	•	如果代谢物名称之间的差异主要表现为字符交换、插入或删除操作，那么 osa 方法可能是一个不错的选择。
# 	3.	“lv” (Levenshtein Distance)
# 	•	如果代谢物名称差异较小，主要是拼写错误或略微的字符错位，那么 lv 方法非常有效。


# 假设 rowData 和 metabolites_df 已经加载，并包含相应的列
result_data <- match_metabolites(rowData, metabolites_df,"lv", distance_threshold = 5)

write.csv(result_data,'data/1_match_metabolites_result_data_lv.csv')
# 假设 rowData 和 metabolites_df 已经加载，并包含相应的列
result_data <- match_metabolites(rowData, metabolites_df,method="osa", distance_threshold = 5)

write.csv(result_data,'data/1_match_metabolites_result_data_osa.csv')

# 假设 rowData 和 metabolites_df 已经加载，并包含相应的列
result_data_jaccard <- match_metabolites(rowData, metabolites_df,method="jaccard", distance_threshold = 0.5)

write.csv(result_data_jaccard,'data/1_match_metabolites_result_data_jaccard.csv')

# 假设 rowData 和 metabolites_df 已经加载，并包含相应的列
result_data_jw <- match_metabolites(rowData, metabolites_df,method="jw", distance_threshold = 0.5)

write.csv(result_data_jw,'data/1_match_metabolites_result_data_jw.csv')

# 假设 rowData 和 metabolites_df 已经加载，并包含相应的列
result_data_cos <- match_metabolites(rowData, metabolites_df,method="cosine", distance_threshold = 0.5)

write.csv(result_data_cos,'data/1_match_metabolites_result_data_cosine.csv')

SMPDB_res$Matched_Metabolite <- SMPDB_res$SMPDB_Metabolite.Name
SMPDB_res$Distance <- 0
# 确保 SMPDB_res 的列名一致并调整顺序
SMPDB_res <- SMPDB_res[, colnames(result_data_jw)]
SMPDB_res <- as.data.frame(SMPDB_res)

class(SMPDB_res)
class(result_data_jw)

# 检查列类型是否一致
str(result_data_jw$SMPDB_ChEBI.ID)
str(result_data_cos$SMPDB_ChEBI.ID)
str(result_data_jaccard$SMPDB_ChEBI.ID)
str(SMPDB_res$SMPDB_ChEBI.ID)

# 统一列类型为字符类型
result_data_jw$SMPDB_ChEBI.ID <- as.character(result_data_jw$SMPDB_ChEBI.ID)
result_data_cos$SMPDB_ChEBI.ID <- as.character(result_data_cos$SMPDB_ChEBI.ID)
result_data_jaccard$SMPDB_ChEBI.ID <- as.character(result_data_jaccard$SMPDB_ChEBI.ID)
SMPDB_res$SMPDB_ChEBI.ID <- as.character(SMPDB_res$SMPDB_ChEBI.ID)

result_data_jw$Method <- 'JW'
result_data_cos$Method <- 'Cosine' 
result_data_jaccard$Method <- 'Jaccard'
SMPDB_res$Method <- 'KEGG_HMDB_CAS'

colnames(result_data_jw)
colnames(result_data_cos)
colnames(result_data_jaccard)
colnames(SMPDB_res)

# 合并四个数据集
merged_data <- bind_rows(result_data_jw, result_data_cos, result_data_jaccard, SMPDB_res)

# 按照 `name` 列排序
merged_data <- merged_data %>%
  arrange(name)

# 检查结果
head(merged_data)

write.csv(merged_data,'data/1_Match_4_method_res_Merged.csv')

length(unique(merged_data$SMPDB_Pathway.Name))
colnames(merged_data)

length(unique(merged_data$SMPDB_Pathway.Name))
length(unique(merged_data$SMPDB_Pathway.Subject))
length(unique(merged_data$subPathway))
length(unique(merged_data$superPathway))

unique(merged_data$SMPDB_Pathway.Subject)
head(unique(merged_data$SMPDB_Pathway.Name))

All_Condition_pathway <- read.csv('data/2_All_ct3_Condition_pathway_activity.csv', check.names = FALSE)
sum(colnames(All_Condition_pathway) %in% unique(merged_data$SMPDB_Pathway.Name))
# 查看列名
length(colnames(All_Condition_pathway))
head(colnames(All_Condition_pathway))

# 加载所需包
library(dplyr)
library(stringdist)
library(purrr)

# 定义匹配函数，找到最相似的名称及其距离
find_best_match <- function(pathway_name, pathway_list) {
  if (is.na(pathway_name) || pathway_name == "") {
    return(list(Matched_Pathway = NA, Distance = NA))
  }
  distances <- stringdist(pathway_name, pathway_list, method = "jw") # 计算 Jaro-Winkler 距离
  best_match_index <- which.min(distances) # 找到距离最小的索引
  list(
    Matched_Pathway = pathway_list[best_match_index],
    Distance = distances[best_match_index]
  )
}



new_pathway_name <-  colnames(All_Condition_pathway)

# 对非空的 SMPDB_Pathway.Name 进行匹配
matched_results <- merged_data$SMPDB_Pathway.Name %>%
  map(~ find_best_match(.x, new_pathway_name)) %>%
  transpose() # 将结果转置为列

# 将匹配结果添加到 merged_data
merged_data <- merged_data %>%
  mutate(
    Pathway_Activity = unlist(matched_results$Matched_Pathway),
    Activity_Distance = unlist(matched_results$Distance)
  )

# 查看结果
head(merged_data$Pathway_Activity)
head(merged_data$SMPDB_Pathway.Name)

#Compass数据
Compass_df <- read.csv('data/2_cohens_combined_data.csv')
colnames(Compass_df)

nrow(Compass_df)
sum(unique(Compass_df$subsystem) %in% unique(merged_data$SMPDB_Pathway.Name))
sum(unique(Compass_df$subsystem) %in% colnames(All_Condition_pathway))
length(unique(Compass_df$subsystem))
head(Compass_df$metadata_r_id)
head(unique(Compass_df$subsystem))

new_pathway_name <-  unique(Compass_df$subsystem)

# 对非空的 SMPDB_Pathway.Name 进行匹配
matched_results <- merged_data$SMPDB_Pathway.Name %>%
  map(~ find_best_match(.x, new_pathway_name)) %>%
  transpose() # 将结果转置为列

# 将匹配结果添加到 merged_data
merged_data <- merged_data %>%
  mutate(
    Pathway_Compass = unlist(matched_results$Matched_Pathway),
    Compass_Distance = unlist(matched_results$Distance)
  )

# 查看结果
head(merged_data$Pathway_Compass)
head(merged_data$SMPDB_Pathway.Name)

write.csv(merged_data,'data/2_Match_4_method_res_Merged_with_info.csv')

colnames(merged_data)

# 选择需要的列
selected_data <- merged_data %>%
  select('name','SMPDB_Metabolite.Name',
         'SMPDB_Pathway.Name',
         'Pathway_Activity', 
         'Pathway_Compass', 
         'subPathway', 'superPathway','Method',
         'Distance','Compass_Distance','Activity_Distance')
write.csv(selected_data,'data/2_selected_Merged_with_info.csv')



#https://www.jianshu.com/p/80af83c2f630
#剔除缺失值，计算代谢物的相对丰度。
met = knn_impute(met,cutoff=0.25)
met <- normalize_met(met)

#对单个代谢物进行差异分析，主要用T检验和ANOVA分析。
met = diff_test(met,
                 group_factors = c("adaptive_group","innate_group","Group ID")) #,"innate_group","Group ID"
 str(metadata(met), max.level=2)


length(rownames(rowData))
length(unique(rownames(rowData)))
length(metadata(met)$ttest_adaptive_group_before_B1_vs_B1_B2$metabolite
)

# 获取代谢物列表并排序
metabolites_of_interest <- sort(metadata(met)$ttest_adaptive_group_before_B1_vs_B1_B2$metabolite)

length(metabolites_of_interest)

# 获取代谢物在 norm_imputed 数据中的行号并排序
row_indices <- which(rowData(met[["norm_imputed"]])$name %in% metabolites_of_interest)

# 对行号进行排序，确保它们按照 metabolites_of_interest 排序
row_indices_sorted <- row_indices[order(match(rowData(met[["norm_imputed"]])$name[row_indices], metabolites_of_interest))]

# 查看排序后的行号
head(row_indices_sorted)

head(metabolites_of_interest)

plot_metabolic_barplots(met, metabolites_of_interest,
                        output_pdf = 'figures/1_All_metabolic_barplot_0.05.pdf', 
                                    groups = c("A0", "A1", "A2", "B0", "B1", "B2"),
                                    colors = c("grey", "red","#C1E6F3", "#C1E6F3", "orange", "#C1E6F3"),
                                    plots_per_page = 20)

##精氨酸
# 获取代谢物名称
metabolites_of_interest <- c("2-oxoarginine*","argininate*","arginine",
                             "dimethylarginine (SDMA + ADMA)",
                             "homoarginine","N-acetylarginine")
plot_metabolic_barplots(met, metabolites_of_interest,
                        output_pdf = 'figures/1_metabolic_barplot_arginine_0.05.pdf', 
                                    groups = c("A0", "A1", "A2", "B0", "B1", "B2"),
                                    # colors = c("grey", "red","lightblue", "lightblue", 
                                    #            "orange", "lightblue"),
                        colors = c("grey", "red","#C1E6F3", "#C1E6F3", "orange", "#C1E6F3"),
                                    plots_per_page = 20)

#adaptive
# adaptive_updated[adaptive_updated$pval < 0.05 & adaptive_updated$dm > 0.25, ]
metabolites_of_interest <- c('2-aminoheptanoate','2-hydroxyheptanoate*',
                             '3-amino-2-piperidone','4-hydroxyphenylacetate',
                             'glucuronate','tryptophan betaine')
plot_metabolic_barplots(met, metabolites_of_interest,
                        output_pdf = 'figures/1_metabolic_barplot_Adaptive_0.05.pdf', 
                                    groups = c("A0", "A1", "A2", "B0", "B1", "B2"),
                                    colors = c("grey", "red","#C1E6F3", "#C1E6F3", "orange", "#C1E6F3"),
                                    plots_per_page = 20)

#innate
#innate[innate$pval < 0.05 & innate$dm > 0.5, ],12
metabolites_of_interest <- c('1-carboxyethylisoleucine',
                             '1-palmitoleoyl-2-linoleoyl-GPC (16:1/18:2)*',
                             '2\'-deoxyuridine',
                             'androstenediol (3beta,17beta) monosulfate (2)',
                             'galactitol (dulcitol)',
                             'gamma-glutamylalanine','gamma-glutamylserine',
                             'N-acetylaspartate (NAA)',
                             'N-palmitoylserine','retinol (Vitamin A)',
                             'sarcosine','tyramine O-sulfate')
plot_metabolic_barplots(met, metabolites_of_interest,
                        output_pdf = 'figures/1_metabolic_barplot_Innate_0.05.pdf', 
                                    groups = c("A0", "A1", "A2", "B0", "B1", "B2"),
                                    colors = c("grey", "red","#C1E6F3", "#C1E6F3", "orange", "#C1E6F3"),
                                    plots_per_page = 20)





find_matching_elements <- function(data, pattern) {
  # 确保列名存在
  required_columns <- c("Pathway_Activity", "Pathway_Compass", "subPathway", "SMPDB_Pathway.Name")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("缺失以下列名:", paste(missing_columns, collapse = ", ")))
  }
  
  # 查找每列中包含 pattern 的唯一值
  result <- list(
    Pathway_Activity = unique(data$Pathway_Activity[grep(pattern, data$Pathway_Activity, ignore.case = TRUE)]),
    Pathway_Compass = unique(data$Pathway_Compass[grep(pattern, data$Pathway_Compass, ignore.case = TRUE)]),
    subPathway = unique(data$subPathway[grep(pattern, data$subPathway, ignore.case = TRUE)]),
    SMPDB_Pathway.Name = unique(data$SMPDB_Pathway.Name[grep(pattern, data$SMPDB_Pathway.Name, ignore.case = TRUE)])
  )
  
  # 返回结果
  return(result)
}

# unique(selected_data$Pathway_Activity)
# unique(selected_data$Pathway_Compass)
# unique(selected_data$subPathway)
# unique(selected_data$SMPDB_Pathway.Name)

plot_metabolic_barplots <- function(met, metabolites_of_interest,output_pdf = "output.pdf", 
                                    groups = c("A0", "A1", "A2", "B0", "B1", "B2"),
                                    colors = c("grey", "red","#C1E6F3", "#C1E6F3", "orange", "#C1E6F3"),
                                    plots_per_page = 20) {
  # 获取代谢物列表并排序
  metabolites_of_interest <- sort(metabolites_of_interest)
  print(unlist(metabolites_of_interest))
  # 获取代谢物在 norm_imputed 数据中的行号并排序
  row_indices <- which(rowData(met[["norm_imputed"]])$name %in% metabolites_of_interest)
  row_indices_sorted <- row_indices[order(match(rowData(met[["norm_imputed"]])$name[row_indices], metabolites_of_interest))]
  
  # 检查代谢物数量
  if (length(row_indices_sorted) == 0) {
    stop("未找到符合条件的代谢物。")
  }
  
  # 确保组因子是有序的因子
  ids <- met$'Group ID' %in% groups
  group_factor <- factor(met$'Group ID'[ids], levels = groups)
  
  # 打开 PDF 输出
  pdf(output_pdf, width = 12.5, height = 10)
  
  # 计算页数
  num_pages <- ceiling(length(row_indices_sorted) / plots_per_page)
  
  # 绘图
  for (page in 1:num_pages) {
    # 设置当前页的绘图窗口
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(row_indices_sorted))
    
    # 更新绘图参数：每页最多显示 20 个图
    par(mfrow = c(4, 5)) # 每页 4 行 5 列
    
    for (i in start_idx:end_idx) {
      idx <- row_indices_sorted[i] # 获取代谢物的行号
      y_data <- assay(met[["norm_imputed"]])[idx, ids] # 获取数据
      
      # 绘制箱线图，不显示异常值
      boxplot(y_data ~ group_factor, 
              main = rowData(met[["norm_imputed"]])$name[idx],  # 代谢物名称作为标题
              xlab = "", ylab = "Normalized values", 
              frame = FALSE, col = colors, 
              outline = FALSE)
      
      # 计算新的中位数（不包括异常值）
      median_values <- tapply(y_data, group_factor, function(x) median(x, na.rm = TRUE))
      
      # 重新绘制中位线
      for (j in seq_along(median_values)) {
        segments(j - 0.4, median_values[j], j + 0.4, median_values[j], lwd = 2, col = "black")
      }
      
      # ANOVA 分析和 Tukey's HSD 检验
      anova_result <- aov(y_data ~ group_factor)                   
      tukey_result <- TukeyHSD(anova_result)
      
      # 获取显著性标记
      significant_comparisons <- tukey_result$group_factor[, "p adj"] < 0.05
      if (any(significant_comparisons)) {
        y_max <- max(y_data[!is.na(y_data) & y_data != Inf], na.rm = TRUE)
        y_offset <- 0.05 * (y_max - min(y_data, na.rm = TRUE))
        text_y_pos <- y_max
        
        for (comp_idx in which(significant_comparisons)) {
          p_value <- tukey_result$group_factor[comp_idx, "p adj"]
          sign_marker <- ifelse(p_value < 0.001, "***",
                                ifelse(p_value < 0.01, "**",
                                       ifelse(p_value < 0.05, "*", "")))
          if (sign_marker != "") {
            text(x = 3, y = text_y_pos, labels = sign_marker, cex = 0.8, col = "red")
            text_y_pos <- text_y_pos + y_offset # 更新文本位置避免重叠
          }
        }
      }
    }
  }
  
  # 关闭 PDF 输出
  dev.off()
}



find_top3_per_column <- function(pathway_data, match_result_column,top_n=3) {
  # 确保匹配元素存在于列名中
  matching_columns <- colnames(pathway_data) %in% match_result_column
  
  if (!any(matching_columns)) {
    stop("No matching columns found.")
  }
  
  # 创建结果列表存储每列的 top 3
  results <- list()
  
  for (col in colnames(pathway_data)[matching_columns]) {
    # 获取列数据
    col_data <- pathway_data[, col]
    
    # 获取对应的行名
    row_names <- rownames(pathway_data)
    
    # 按值降序排序并取前n
    top_n_indices <- order(col_data, decreasing = TRUE)[1:top_n]
    top_n_values <- col_data[top_n_indices]
    top_n_row_names <- row_names[top_n_indices]
    
    # 存储结果
    results[[col]] <- data.frame(
      Row_Name = top_n_row_names,
      Value = top_n_values,
      stringsAsFactors = FALSE
    )
  }
  
  return(results)
}

plot_Compass_heatmap_function <- function(combined_data, condition_columns, 
                                  output_pdf_path ) {
  library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(reshape2)
  # 初始化 plot_list 用于存储每个细胞类型的图
  plot_list <- list()
  
  # 颜色设置
  color.2 <- colorRampPalette(c('#24A1FA','#24A1FA','#FCFBFD', '#E83E15','#E20612'))(1000)
  
  # 循环生成每个细胞类型的图并存入 plot_list
  for (cell in sort(unique(combined_data$cell_type))) {
    
    cell_data <- combined_data %>%
      filter(cell_type == cell)
    
    # 转换数据格式，便于计算均值
    plot_data <- cell_data %>%
      pivot_longer(cols = all_of(condition_columns), names_to = "Condition", values_to = "Value")
    
    # 计算各个 Condition 列的均值
    condition_means <- plot_data %>%
      group_by(Condition) %>%
      summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
    
    # 计算每个 metadata_r_id 的 total_deviation
    cell_data <- cell_data %>%
      rowwise() %>%
      mutate(
        total_dev = sum(abs(c_across(all_of(condition_columns)) - unlist(condition_means[match(condition_columns, condition_means$Condition), "mean_value"])), na.rm = TRUE)
      ) %>%
      ungroup()
    
    max_total_dev <- max(cell_data$total_dev, na.rm = TRUE)  # 计算最高值
    filtered_data <- cell_data %>%
      filter(total_dev < max_total_dev)  # 筛选出 total_dev 小于最高值的数据
    
    # 转换数据格式，便于绘图
    plot_data <- filtered_data %>%
      pivot_longer(cols = all_of(condition_columns), names_to = "Condition", values_to = "Value") %>%
      left_join(condition_means, by = "Condition")  # 加入均值 # total_dev 范围
    total_dev_range <- range(filtered_data$total_dev, na.rm = TRUE)
    
    # 创建 ggplot 对象
    p <- ggplot(plot_data, aes(x = Condition, y = Value, group = metadata_r_id, color = total_dev)) +
      geom_point(size = 0.3) +
      geom_line(size = 0.1) +
      scale_color_gradientn(colors = color.2, values = scales::rescale(total_dev_range)) +
      geom_line(aes(y = mean_value), color = "red", size = 0.5, linetype = "dashed") +
      labs(title = paste(cell),
           x = "Condition",
           y = "Cohen's d",
           color = "total_dev") +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
      ) +
      guides(color = guide_colorbar(barwidth = 0.5, barheight = 4))
    
    # 将图表添加到列表
    plot_list[[cell]] <- p
  }
  
  # 生成热图 PDF 合并图
  pdf(output_pdf_path, width = 25, height = 21)
  combined_plot <- wrap_plots(plot_list, ncol = 5)  # 选择合适的列数
  print(combined_plot)
  dev.off()
  
  # 返回 plot_list 如果需要进一步处理
  return(plot_list)
}

plot_Activity_pathway_heatmap_only_1 <- function(pathway_ct3_merged,Pathway_Activity,margin_size) {
    # 定义处理和绘制热图的函数
    library(tidyr)
    library(dplyr)
    library(ggplot2)
  # 数据处理部分
  # 提取 '.' 后部分为 Condition，并将宽表转换成长表
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
    
  pathway_ct3_long <-  pathway_ct3_long%>% filter(Pathway %in% Pathway_Activity)
  # 使用 ggplot 绘制热图
  ggplot(pathway_ct3_long, aes(x = Condition, y = cell_type, fill = Value)) +
    geom_tile() +  # 使用 geom_tile 绘制热图
    facet_wrap(~ Pathway, scales = "free") +  # 按 Pathway 分面
    scale_fill_gradient(low = "white", high = "red") +  # 颜色渐变
    theme_minimal() + 
    labs(title = "Pathway Activity Heatmap", 
         x = "Condition", 
         y = "Cell Type", 
         fill = "Activity Value") + 
    theme(#axis.text.x = element_text(hjust = 0),  # 调整 X 轴标签的角度,#angle = 45, 
          strip.text = element_text(size = 10),# 调整分面标题的字体大小
       plot.margin = unit(margin_size, "char"))  #margin_size=(上，右，下，左）
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

# 定义一个完整的函数，封装所有步骤
Process_Find_Pathway_analysis <- function(
  merged_data, pathway_ct1_merged, pathway_ct2_merged, pathway_ct3_merged,
  Compass_df, metabolites_df, met, pattern,save_name) {
  library(cowplot)
  # Helper function: 检查数据是否为空
  check_empty <- function(data) {
    if (is.null(data) || length(data) == 0 || all(is.na(data))) {
      return(TRUE)
    }
    return(FALSE)
  }
  
  # Step 1: 查找匹配的通路
  if (!check_empty(merged_data)) {
    match_result <- find_matching_elements(merged_data, pattern)
    print(match_result)
    print('Pathway_Activity:')
    if (!check_empty(pathway_ct1_merged)) {
      print(colnames(pathway_ct1_merged)[grep(pattern, colnames(pathway_ct1_merged), ignore.case = TRUE)])
    }
    
    print('Pathway_Compass:')
    if (!check_empty(Compass_df)) {
      print(unique(Compass_df$subsystem[grep(pattern, Compass_df$subsystem, ignore.case = TRUE)]))
    }
    
    print('SMPDB_Pathway.Name:')
    if (!check_empty(metabolites_df)) {
      print(length(unique(metabolites_df$Pathway.Name[grep(pattern, metabolites_df$Pathway.Name, ignore.case = TRUE)])))
    }
  } else {
    match_result <- NULL
  }
  
  # Step 2: 筛选数据并生成柱状图
  if (!check_empty(merged_data) && !check_empty(match_result)) {
    subpathway_names <- unique(merged_data$name[merged_data$subPathway %in% match_result$subPathway])
    smpdb_names <- unique(merged_data$name[merged_data$SMPDB_Pathway.Name %in% match_result$SMPDB_Pathway.Name & 
                                             merged_data$Method == "KEGG_HMDB_CAS"])
    
    if (!check_empty(subpathway_names)) {
      plot_metabolic_barplots(met, subpathway_names, output_pdf = paste0("figures/3_",save_name,"_Subpathway_barplot_0.05.pdf"))
    }
    if (!check_empty(smpdb_names)) {
      plot_metabolic_barplots(met, smpdb_names, output_pdf = paste0("figures/3_",save_name,"_SMPDB_barplot_0.05.pdf"))
    }
  }
  
  # Step 3: 绘制通路活动的热图
   Pathway_activity  <- colnames(pathway_ct1_merged)[grep(pattern, colnames(pathway_ct1_merged), ignore.case = TRUE)]
  if (!check_empty(pathway_ct1_merged) && !check_empty(Pathway_activity)) {
    n_col <- length(Pathway_activity)
    pdf(paste0("figures/3_",save_name,"_pathway_activity_heatmap.pdf"), width = 5 * n_col, height = 6)
      pathway_heatmaps <- plot_Activity_pathway_heatmap(pathway_ct1_merged, Pathway_activity, margin_size = c(10, 3.5, 10, 3.5))
      print(plot_grid(plotlist = pathway_heatmaps, ncol = n_col))
      pathway_heatmaps <- plot_Activity_pathway_heatmap(pathway_ct2_merged, Pathway_activity, margin_size = c(9, 3.5, 9, 3.5))
      print(plot_grid(plotlist = pathway_heatmaps, ncol = n_col))
      pathway_heatmaps <- plot_Activity_pathway_heatmap(pathway_ct3_merged, Pathway_activity, margin_size = c(0.5, 0.5, 0.5, 0.5))
      print(plot_grid(plotlist = pathway_heatmaps, ncol = n_col))
    dev.off()
      print('绘制通路活动的热图 Pathway_activity_heatmap:Done!')
  }
  
  # Step 4: 绘制Compass趋势图
  if (!check_empty(Compass_df) && !check_empty(match_result$Pathway_Compass)) {
    plot_data <- Compass_df %>% filter(subsystem %in% match_result$Pathway_Compass)
      plot_results <- plot_Compass_heatmap_function(
        combined_data = plot_data, 
        condition_columns = c("A0", "A1", "A2", "B0", "B1", "B2", "C0", "C1", "C2"),
        output_pdf_path = paste0("figures/3_",save_name,"_Cohens_d_plots.pdf")
      )
    print('绘制Compass趋势图:Done!')
  }
}

pathway_ct1_merged <- read.csv('data/3_All_ct1_merged_pathway_activity.csv',check.names = FALSE,row.names = 1)
pathway_ct2_merged <- read.csv('data/3_All_ct2_merged_pathway_activity.csv',check.names = FALSE,row.names = 1)
pathway_ct3_merged <- read.csv('data/3_All_ct3_merged_pathway_activity.csv',check.names = FALSE,row.names = 1)
head(pathway_ct1_merged)

pathway_ct1 <- read.csv('data/3_All_ct1_celltypel1_pathway_activity.csv',check.names = FALSE,row.names = 1)
pathway_ct2 <- read.csv('data/3_All_ct2_celltypel2_pathway_activity.csv',check.names = FALSE,row.names = 1)
pathway_ct3 <- read.csv('data/3_All_ct3_celltypel3_pathway_activity.csv',check.names = FALSE,row.names = 1)
head(colnames(pathway_ct1))
rownames(pathway_ct1)
rownames(pathway_ct2)
rownames(pathway_ct3)

#Compass数据
Compass_df <- read.csv('data/2_cohens_combined_data.csv')
colnames(Compass_df)

## Sphingomyelins
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Bile",
  save_name= "Bile"
)

## Sphingomyelins
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Sphing",
  save_name= "Sphingomyelins"
)

Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Steroid",
  save_name= "Steroid"
)

Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Androgenic Steroids",
  save_name= "Androgenic Steroids"
)

Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Lysophospholipid",
  save_name= "Lysophospholipid"
)

Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Phosphatidylcholine",
  save_name= "Phosphatidylcholine"
)

Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Amino Acid",
  save_name= "Amino_Acid"
)

Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Tryptophan",
  save_name= "Tryptophan"
)



# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Bile",
  save_name= "Bile"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Arginine",
  save_name= "Arginine"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Ascorbate",
  save_name= "Ascorbate"
)

# Primary bile acid biosynthesis
# Drug metabolism − cytochrome P450
# Pentose and glucuronate interconversions
# Metabolism of xenobiotics by cytochrome P450
# Lipoic acid metabolism

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Lipoic",
  save_name= "Lipoic"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Glucuronate",
  save_name= "Glucuronate"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Cytochrome",
  save_name= "Cytochrome"
)





# Mucin type O−glycan biosynthesis
#Phosphonate and phosphinate metabolism
#Inositol phosphate metabolism
#Lysine degradation
#Glycosaminoglycan biosynthesis − chondroitin sulfate / dermatan sulfate 
#Starch and sucrose metabolism

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Mucin",
  save_name= "Mucin"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Starch",
  save_name= "Starch"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Sphingolipid",
  save_name= "Sphingolipid"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Glycosaminoglycan",
  save_name= "Glycosaminoglycan"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Glutamate",
  save_name= "Glutamate"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Phosphate",
  save_name= "Phosphate"
)

# 使用函数
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Lysine",
  save_name= "Lysine"
)





# CD4_Treg #Other types of O−glycan biosynthesis
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Glycan",
  save_name= "Glycan"
)



# CD8_Temra
#Synthesis and degradation of ketone bodies
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "ketone",
  save_name= "Ketone"
)

# CD8_NELL2
#Thiamine metabolism
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Thiamine",
  save_name= "Thiamine"
)

# CD8_NELL2
# Biotin metabolism
Process_Find_Pathway_analysis(
  merged_data = merged_data, 
  pathway_ct1_merged = pathway_ct1_merged, 
  pathway_ct2_merged = pathway_ct2_merged, 
  pathway_ct3_merged = pathway_ct3_merged, 
  Compass_df = Compass_df, 
  metabolites_df = metabolites_df, 
  met = met, 
  pattern = "Biotin",
  save_name= "Biotin"
)



pathways <- colnames(pathway_ct1_merged)[1:85]
length(pathways)
pathways

library(cowplot)  # 用于拼接多图
# 调用函数绘制Activity_pathway热图
pdf("figures/4_ALL_pathway_activity_heatmap.pdf", width = 40, height = 66)
pathway_heatmaps <-plot_Activity_pathway_heatmap(pathway_ct1_merged,pathways,margin_size=c(10,3.5,10,3.5))
plot_grid(plotlist = pathway_heatmaps, ncol = 8) 
pathway_heatmaps <-plot_Activity_pathway_heatmap(pathway_ct2_merged,pathways,margin_size=c(9,3.5,9,3.5))
plot_grid(plotlist = pathway_heatmaps, ncol = 8)  
pathway_heatmaps <-plot_Activity_pathway_heatmap(pathway_ct3_merged,pathways,margin_size=c(0.5,0.5,0.5,0.5))
plot_grid(plotlist = pathway_heatmaps, ncol = 8)  
dev.off()









nrow(rowData)
colnames(rowData)

head(rowData)

table(rowData$superPathway)
table(rowData$subPathway)

# 加载必要的包
library(ggplot2)
library(dplyr)
pdf("figures/5_Superpathway_num_barplot.pdf", width = 3, height = 1.35)
# 计算 superPathway 的频数
super_pathway_table <- as.data.frame(table(rowData$superPathway)) %>%
  arrange(desc(Freq))
colnames(super_pathway_table) <- c("SuperPathway", "Frequency")

# 绘制条形图
ggplot(super_pathway_table, aes(x = reorder(SuperPathway, -Frequency),
                                y = Frequency)) +
  geom_bar(stat = "identity", width = 0.75,
           fill ="#f18800",color="#f18800") +
  coord_flip() +  # 横向条形图
  theme_minimal() +
  labs(
    #title = "Superpathway Number",
    x = "Superpathway",
    y = "Frequency"
  ) +
  theme(
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    panel.grid = element_blank(), # 去除背景格线
    text = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )
dev.off()

# 加载必要的包
library(ggplot2)
library(dplyr)

# 保存为 PDF
pdf("figures/5_Subpathway_num_barplot.pdf", width = 4, height = 3.8)

# 计算 subPathway 的频数
sub_pathway_table <- as.data.frame(table(rowData$subPathway)) %>%
  arrange(desc(Freq)) %>%
  head(30)  # 仅保留前 30 项
colnames(sub_pathway_table) <- c("SubPathway", "Frequency")

# 绘制条形图
ggplot(sub_pathway_table, aes(x = reorder(SubPathway, -Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "#F1CC92", width = 0.8) + # 设置柱状条宽度
  coord_flip() +  # 横向条形图
  theme_minimal() +
  labs(
    title = "Top 30 Subpathway",
    x = "SubPathway",
    y = "Frequency"
  ) +
  theme(
    panel.background = element_blank(), # 去掉背景
    strip.background = element_blank(), # 去掉标题的背景框
    panel.grid = element_blank(),       # 去除背景格线
    text = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )

# 关闭 PDF 设备
dev.off()

# 加载必要的包
library(ggplot2)
library(dplyr)

# 创建输出目录
dir.create("figures/SuperPathway_SubPathway_Barplots", showWarnings = FALSE)

# 获取所有 unique 的 SuperPathway
super_pathways <- unique(rowData$superPathway)

# 遍历每个 SuperPathway
for (super_pathway in super_pathways) {
  # 筛选对应 SuperPathway 下的 SubPathway 数据
  sub_pathway_data <- rowData %>%
    filter(superPathway == super_pathway)
  
  # 计算 SubPathway 的频数
  sub_pathway_table <- as.data.frame(table(sub_pathway_data$subPathway)) %>%
    arrange(desc(Freq))%>%
     head(20)  # 仅保留前 30 项
  colnames(sub_pathway_table) <- c("SubPathway", "Frequency")
  
  # 如果没有数据，跳过此 SuperPathway
  if (nrow(sub_pathway_table) == 0) {
    next
  }
  
  # 保存为 PDF
  pdf(paste0("figures/SuperPathway_SubPathway_Barplots/", super_pathway, "_Subpathway_Barplot.pdf"), 
      width = 4, height = 0.65+nrow(sub_pathway_table)*0.105)
  
  # 绘制条形图
  print(ggplot(sub_pathway_table, aes(x = reorder(SubPathway, -Frequency), y = Frequency)) +
    geom_bar(stat = "identity", fill = "#F1CC92", width = 0.8) + # 设置柱状条宽度
    coord_flip() +  # 横向条形图
    theme_minimal() +
    labs(
      title = paste0("SubPathways in ", super_pathway),
      x = "SubPathway",
      y = "Frequency"
    ) +
    theme(
      panel.background = element_blank(), # 去掉背景
      strip.background = element_blank(), # 去掉标题的背景框
      panel.grid = element_blank(),       # 去除背景格线
      text = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    )
   )
  # 关闭 PDF 设备
  dev.off()
}

# 加载必要的包
library(ggplot2)
library(dplyr)

# 创建输出目录
dir.create("figures/SuperPathway_All_SubPathway_Barplots", showWarnings = FALSE)

# 获取所有 unique 的 SuperPathway
super_pathways <- unique(rowData$superPathway)

# 遍历每个 SuperPathway
for (super_pathway in super_pathways) {
  # 筛选对应 SuperPathway 下的 SubPathway 数据
  sub_pathway_data <- rowData %>%
    filter(superPathway == super_pathway)
  
  # 计算 SubPathway 的频数
  sub_pathway_table <- as.data.frame(table(sub_pathway_data$subPathway)) %>%
    arrange(desc(Freq))#%>%
     #head(20)  # 仅保留前 30 项
  colnames(sub_pathway_table) <- c("SubPathway", "Frequency")
  
  # 如果没有数据，跳过此 SuperPathway
  if (nrow(sub_pathway_table) == 0) {
    next
  }
  
  # 保存为 PDF
  pdf(paste0("figures/SuperPathway_All_SubPathway_Barplots/", super_pathway, "_Subpathway_Barplot.pdf"), 
      width = 4, height = 0.65+nrow(sub_pathway_table)*0.105)
  
  # 绘制条形图
  print(ggplot(sub_pathway_table, aes(x = reorder(SubPathway, -Frequency), y = Frequency)) +
    geom_bar(stat = "identity", fill = "#F1CC92", width = 0.8) + # 设置柱状条宽度
    coord_flip() +  # 横向条形图
    theme_minimal() +
    labs(
      title = paste0("SubPathways in ", super_pathway),
      x = "SubPathway",
      y = "Frequency"
    ) +
    theme(
      panel.background = element_blank(), # 去掉背景
      strip.background = element_blank(), # 去掉标题的背景框
      panel.grid = element_blank(),       # 去除背景格线
      text = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    )
   )
  # 关闭 PDF 设备
  dev.off()
}


