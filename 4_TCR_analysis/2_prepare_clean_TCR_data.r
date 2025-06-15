# 0.环境准备 --------------------------------------------------------------------
library("tidyverse")


getwd()

# 设置路径
results_path <- getwd()
data_save_path <- getwd()

CD8T <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD8T_TCR.rds")
CD4T <- readRDS("/share/home/qlab/projects/project_cvv/yyp_results_3/6_TCR_new/1_TCR_info/data/1_CD4T_TCR.rds")

colnames(CD8T@meta.data)[c(19:25,80,81,85:123)]
colnames(CD4T@meta.data)[c(19:25,80:120)]

CD8T.df <- as.data.frame(CD8T@meta.data[ ,colnames(CD8T@meta.data)[c(19:25,80,81,85:123)]])
write.csv(CD8T.df,  "data/0_CD8T_TCR_df.csv")
table(CD8T.df$tcr_chain)

CD4T.df <- as.data.frame(CD4T@meta.data[ ,colnames(CD4T@meta.data)[c(19:25,80:120)]])
write.csv(CD4T.df,  "data/0_CD4T_TCR_df.csv")
table(CD4T.df$tcr_chain)

head(CD4T.df$tcr_cdr3s_aa)

nrow(CD8T.df)
nrow(CD4T.df)
T_cell_df <- rbind(CD8T.df, CD4T.df)
nrow(T_cell_df)
colnames(T_cell_df)

# 提取第一条TRB链的函数
extract_first_trb <- function(cdr3_string) {
  # 使用正则表达式查找所有TRB链
  trb_matches <- unlist(regmatches(cdr3_string, gregexpr("TRB:[^;]+", cdr3_string)))
  # 如果存在TRB链，返回第一条；否则返回NA
  if (length(trb_matches) > 0) {
    return(trb_matches[1])
  } else {
    return(NA)
  }
}

# 应用函数提取所有数据中的第一条TRB链
T_cell_df$tcr_TRB_cdr3<- sapply(T_cell_df$tcr_cdr3s_aa, extract_first_trb)

# 查看结果
head(T_cell_df$tcr_TRB_cdr3)

T_cell_df$cdr3_b_aa <-sub("TRB:", "", T_cell_df$tcr_TRB_cdr3)
T_cell_df$j_b_gene <- T_cell_df$tcr_j_gene
T_cell_df$v_b_gene <- T_cell_df$tcr_v_gene
nrow(T_cell_df)
T_cell_df <- T_cell_df%>% 
           filter(!is.na(cdr3_b_aa))
nrow(T_cell_df)
head(T_cell_df$cdr3_b_aa)

colnames(T_cell_df) <-  gsub("tcr_","",colnames(T_cell_df))
colnames(T_cell_df)

write.csv(T_cell_df,  "data/1_Tcell_tcr_meta_info.csv")

#每一种序列的统计数
clonotype_num <- as.data.frame(table(T_cell_df$cdr3_b_aa))
table(clonotype_num$Freq)
T_cell_df$Clonotype_num_b <- clonotype_num$Freq[match(T_cell_df$cdr3_b_aa, clonotype_num$Var1)]

#每个 Condition 下每一种序列的统计数
clonotype_num_Condition <- T_cell_df[ ,c("cdr3_b_aa","Condition")] %>% 
                  group_by(cdr3_b_aa,Condition) %>%
                  mutate(Clonotype_num_b_Condition = n())
nrow(T_cell_df)
nrow(clonotype_num_Condition)
head(clonotype_num_Condition)

## 筛选
nrow(T_cell_df)
T_cell_df$Clontype_num_b_Condition <- clonotype_num_Condition$Clonotype_num_b_Condition
T_cell_df <-  T_cell_df%>%  filter (!(Condition =="A0" & Clontype_num_b_Condition >1))
nrow(T_cell_df)

write.csv(T_cell_df,  "data/1_Tcell_tcr_meta_filter.csv")




