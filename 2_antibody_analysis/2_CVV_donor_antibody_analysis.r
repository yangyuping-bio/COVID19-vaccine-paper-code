library(tidyverse)
library(pheatmap)
library(Hmisc)
library(corrplot)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(openxlsx)

getwd()
results_path <- getwd()

ab_titer <- read.xlsx("donor_antibody_data.xlsx")
unique(ab_titer$Donor)
nrow(ab_titer)
head(ab_titer)

ab_titer$Donor <- factor(
  ab_titer$Donor,
  levels = c("Donor1", "Donor2", "Donor3",
             "Donor4", "Donor5", "Donor6",
            "Donor7", "Donor8", "Donor9",
            "Donor10", "Donor11", "Donor12"))
colnames(ab_titer) <- c("Donor","Sex","Condition","Total_Ig_new","Total_Ig","IgG","IgM")
head(ab_titer)

#尝试画相关性的那种分图（更好看？）
new_ab_titer <- ab_titer
new_ab_titer$Total_Ig <- new_ab_titer$Total_Ig_new
new_ab_titer$Total_Ig_new <- NULL
# 使用 pivot_longer 将 Total_Ig, IgG, IgM 归到 IG 列
new_ab_titer <- new_ab_titer %>%
  pivot_longer(
    cols = c(Total_Ig, IgG, IgM),
    names_to = "IG",
    values_to = "Ig_value"
  )

new_ab_titer$IG <- factor(new_ab_titer$IG,level=c('IgG','IgM','Total_Ig'))
unique(new_ab_titer$IG)
# 查看结果
head(new_ab_titer)

my_comparisons1 <- list(c("A0","A1"),
                        c("A1","A2"),
                        c("A2","B0"),
                        c("B0","B1"),
                        c("B1","B2"),
                        c("B2","C0"),
                        c("C0","C1"),
                        c("C1","C2")
                       )

ncols <- c( "#51A3CCFF", "#f18800","#85B22CFF","#B22C2CFF",
            "#E57E7EFF", "#BFB2FFFF","#e20612", "#ffd401", 
             "#E5B17EFF", "#942d8d","#573333FF", "#2e409a")

head(ab_titer)

pdf('2_ab_titer_IG_Total_Ig_new_2.pdf', 4.5, 3)
 ggplot(ab_titer, aes(x = Condition, y = log(Total_Ig_new))) +
    geom_boxplot(alpha = 0.9, size = 0.5, outlier.shape = NA, color = "#585658") +
    geom_point(aes(fill = Donor),shape = 21, color = "white", size = 2, stroke = 0.2,
               alpha = 1, position = position_jitter(width = 0.23)) +  
   scale_fill_manual(values = ncols) +
    stat_compare_means(comparisons = my_comparisons1, 
                       step.increase =0.07,bracket.size = 0.3,tip.length = 0.018,
                       method = "wilcox.test", size = 2.5) + 
    theme_bw(base_size = 10) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      axis.text.x = element_text(size = 10,color="black"),
      axis.text.y = element_text(size = 10,color="black"),
      plot.title = element_text(size = 12, hjust = 0.5,face='bold'),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
        legend.title = element_text(size = 10, hjust = 0.5,face='bold'),
        legend.text = element_text(size = 8),
        legend.position = 'right',
        legend.key.size = unit(0.15, "inch")
    ) +
    labs(title="Total_Ig",y = "lg(Total_Ig titer)", x = 'Condition')

dev.off()

pdf('2_ab_titer_IG_Total_Ig_new.pdf', 4.5, 3)
 ggplot(ab_titer, aes(x = Condition, y = log(Total_Ig_new))) +
    geom_boxplot(alpha = 0.95, size = 0.3, outlier.shape = NA, color = "#585658") +
    geom_point(aes(fill = Donor),shape = 21, color = "white",  size = 2.3, stroke = 0.3,
               alpha = 1, position = position_jitter(width = 0.23)) +  
    scale_fill_manual(values = ncols) +
    stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test",
                       bracket.size = 0.3,tip.length = 0.018,step.increase =0.028,size = 2.2) + 
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 10,color="black"),
      axis.text.y = element_text(size = 10,color="black"),
      plot.title = element_text(size = 12, hjust = 0.5,face='bold'),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
      panel.border = element_rect(color = "black", size = 1.1),  # 边框加粗
      strip.text = element_text(size = 10),
        legend.title = element_text(size = 10, hjust = 0.5,face='bold'),
        legend.text = element_text(size = 8),
        legend.position = 'right',
        legend.key.size = unit(0.15, "inch")
    ) +
    labs(title="Total_Ig",y = "lg(Total_Ig titer)", x = 'Condition')

dev.off()

pdf('2_ab_titer_IG_log.pdf', 24, 6)
ggplot(new_ab_titer, aes(x = Condition, y = log(Ig_value))) +
  geom_boxplot(alpha = 1, size = 1, outlier.shape = NA, color = "#585658") +
  geom_point(aes(fill = Donor), shape = 21, color = "white", size = 4, stroke = 0.6,
             alpha = 1, position = position_jitter(width = 0.23)) + 
  scale_fill_manual(values = ncols) +
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test",
                     bracket.size = 0.8, tip.length = 0.03,
                     step.increase = 0.05,
                      size = 5) + 
  theme_bw(base_size = 26) +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 30, hjust = 0.5, face = 'bold'),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),  # 修改为 cm
    panel.border = element_rect(color = "black", size = 2.5),  # 边框加粗
    strip.text = element_text(size = 20)#,
    #legend.position = "none"
  ) +
    ylim(c(-5,12))+
  labs(y = "log(Ig titer)", x = 'Condition') +
  facet_wrap(~IG, ncol = 3,scales  ='free')
dev.off()

pdf(file.path(results_path, "ab_titer_IG_Total.pdf"), 8,5)
ggplot(ab_titer, aes(x = Condition, y = Total_Ig_new))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("Total_Ig_new titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(file.path(results_path, "ab_titer_IG.pdf"), 8,5)

ggplot(ab_titer, aes(x = Condition, y = IgG))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(1.5,20))+
  ylab("IgG titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = IgM))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("IgM titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = Total_Ig_new))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("Total_Ig_new titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = Total_Ig))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,100))+
  ylab("Total_Ig")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

head(ab_titer)

log10(10)

ab_titer <- ab_titer %>% 
      mutate(log_Total_Ig_new=log10(Total_Ig_new),
             log_IgG=log10(IgG),
             log_IgM=log10(IgM)
            )

head(ab_titer)

pdf(file.path(results_path, "ab_titer_IG_log10.pdf"), 8,5)

ggplot(ab_titer, aes(x = Condition, y = log_IgG))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(1.5,20))+
  ylab("log_IgG titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = log_IgM))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("log_IgM titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = log_Total_Ig_new))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  # scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  # scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("log_Total_Ig_new titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(file.path(results_path, "ab_titer_IG_Total_log10.pdf"), 8,5)
ggplot(ab_titer, aes(x = Condition, y = log_Total_Ig_new))+
  geom_boxplot(alpha = 0.4, outlier.shape = NA)+
  geom_point(aes(fill=Donor), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  geom_line(aes(group = Donor), color = "lightgrey", size = 0.3) + # 添加数据连接线
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  ylab("log_Total_Ig_new titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




