library(tidyverse)
library(pheatmap)
library(Hmisc)
library(corrplot)
library(dplyr)
library(tidyverse)
library(ggpubr)
# 加载必要的包
library(tidyr)
library(dplyr)

ab_titer <- read_csv("./vaccine_data.csv")

getwd()
results_path <- getwd()
nrow(ab_titer)
head(ab_titer)

unique(ab_titer$Age)
summary(ab_titer$Age)

ab_titer$Age_group <- "Not detected"
ab_titer$Age_group[ab_titer$Age>  20 & ab_titer$Age <= 29] <- "20s"#"Young Adults"
ab_titer$Age_group[ab_titer$Age>=  30 & ab_titer$Age <= 39] <- "30s"#"Early 30s"
ab_titer$Age_group[ab_titer$Age>=  40 & ab_titer$Age <= 49] <- "40s"#"Mid 30s to Early 40s"
ab_titer$Age_group[ab_titer$Age>=  50 & ab_titer$Age <= 59] <- "50s"#"Mid 40s to Early 50s"
ab_titer$Age_group[ab_titer$Age>=  60 & ab_titer$Age <= 69] <- "60s"#"Late 50s to Early 60s"
table(ab_titer$Age_group)
ab_titer$Age_group <- factor(ab_titer$Age_group,
                             level=c('30s','50s','40s','20s','60s'))
unique(ab_titer$Age_group)

ab_titer$Age_group <- "Not detected"
ab_titer$Age_group[ab_titer$Age>  20 & ab_titer$Age <= 29] <- "YA"#"Young Adults"
ab_titer$Age_group[ab_titer$Age>=  30 & ab_titer$Age <= 34] <- "EA"#"Early 30s"
ab_titer$Age_group[ab_titer$Age>=  35 & ab_titer$Age <= 44] <- "MAE"#"Mid 30s to Early 40s"
ab_titer$Age_group[ab_titer$Age>=  45 & ab_titer$Age <= 54] <- "MEF"#"Mid 40s to Early 50s"
ab_titer$Age_group[ab_titer$Age>=  55 & ab_titer$Age <= 65] <- "LE"#"Late 50s to Early 60s"
table(ab_titer$Age_group)

ab_titer$Sex[ab_titer$Sex == "F"] <- "Female"
ab_titer$Sex[ab_titer$Sex == "M"] <- "Male"

ab_titer$Condition <- "Not detected"
ab_titer$Condition[ab_titer$Day == 0] <- "A0"
ab_titer$Condition[ab_titer$Day>=  5 & ab_titer$Day <= 9] <- "A1"
ab_titer$Condition[ab_titer$Day>=  12 & ab_titer$Day <= 16] <- "A2"
ab_titer$Condition[ab_titer$Day>=  19 & ab_titer$Day <= 23] <- "B1"
ab_titer$Condition[ab_titer$Day>=  26 & ab_titer$Day <= 30] <- "B2"

write.csv(ab_titer, file.path(results_path, "ab_titer_week.csv"))

ab_titer <- subset(ab_titer, Condition != "Not detected")

head(ab_titer)

table(ab_titer$Day)
table(ab_titer$Condition)

ab_titer_without_A0 <- ab_titer%>% filter(Day!=0)
table(ab_titer_without_A0$Day)

pdf("1_Age_cor_ab_titer_IG_SEX.pdf",18,5)
ggplot(ab_titer, aes(x =Age, y = IgG)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x =Age, y = IgM)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x =Age, y = Total_Ig)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
dev.off()

pdf("1_Age_cor_ab_titer_IG.pdf",10,6)
ggplot(ab_titer, aes(x =Age, y = IgG)) +
  geom_point(color = 'lightblue',size = 4)+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x =Age, y = IgM)) +
  geom_point(color = 'lightblue',size = 4)+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x =Age, y = Total_Ig)) +
  geom_point(color = 'lightblue',size = 4)+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  theme_bw(base_size = 25)

dev.off()

my_comparisons1 <- list(c("A0","A1"),
                        c("A1","A2"),
                        c("A2","B1"),
                        c("B1","B2")#,
                        # c("A0","A2"),
                        #  c("A0","B1"),
                        # c("A0","B2")
                       )

#尝试画相关性的那种分图（更好看？）
new_ab_titer <- ab_titer
new_ab_titer$IgWT <- NULL

# 使用 pivot_longer 将 Total_Ig, IgG, IgM 归到 IG 列
new_ab_titer <- new_ab_titer %>%
  pivot_longer(
    cols = c(Total_Ig, IgG, IgM),
    names_to = "IG",
    values_to = "Ig_value"
  )

# 查看结果
head(new_ab_titer)

new_ab_titer$IG <- factor(new_ab_titer$IG,level=c('IgG','IgM','Total_Ig'))
unique(new_ab_titer$IG)

pdf('2_ab_titer_IG_log.pdf', 20, 6)
ggplot(new_ab_titer, aes(x = Condition, y = log(Ig_value))) +
  geom_boxplot(alpha = 1, size = 1, outlier.shape = NA, color = "#585658") +
  geom_point(aes(fill = IG), shape = 21, color = "white", size = 4, stroke = 0.6,
             alpha = 1, position = position_jitter(width = 0.23)) + 
  scale_fill_manual(values = c("Total_Ig" = "#BD3B29", "IgG" = "#0472B6", "IgM" = "#53A85F")) +
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test",
                     bracket.size = 0.8, tip.length = 0.03,
                     step.increase = 0.1,
                      size = 6) + 
  theme_bw(base_size = 25) +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    plot.title = element_text(size = 30, hjust = 0.5, face = 'bold'),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),  # 修改为 cm
    panel.border = element_rect(color = "black", size = 2.5),  # 边框加粗
    strip.text = element_text(size = 20),
      legend.title = element_text(size = 20,hjust = 0.5, face = 'bold'),
    legend.text = element_text(size = 18),
    legend.position = 'right',
    legend.key.size = unit(0.5, "inch"),
    #legend.position = "none"
  ) +
    ylim(c(-5,11.8))+
  labs(y = "log(Ig titer)", x = 'Condition') +
  facet_wrap(~IG, ncol = 3)
dev.off()

# 创建一个绘图函数，传入不同的抗体类型
plot_antibody_titer <- function(data, antibody, y_label,color1, color2) {
  ggplot(data, aes(x = Condition, y = log(get(antibody)))) +
    geom_boxplot(alpha = 0.95, size = 0.5, outlier.shape = NA, color = color1) +
    geom_point(shape = 21, color = "white", fill = color2, size = 1.8, stroke = 0.2,
               alpha = 1, position = position_jitter(width = 0.21)) +  
    stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test",
                       #symnum.args=list(cutpoints = c(0, 0.0001, 0.001,0.01, 0.05, Inf), symbols = c("****", "***", "**", "*",  " ")),
                       bracket.size = 0.3,tip.length = 0.02,step.increase =0.03,size = 2.2) + 
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 10,color="black"),
      axis.text.y = element_text(size = 10,color="black"),
      plot.title = element_text(size = 12, hjust = 0.5,face='bold'),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
      panel.border = element_rect(color = "black", size = 1),  # 边框加粗
      strip.text = element_text(size = 10),
      legend.position = "none"
    )+ 
    labs(title=antibody,y = y_label, x = 'Condition')
}

# 保存PDF并调用绘图函数
pdf(file.path(results_path, "2_ab_titer_IG_all_Log_without_SEX.pdf"), 3, 3)

plot_antibody_titer(ab_titer, "IgG", "lg(IgG titer)","#585658","#0472B6")#"#9FA3A8",
plot_antibody_titer(ab_titer, "IgM", "lg(IgM titer)","#585658","#53A85F")
plot_antibody_titer(ab_titer, "IgWT", "lg(IgWT titer)","#585658","#BD3B29")
plot_antibody_titer(ab_titer, "Total_Ig", "lg(Total_Ig titer)","#585658","#BD3B29")

dev.off()

# 创建一个绘图函数，传入不同的抗体类型
plot_antibody_titer <- function(data, antibody, y_label,color1, color2) {
  ggplot(data, aes(x = Condition, y = log(get(antibody)))) +
    geom_boxplot(alpha = 0.9, size = 0.5, outlier.shape = NA, color = color1) +
    geom_point(shape = 21, color = "white", fill = color2, size = 1.5, stroke = 0.1,
               alpha = 1, position = position_jitter(width = 0.21)) +  
    stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 2.5) + 
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
      legend.position = "none"
    ) +
    labs(title=antibody,y = y_label, x = 'Condition')
}

# 保存PDF并调用绘图函数
pdf(file.path(results_path, "2_ab_titer_IG_all_Log_without_SEX_2.pdf"), 3, 3)

plot_antibody_titer(ab_titer, "IgG", "lg(IgG titer)","#585658","#0472B6")
plot_antibody_titer(ab_titer, "IgM", "lg(IgM titer)","#585658","#53A85F")
plot_antibody_titer(ab_titer, "IgWT", "lg(IgWT titer)","#585658","#BD3B29")
plot_antibody_titer(ab_titer, "Total_Ig", "lg(Total_Ig titer)","#585658","#BD3B29")

dev.off()

pdf(file.path(results_path, "2_ab_titer_IG_all_Log.pdf"), 8,5)

ggplot(ab_titer, aes(x = Condition, y = log(IgG)))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(1.5,20))+
  ylab("log_IgG titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = log(IgM)))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("log_IgM titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = log(IgWT)))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("log_IgWT titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = log(Total_Ig)))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,100))+
  ylab("log_Total_Ig")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

pdf(file.path(results_path, "2_ab_titer_IG.pdf"), 8,5)

ggplot(ab_titer, aes(x = Condition, y = IgG))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  ylim(c(1.5,20))+
  ylab("IgG titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = IgM))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  ylim(c(0,2))+
  ylab("IgM titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = IgWT))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  ylim(c(0,2))+
  ylab("IgWT titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = Total_Ig))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  ylim(c(0,100))+
  ylab("Total_Ig")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

pdf(file.path(results_path, "2_ab_titer_IG_all.pdf"), 8,5)

ggplot(ab_titer, aes(x = Condition, y = IgG))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(1.5,20))+
  ylab("IgG titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = IgM))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("IgM titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = IgWT))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,2))+
  ylab("IgWT titer")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(ab_titer, aes(x = Condition, y = Total_Ig))+
  geom_boxplot(alpha = 0.4,  outlier.shape = NA)+
  geom_point(aes(fill=Sex), shape = 21, color = "white", size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 3)+
  theme_bw(base_size = 20)+
  #ylim(c(0,100))+
  ylab("Total_Ig")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

pdf(file.path(results_path, "3_ab_titer_IgM_corspearman.pdf"), 10,5)
ggplot(ab_titer, aes(x = IgG, y = Total_Ig)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  ylim(c(0,1))+
  xlim(c(0,20))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x = IgM, y = Total_Ig)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  xlim(c(0,20))+
  ylim(c(0,100))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x = IgWT, y = Total_Ig)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  xlim(c(0,20))+
  ylim(c(0,1))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x = IgM, y = IgG)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  xlim(c(0,20))+
  ylim(c(0,1))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  theme_bw(base_size = 25)
dev.off()

pdf(file.path(results_path, "3_ab_titer_IgM_corspearman_sex.pdf"), 18,5)
ggplot(ab_titer, aes(x = IgG, y = IgM)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  ylim(c(0,1))+
  xlim(c(0,20))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x = IgG, y = Total_Ig)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  xlim(c(0,20))+
  ylim(c(0,100))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x = IgG, y = IgWT)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  xlim(c(0,20))+
  ylim(c(0,1))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x = IgM, y = Total_Ig)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  xlim(c(0,20))+
  ylim(c(0,1))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
ggplot(ab_titer, aes(x = IgM, y = IgWT)) +
  geom_point(aes(color = Sex, shape = Sex),size = 4)+
  scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
  geom_smooth(method = "lm", linetype = "dashed", size = 1, color = "grey30")+
  xlim(c(0,20))+
  ylim(c(0,1))+
  stat_cor(size = 7, method = "spearman", cor.coef.name = "rho")+
  facet_wrap(~Sex, ncol = 2)+ 
  theme_bw(base_size = 25)
dev.off()


