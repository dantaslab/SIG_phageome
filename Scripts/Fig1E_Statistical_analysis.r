library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
library(ggpattern) #remotes::install_github("coolbutuseless/ggpattern")
set.seed(125)

setwd("~/Google Drive/WashU/Phageome/SIG/Manuscript/Figures/Figure1/Others")
df_all = read_excel("SIGisolate_certainAMG_matrix.xlsx")
df_all$Total = rowSums(df_all[,2:25])
df_all_AMG = as.matrix(df_all[,c(2:25,38)])
row.names(df_all_AMG) = df_all$Isolate
data_all_AMG = reshape2::melt(df_all_AMG, c("Isolate", "AMG"), value.name = "Count")
head(data_all_AMG)
#data_all_AMG$AMG_per_phage = rep(df_all_AMG$AMG_per_phage,24)
data_all_AMG$Belongs = rep(c("All"),25*493)
head(data_all_AMG)

df_bac = read_excel("231006_SIGisolate_bacAMG_matrix.xlsx")
df_bac$Total = rowSums(df_bac[,2:25])
df_bac_AMG = as.matrix(df_bac[,c(2:25,38)])
row.names(df_bac_AMG) = df_bac$Isolate
data_bac_AMG = reshape2::melt(df_bac_AMG, c("Isolate", "AMG"), value.name = "Count")
head(data_bac_AMG)
#data_all_AMG$AMG_per_phage = rep(df_bac_AMG$AMG_per_phage,24)
data_bac_AMG$Belongs = rep(c("Bacteria"),25*493)
head(data_bac_AMG)

data = rbind(data_all_AMG, data_bac_AMG)

#--------------------------------------------------------------------------------------------------------------------------------------------
plot_AMGcount <- ggplot(data, aes(x=AMG, y=Count, fill=AMG, pattern=Belongs)) +
  geom_boxplot() +
  geom_boxplot_pattern(color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  scale_pattern_manual(values = c("All" = "stripe", "Bacteria" = "none")) +
  scale_fill_manual(values = c(rep("white",24))) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("AMG count") +
  scale_x_discrete(limits = c('K01047','K01495','K00390','K02788','K00248','K00558','K01951','K06167','K10026','K06920','K01737','K01915','K00657','K22477','K03800','K00382','K01220','K20895','K01069','K01776','K00459','K11212','K15634','K03476')) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_compare_means(method = "wilcox.test", label ="p.format", label.x = 1.5) #p.signif
  # ylim(0,12.5)
  #stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_AMGcount
ggsave("231009_Box_AMGcount_all_vs_bacteria.pdf", plot_AMGcount, width=18, height=3)

#--------------------------------------------------------------------------------------------------------------------------------------------

data$Type  = rep(df_all$Type,25)
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}

fun_median <- function(x){
  return(data.frame(y=median(x),label=median(x,na.rm=T)))}
# 
# data$Type[data$Type=="Colonizing"] = "Household"
# data$Type[data$Type=="Environment"] = "Household"

plot_Type <- ggplot(data, aes(x=Type, y=Count, fill=Type, pattern = Belongs)) +
  # geom_boxplot() +
  geom_violin(position = position_dodge(width=0.9))+
  geom_violin_pattern(color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c("All" = "stripe", "Bacteria" = "none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  # geom_point(color="gray", size=0.8, alpha=1, position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0.25)) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10),legend.position="right") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5, size = 4) + #label.x = 1.5,
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Type)), 2)), seq(nrow(t(combn(levels(as.factor(data$Type)), 2))))), #stat_pvalue$Types, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = T, textsize=3, step_increase = .1, vjust=0.2) +
  stat_summary(fun.y = mean, geom="point",colour="black", size=2, position=position_dodge(width=0.9)) +
  # stat_summary(fun.data = fun_mean, geom="text", position=position_dodge(width=1))
  facet_wrap(~AMG, scale="free")
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_Type 
ggsave("231009_Box_individual_AMGcount_Type.pdf", plot_Type, width=10, height=10)

#--------------------------------------------------------------------------------------------------------------------------------------------

data = filter(data, AMG == "K00382" | AMG == "K00558" | AMG == "K01915" | AMG == "K03800")
data$Type  = rep(df_all$Type,4)
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}

fun_median <- function(x){
  return(data.frame(y=median(x),label=median(x,na.rm=T)))}
# 
# data$Type[data$Type=="Colonizing"] = "Household"
# data$Type[data$Type=="Environment"] = "Household"

plot_Type <- ggplot(data, aes(x=Type, y=Count, fill=Type, pattern = Belongs)) +
  # geom_boxplot() +
  geom_violin(position = position_dodge(width=0.9))+
  geom_violin_pattern(color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c("All" = "stripe", "Bacteria" = "none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  # geom_point(color="gray", size=0.8, alpha=1, position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0.25)) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10),legend.position="none") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5, size = 4) + #label.x = 1.5,
  # geom_signif(comparisons = split(t(combn(levels(as.factor(data$Type)), 2)), seq(nrow(t(combn(levels(as.factor(data$Type)), 2))))), #stat_pvalue$Types, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
  #             map_signif_level = T, textsize=3, step_increase = .1, vjust=0.2) +
  stat_summary(fun.y = mean, geom="point",colour="black", size=2, position=position_dodge(width=0.9)) +
  # stat_summary(fun.data = fun_mean, geom="text", position=position_dodge(width=1))
  facet_wrap(~AMG, scale="free",nrow = 1)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_Type 
ggsave("231009_Box_individual_AMGcount_Type.pdf", plot_Type, width=11, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------

data = filter(data, AMG == "Total")
data$Type = c(df_all$Type, df_bac$Type)

fun_mean <- function(x){
  return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}

fun_median <- function(x){
  return(data.frame(y=median(x),label=median(x,na.rm=T)))}
# 
# data$Type[data$Type=="Colonizing"] = "Household"
# data$Type[data$Type=="Environment"] = "Household"

plot_Type2 <- ggplot(data, aes(x=Type, y=Count, fill=Type, pattern = Belongs)) +
  # geom_boxplot() +
  geom_violin(position = position_dodge(width=0.9))+
  geom_violin_pattern(color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c("All" = "stripe", "Bacteria" = "none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  # geom_point(color="gray", size=0.8, alpha=1, position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0.25)) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10),legend.position="right") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5, size = 4) + #label.x = 1.5,
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Type)), 2)), seq(nrow(t(combn(levels(as.factor(data$Type)), 2))))), #stat_pvalue$Types, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = T, textsize=3, step_increase = .1, vjust=0.2) +
  stat_summary(fun.y = mean, geom="point",colour="black", size=2, position=position_dodge(width=0.9))
  # stat_summary(fun.data = fun_mean, geom="text", position=position_dodge(width=1))

#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_Type2
ggsave("231009_Box_total_AMGcount_Type.pdf", plot_Type2, width=5, height=5)

