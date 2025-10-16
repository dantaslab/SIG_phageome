library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
set.seed(125)
setwd("~/Google Drive/WashU/Phageome/SIG/Manuscript/Figures/Figure1/Fig1D-F")
data = read_excel("230831_isolate_info.xlsx")
#--------------------------------------------------------------------------------------------------------------------------------------------
# stat_pvalue <- data %>% 
#   rstatix::wilcox_test(Pop_Count ~ Type) %>%
#   filter(p < 0.05) %>% 
#   rstatix::add_significance("p") %>% 
#   rstatix::add_y_position() %>% 
#   mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))


plot_type <- ggplot(data, aes(x=Type, y=Pop_Count, fill=Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("vOTU richness") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Type)), 2)), seq(nrow(t(combn(levels(as.factor(data$Type)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
  # ylim(0,12.5)
  #stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_type 
ggsave("250804_Box_vOTU_richness_type.pdf", plot_type, width=5, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------
data_group = filter(data, Group!="Houehold")
plot_group <- ggplot(data_group, aes(x=Group, y=Pop_Count, fill=Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("AI"="#F1948A", "HI"="#CB4335","Household"="#48C9B0")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("vOTU richness") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data_group$Group)), 2)), seq(nrow(t(combn(levels(as.factor(data_group$Group)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
# ylim(0,12.5)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_group 
ggsave("250804_Box_vOTU_richness_group.pdf", plot_group, width=5, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------
data = filter(data, Source == "Human")
data = filter(data, Bodysite == "nares"|Bodysite == "respiratory"|Bodysite == "tissue"|Bodysite == "wound")

plot_bodysite <- ggplot(data, aes(x=Bodysite, y=Pop_Count, fill=Bodysite)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=1, option = "E") +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("vOTU richness") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Bodysite)), 2)), seq(nrow(t(combn(levels(as.factor(data$Bodysite)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1) +
  scale_y_continuous(breaks = seq(0, 11, by = 2))
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_bodysite
ggsave("250804_Box_vOTU_richness_bodysite.pdf", plot_bodysite, width=5, height=5)
