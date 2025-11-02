## Compare the Matched spacer pops against phages in each metadata groups

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
set.seed(125)

data_all = read_excel("250722_phage_spacer_hit.xlsx", sheet="Phage_hits")
data_spacer_count = as.data.frame(table(data_all$Phage))
colnames(data_spacer_count) = c("Phage", "Spacer Pops")
phage_info = read_excel("Data_S2_v1.xlsx")
data = merge(data_spacer_count, phage_info, all = TRUE,  by = "Phage")
data$`Spacer Pops`[is.na(data$`Spacer Pops`)] <- 0

#--------------------------------------------------------------------------------------------------------------------------------------------
plot_Type <- ggplot(data, aes(x=`Clinical status`, y=`Spacer Pops`, fill=`Clinical status`)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environmental"="#B3B6B7")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  #theme_ipsum() +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Matched spacer population count") +
  xlab("Clinical status") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$`Clinical status`)), 2)), seq(nrow(t(combn(levels(as.factor(data$`Clinical status`)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=4, step_increase = .1)
# ylim(0,12.5)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_Type
ggsave("250714_Box_matched_spacerPOPCount_phage_Type.pdf", plot_Type, width=5, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------
data1 = filter(data, Group!="Household")
plot_Group <- ggplot(data1, aes(x=Group, y=`Spacer Pops`, fill=Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("AI"="#F1948A", "HI"="#CB4335","Household"="#48C9B0")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  #theme_ipsum() +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Matched spacer population count") +
  xlab("Group") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data1$Group)), 2)), seq(nrow(t(combn(levels(as.factor(data1$Group)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=4, step_increase = .1)
# ylim(0,12.5)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_Group
ggsave("250813_Box_matched_spacerPOPCount_phage_Group.pdf", plot_Group, width=2.8, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------
data = filter(data, Host == "Human")
data = filter(data, `Body site` == "nares"|`Body site` == "respiratory"|`Body site` == "tissue"|`Body site` == "wound")
plot_Bodysite <- ggplot(data, aes(x=`Body site`, y=`Spacer Pops`, fill=`Body site`)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=1, option = "E") +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  #theme_ipsum() +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Matched spacer population count") +
  xlab("Bodysite") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$`Body site`)), 2)), seq(nrow(t(combn(levels(as.factor(data$`Body site`)), 2))))), 
              map_signif_level = F, textsize=4, step_increase = .1)
# ylim(0,12.5)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_Bodysite
ggsave("250813_Box_matched_spacerPOPCount_phage_Bodysite.pdf", plot_Bodysite, width=5, height=5)

