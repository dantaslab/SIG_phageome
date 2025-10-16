library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
set.seed(125)

data = read_excel("Data_S1_v1.xlsx")
#--------------------------------------------------------------------------------------------------------------------------------------------
# stat_pvalue <- data %>% 
#   rstatix::wilcox_test(`vOTU richness` ~ `Clinical status`) %>%
#   filter(p < 0.05) %>% 
#   rstatix::add_significance("p") %>% 
#   rstatix::add_y_position() %>% 
#   mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))


plot_type <- ggplot(data, aes(x=`Clinical status`, y=`vOTU richness`, fill=`Clinical status`)) +
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
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$`Clinical status`)), 2)), seq(nrow(t(combn(levels(as.factor(data$`Clinical status`)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
  # ylim(0,12.5)
  #stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_type 
ggsave("250804_Box_vOTU_richness_type.pdf", plot_type, width=5, height=5)
