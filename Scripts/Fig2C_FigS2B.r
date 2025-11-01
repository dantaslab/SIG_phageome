library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
library(reshape2)
set.seed(125)

data_all= read_excel("Data_S4_v1.xlsx", sheet = "defense_sys_sum")
data = as.data.frame.matrix(table(data_all$isolate, data_all$type))
data$Isolate = rownames(data)

data$Gabija[data$Gabija>1] = "2>1"
data$Cas[data$Cas>1] = "2>1"

data_prophage= read_excel("Data_S1_v1.xlsx")
data = merge(data, data_prophage[,c(1,12)], by = "Isolate")

data_melt = melt(data, id=c("Isolate", "Prophage count"), variable.name = "Defense_system", value.name = "Defense_Count")

# data_melt = filter(data_melt, Defense_system!="Cas" & Defense_system!="RosmerTA")

# Each defense system vs phage count
plot_Phagecount <- ggplot(data_melt, aes(x=as.factor(Defense_Count), y=`Prophage count`, fill=as.factor(Defense_Count))) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.5) +
  # geom_violin(adjust = 2) +
  scale_fill_manual(values = c("0"="#E0F7E1", "1"="#48C2B4", "2"="#007FA8", "3"="#423460", "4"="#070707")) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 2)+
  # geom_jitter(color="gray", size=0.3, alpha=0.3, height = 0) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12),legend.position="none") +
  ylab("Phage count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons =  list(c("0","1")), # list(split(t(combn(levels(as.factor(data_melt$Defense_Count)), 2)), seq(nrow(t(combn(levels(as.factor(data_melt$Defense_Count)), 2)))))),
              map_signif_level = F, textsize=3, step_increase = .06) +
  facet_wrap(~Defense_system, scales = "free", nrow = 3) +
  ylim(0,7)
plot_Phagecount
ggsave("250813_defense_phageCount_boxplot_v1.pdf", plot_Phagecount , width=16, height=7.5)
