##Recursive feature elimination using caret
library(ggplot2)
library(caret)
library(vegan)
library(dendextend)
library(gplots)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library("compositions")
library(ggstatsplot)
library(readxl)

data = read_excel("Data_S4_v1.xlsx", sheet = "defense_sys_sum")
data$Source[data$`Phage?`==1] = "Phage"
data$Source[data$`Phage?`==0] = "Bacteria"
data = filter(data, type=="Abi2" | type=="AbiD" | type=="AbiJ" | type=="AbiN" | type=="Gabija" | type=="Kiwa"
              | type=="PD-Lambda-1" | type=="PD-T4-9")

##-------------------------------------------------------------------------------------------------------------------------

plot_defense <- ggbarstats(
  data = data,
  x = Source,
  y = type,
  legend.title ="Source",
  xlab = "Antiphage defense system",
  ylab = "Percentage",
  package = "ggsci",
  palette = "category10_d3") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="top") +
  scale_fill_manual(values = c("Bacteria"="#8ECDDD", "Phage"="#FFBB5C"))
plot_defense
ggsave("250807_phage_defense.pdf", plot_defense, width=6, height=5)




