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

data = read_excel("Data_S5_v2.xlsx", sheet = "Spacer_prophage_match")
data = filter(data, Isolates>=5)

data$Diagnostic[data$`Diagnostic_isolate%` > 75] = "Diagnostic"
data$Diagnostic[data$`Diagnostic_isolate%` <= 75] = "Non-diagnostic"

##-------------------------------------------------------------------------------------------------------------------------
## Fig. S5G
# basic histogram
p_hist <- ggplot(data, aes(x=`Diagnostic_isolate%`, fill = Diagnostic)) + 
  geom_histogram(bins=21,breaks=c(seq(0,100,5)),color="black", alpha=0.8) +
  # scale_fill_brewer(palette="RdYlBu") +
  scale_fill_manual(values = c("Diagnostic"="#EF8A62","Non-diagnostic"="#B3B6B7")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position=c(0.8,0.8)) +
  # ylab("Count of phage population") +
  # xlab("% of spacer pop from diagnostic isolates") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
p_hist
ggsave("250722_Histogram_CasSpacer_diagnostic_25_v1.pdf", p_hist, width=4.2, height=5)

##-------------------------------------------------------------------------------------------------------------------------
## Fig. 4D
data$Cas_kind_50[data$Cas_kind_50=="CasII"] = "1CasII"
data$Cas_kind_50[data$Cas_kind_50=="Both"] = "2Both"
data$Cas_kind_50[data$Cas_kind_50=="CasIII"] = "3CasIII"

plot_spacer_diagnostic <- ggbarstats(
  data = data,
  x = Cas_kind_50,
  y = Diagnostic,
  legend.title ="Cas type",
  xlab = "Spacer population source",
  ylab = "Percentage",
  package = "ggsci",
  palette = "category10_d3") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="right") +
  scale_fill_manual(values = c("1CasII"="#618264", "2Both"="#FFBB5C","3CasIII"="#9A3B3B"))
plot_spacer_diagnostic
ggsave("250807_Spacer_diagnostic_25_v1.pdf", plot_spacer_diagnostic, width=4, height=5)




