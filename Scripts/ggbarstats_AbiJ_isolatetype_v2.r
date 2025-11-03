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

data = read_excel("231012_phagebacteria_defense.xlsx", sheet = "AbiJ")
data$Allele[data$Allele=="NA"] = "0NA"

plot_AbiJ_type <- ggbarstats(
  data = data,
  x = Allele,
  y = Isolate_type,
  legend.title ="Belongs",
  xlab = "AbiJ allele",
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
  # scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7"))
  scale_fill_manual(values = c("#3d5941","#42B7B9","#ca562c","#edbb8a","#cccccc"))
plot_AbiJ_type
ggsave("240410_AbiJ_Allele_type_v2.pdf", plot_AbiJ_type, width=5, height=5)


####################################################################

data_chisqt = table(data$Isolate_type,data$Allele)
data_chisqt

pvalue<-matrix(0,3,3)
for(row1 in 1:nrow(data_chisqt)){
  for(row2 in 1:nrow(data_chisqt)){
    print(row1)
    print(row2)
    pvalue[row1,row2] = chisq.test(rbind(data_chisqt[row1,],data_chisqt[row2,]),correct = TRUE,sim=T)$p.value
  }
}
rownames(pvalue) = rownames(data_chisqt)
colnames(pvalue) = rownames(data_chisqt)


library("gplots")
library(RColorBrewer)
Colors=brewer.pal(5,"Spectral")
breaks <- c(0.01, 0.05, 0.1,0.2, 0.5, 1)
heatmap_ast_lineage <- heatmap.2(x=pvalue,
                                 cellnote = round(pvalue,4),
                                 notecol="black",
                                 density.info="none",
                                 Colv=FALSE,
                                 Rowv=FALSE,
                                 scale="none",
                                 breaks = breaks,
                                 col=Colors, #"bluered",
                                 trace="none",
                                 main="AbiJ-Type",
                                 ylab="Type",
                                 xlab="Type")