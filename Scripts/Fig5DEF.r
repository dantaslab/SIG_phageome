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
library(ggpubr) 

data = read_excel("Data_S6_v1.xlsx")
data$`AbiJ allele`[data$`AbiJ allele`=="NA"] = "0NA"

##------------------------------------------------------------------------------------------------------------
## Fig. 5D
plot_AbiJ_belongs3 <- ggbarstats(
  data = data,
  x = `AbiJ allele`,
  y = Belongs,
  legend.title ="Belongs",
  xlab = "Belongs",
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
  scale_fill_manual(values = c("#3d5941","#42B7B9","#ca562c","#edbb8a","#cccccc"))
  # scale_fill_manual(values = c("#FFBB5C","#8ECDDD"))
plot_AbiJ_belongs3
ggsave("231013_AbiJ_Allele_belongs_v3.pdf", plot_AbiJ_belongs3, width=5, height=5)

##------------------------------------------------------------------------------------------------------------
## Fig. 5E
Isolate_info = read_excel("Data_S1_v1.xlsx")
data = merge(data, Isolate_info, by = "Isolate")
data$`Prophage count`[data$Belongs=="Phage"] = data$`Prophage count`[data$Belongs=="Phage"] -1
plot_AbiJ_allele <- ggplot(data, aes(x=`AbiJ allele`, y=`Prophage count`, fill=`AbiJ allele`)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#cccccc","#edbb8a","#ca562c","#42B7B9","#3d5941")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Count of phages without AbiJ") +
  xlab("AbiJ allele") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$`AbiJ allele`)), 2)), seq(nrow(t(combn(levels(as.factor(data$`AbiJ allele`)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=4,step_increase = .1)
# facet_wrap(~Isolate_source)

plot_AbiJ_allele 
ggsave("231016_Box_phageCount_AbiJ_allele_v3.pdf", plot_AbiJ_allele, width=6, height=5)

##------------------------------------------------------------------------------------------------------------
## Fig. 5F
plot_AbiJ_type <- ggbarstats(
  data = data,
  x = `AbiJ allele`,
  y = `Clinical status`,
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


data_chisqt = table(data$`Clinical status`,data$`AbiJ allele`)
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


