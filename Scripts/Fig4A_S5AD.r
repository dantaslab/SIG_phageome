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
library(reshape2)

data = read_excel("Data_S1_v1.xlsx")
data$Cas_Kind[data$Cas_Kind=="None"] = "aNone"

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_type <- 
  ggbarstats(
  data = data,
  x = Cas_Kind,
  y = "Clinical status",
  legend.title ="Cas type",
  xlab = "Clinical status of Sp isolates",
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
  scale_fill_manual(values = c("CasII"="#618264", "Both"="#FFBB5C","CasIII"="#9A3B3B", "aNone"="#c1d6e9"))
  # facet_wrap(~variable, scales = "free")
plot_type
ggsave("250807_ggbar_Cas_kind_Type.pdf", plot_type, width=5, height=5)


# Type
data_chisqt = table(data$Type, data$Cas_Kind)
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
                                 main="Cas kind",
                                 ylab="",
                                 xlab="")

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------
data_group = filter(data, Group!="Houehold")
plot_group <- 
  ggbarstats(
    data = data_group,
    x = Cas_Kind,
    y = Group,
    legend.title ="Cas type",
    xlab = "Group",
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
  scale_fill_manual(values = c("CasII"="#618264", "Both"="#FFBB5C","CasIII"="#9A3B3B", "aNone"="#c1d6e9"))
# facet_wrap(~variable, scales = "free")
plot_group
ggsave("250809_ggbar_Cas_kind_Group.pdf", plot_group, width=3.6, height=5)

##-----------------------------------------------------------------------------------------------------------------------------------------------------------------
data_bodysite = filter(data, Host=="Human")
data_bodysite = filter(data_bodysite, `Body site`=="nares" | `Body site`=="tissue" | `Body site`=="wound" | `Body site`=="respiratory")

plot_bodysite <- 
  ggbarstats(
    data = data_bodysite,
    x = Cas_Kind,
    y = `Body site`,
    legend.title ="Cas type",
    xlab = "Bodysite",
    ylab = "Percentage",
    package = "ggsci",
    palette = "category10_d3") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="bottom") +
  scale_fill_manual(values = c("CasII"="#618264", "Both"="#FFBB5C","CasIII"="#9A3B3B", "aNone"="#c1d6e9"))
# facet_wrap(~variable, scales = "free")
plot_bodysite
ggsave("250809_ggbar_Cas_kind_Bodysite.pdf", plot_bodysite, width=5, height=5)




