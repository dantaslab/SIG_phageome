library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
library(labdsv)
library(vegan)
set.seed(125)

data_all = read_excel("250722_SIGisolate_spacerPOP_mtx.xlsx")
data = filter(data_all, Isolate_Cas_kind!="None")
data$spacerPOP_count = rowSums(data[,2:606])

summary_by_group <- data %>%
  group_by(Isolate_Cas_kind) %>%
  summarize(
    Mean = mean(spacerPOP_count, na.rm = TRUE),
    Median = median(spacerPOP_count, na.rm = TRUE)
  )


##--------------------------------------------------------------------------------------------------------------------------------------------
## Fig. 3C
plot_spacer_count <- ggplot(data, aes(x=Isolate_Cas_kind, y=spacerPOP_count, fill=Isolate_Cas_kind)) +
  geom_boxplot() +
  scale_fill_manual(values = c("CasII"="#618264", "Both"="#FFBB5C","CasIII"="#9A3B3B")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 2)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Spacer population count") +
  xlab("Isolate Cas kind") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Isolate_Cas_kind)), 2)), seq(nrow(t(combn(levels(as.factor(data$Isolate_Cas_kind)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
# ylim(0,12.5)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_spacer_count
ggsave("250722_Box_spacerPOPCount_CasType.pdf", plot_spacer_count, width=4, height=5)


##--------------------------------------------------------------------------------------------------------------------------------------------
## Fig. 3D
data = filter(data_all, Isolate_Cas_kind!="None")
k <- rowSums(data[,2:606]) > 0
data <- data[k,]

df_spacerpop = as.data.frame(data[, c(2:606)])
row.names(df_spacerpop)=data$Isolate

metadata_df <- as.data.frame(data[, c(607)])
row.names(metadata_df)=data$Isolate

#Calculate jaccard distance. Set to matrix
#df_spacerpop[df_spacerpop!= 0] <- 1

JAC_spacerpop<-vegdist(df_spacerpop, method='jaccard', binary=T)
JAC_spacerpop<-as.matrix(JAC_spacerpop)

# #Make PCoA (correction is to account for negative eigenvalues)
pcoa_spacerpop_corr <- pco(JAC_spacerpop, k = 4)
pcoavectors_spacerpop_corr <- as.data.frame(pcoa_spacerpop_corr$points)

#Get % variance captured by each axis. Rel_corr_eig = Relative eigenvalues following correction method. Sums to 100
eigen_BC_spacerpop <- pcoa_spacerpop_corr$eig
eigen_BC_spacerpop[eigen_BC_spacerpop < 0] <- 0
rel_eigen_spacerpop_corr<-as.data.frame(eigen_BC_spacerpop*100/sum(eigen_BC_spacerpop))
rel_eigen_spacerpop_corr[1:2,]

#Scree plot
# plot_bar <- barplot(rel_eigen_spacerpop_corr[1:10,], ylab="Variance Explained", xlab="Principal Components")
# plot_bar

pcoavectors_spacerpop_corr = merge(pcoavectors_spacerpop_corr, metadata_df, by = "row.names")

#plot first two principal coordinate axes
plot_PCoA <- ggplot(pcoavectors_spacerpop_corr,aes(x=V1, y=V2,color=Isolate_Cas_kind))+
  geom_point(size = 2)+
  scale_color_manual(values = c("CasII"="#618264", "Both"="#FFBB5C","CasIII"="#9A3B3B")) +
  stat_ellipse(aes(color=Isolate_Cas_kind), type = "t")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position="top")+
  xlab("PCoA1 (13.20%)")+ylab("PCoA2 (8.26%)") + #jaccard T
  coord_fixed((ratio=1))
plot_PCoA
ggsave("250722_spacerpop_PCoA_JAC_Cas.pdf", plot_PCoA , width=5, height=5)

#Adonis test for significance (permANOVA)
df_spacerpop_meta = merge(df_spacerpop, metadata_df, by = "row.names")
adonis2(df_spacerpop ~ df_spacerpop_meta$Isolate_Cas_kind,
        permutations = 999, method = 'jaccard', binary=T)





