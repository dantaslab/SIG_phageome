library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
library(vegan)
library(labdsv)
set.seed(126)

data_all= read_excel("Data_S4_v1.xlsx", sheet = "defense_sys_sum")
data = as.data.frame.matrix(table(data_all$isolate, data_all$type))
data$Total = rowSums(data)
data$Isolate = rownames(data)

data_meta = read_excel("Data_S1_v1.xlsx")
data = merge(data, data_meta[,c(1,3)], by = "Isolate")
colnames(data)[33] = "Type"
rownames(data) = data$Isolate

##--------------------------------------------------------------------------------------------------------------------------------------------
## Fig. 2C
df_defenseSys = data[, c(2:31)]
metadata_df = data[, c(1,33)]

#Calculate jaccard distance. Set to matrix
k <- rowSums(df_defenseSys) > 0
df_defenseSys <- df_defenseSys[k,]

df_defenseSys[df_defenseSys!= 0] <- 1

JAC_defenseSys<-vegdist(df_defenseSys, method='jaccard',binary=F)
JAC_defenseSys<-as.matrix(JAC_defenseSys)

# #Make PCoA (correction is to account for negative eigenvalues)
pcoa_defenseSys_corr <- pco(JAC_defenseSys, k = 4)
pcoavectors_defenseSys_corr <- as.data.frame(pcoa_defenseSys_corr$points)

#Get % variance captured by each axis. Rel_corr_eig = Relative eigenvalues following correction method. Sums to 100
eigen_BC_defenseSys <- pcoa_defenseSys_corr$eig
eigen_BC_defenseSys[eigen_BC_defenseSys < 0] <- 0
rel_eigen_defenseSys_corr<-as.data.frame(eigen_BC_defenseSys*100/sum(eigen_BC_defenseSys))
rel_eigen_defenseSys_corr[1:2,]

#Scree plot
# plot_bar <- barplot(rel_eigen_defenseSys_corr[1:10,], ylab="Variance Explained", xlab="Principal Components")
# plot_bar

pcoavectors_defenseSys_corr = merge(pcoavectors_defenseSys_corr,metadata_df, by = "row.names")

#plot first two principal coordinate axes
plot_PCoA <- ggplot(pcoavectors_defenseSys_corr,aes(x=V1, y=V2,color=as.factor(Type)))+
  geom_point(size = 2)+
  scale_color_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  stat_ellipse(aes(color=as.factor(Type)), type = "t")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position=c(0.85,0.85))+
  # xlab("PCoA1 (16.10%)")+ylab("PCoA2 (11.18%)") + # bray
  xlab("PCoA1 (15.21%)")+ylab("PCoA2 (11.86%)") + # jaccard
  coord_fixed((ratio=1))
plot_PCoA
ggsave("230903_DefenseSys_PCoA_JAC_type.pdf", plot_PCoA , width=5, height=5)


df_defenseSys_meta = merge(df_defenseSys,metadata_df, by = "row.names")
adonis2(df_defenseSys ~ df_defenseSys_meta$Type,
        permutations = 999, method = "jaccard")

##--------------------------------------------------------------------------------------------------------------------------------------------
## Fig. 2D
stat_pvalue <- data %>%
  rstatix::wilcox_test(Total ~ Type) %>%
  filter(p < 0.05) %>%
  rstatix::add_significance("p") %>%
  rstatix::add_y_position() %>%
  mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

plot_type <- ggplot(data, aes(x=Type, y=Total, fill=Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environmental"="#B3B6B7")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Phage defense system count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Type)), 2)), seq(nrow(t(combn(levels(as.factor(data$Type)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), # stat_pvalue$groups, #
              map_signif_level = F, textsize=5, step_increase = .1)
# ylim(0,12.5)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_type 
ggsave("250804_Box_sysCount_type_v2.pdf", plot_type, width=4.2, height=5)
