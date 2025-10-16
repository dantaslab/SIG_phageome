#-----------------------------------------Phage_population_PCoA-------------------------------------------------------------------

#load vegan for jaccard distance, ape for pcoa, ggplot2 for visualization
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(labdsv)
library(ggpubr)
library(ggExtra)
library(readxl)

data = read_excel("Data_S2_v1.xlsx")
df_phagepop = as.data.frame.matrix(table(data$Isolate, data$vOTU))

metadata_df <- data[, c(2,4,5,6)] %>% group_by(Isolate,`Clinical status`,Host,Cohort) %>%
  summarise(.groups = 'drop') %>%
  as.data.frame()
row.names(metadata_df) = metadata_df$Isolate
head(metadata_df)

df_phagepop_meta = merge(metadata_df, df_phagepop, by = "row.names")
row.names(df_phagepop_meta) = df_phagepop_meta$Row.names
df_phagepop = select(df_phagepop_meta, -c(Row.names,Isolate,`Clinical status`,Host,Cohort))
metadata_df = select(df_phagepop_meta, c(`Clinical status`,Host,Cohort))

#Calculate jaccard distance. Set to matrix
k <- rowSums(df_phagepop) > 0
df_phagepop <- df_phagepop[k,]

df_phagepop[df_phagepop!= 0] <- 1

JAC_phagepop<-vegdist(df_phagepop, method='jaccard',binary=F)
JAC_phagepop<-as.matrix(JAC_phagepop)

# #Make PCoA (correction is to account for negative eigenvalues)
pcoa_phagepop_corr <- pco(JAC_phagepop, k = 4)
pcoavectors_phagepop_corr <- as.data.frame(pcoa_phagepop_corr$points)

#Get % variance captured by each axis. Rel_corr_eig = Relative eigenvalues following correction method. Sums to 100
eigen_BC_phagepop <- pcoa_phagepop_corr$eig
eigen_BC_phagepop[eigen_BC_phagepop < 0] <- 0
rel_eigen_phagepop_corr<-as.data.frame(eigen_BC_phagepop*100/sum(eigen_BC_phagepop))
rel_eigen_phagepop_corr[1:2,]

#Scree plot
# plot_bar <- barplot(rel_eigen_phagepop_corr[1:10,], ylab="Variance Explained", xlab="Principal Components")
# plot_bar

pcoavectors_phagepop_corr = merge(pcoavectors_phagepop_corr,metadata_df, by = "row.names")

#plot first two principal coordinate axes
plot_PCoA <- ggplot(pcoavectors_phagepop_corr,aes(x=V1, y=V2,color=as.factor(`Clinical status`)))+
  geom_point(size = 2)+
  scale_color_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  stat_ellipse(aes(color=as.factor(`Clinical status`)), type = "t")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position=c(0.85,0.85))+
  # xlab("PCoA1 (11.40%)")+ylab("PCoA2 (8.73%)") + # bray
  xlab("PCoA1 (10.11%)")+ylab("PCoA2 (8.14%)") + # jaccard
  coord_fixed((ratio=1))
plot_PCoA
ggsave("230831_PhagePop_PCoA_JAC_type.pdf", plot_PCoA , width=5, height=5)


#Adonis test for significance (permANOVA)

df_phagepop_meta = merge(df_phagepop,metadata_df, by = "row.names")
adonis2(df_phagepop ~ df_phagepop_meta$`Clinical status`,
        permutations = 999, method = "jaccard")

