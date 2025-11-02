#-----------------------------------------Spacer_population_PCoA-------------------------------------------------------------------

#load vegan for jaccard distance, ape for pcoa, ggplot2 for visualization
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(labdsv)
library(ggpubr)
library(ggExtra)
library(readxl)
library(viridis)

data = read_excel("250722_SIGisolate_spacerPOP_mtx.xlsx")
data = filter(data, Isolate_Cas_kind!="None")
k <- rowSums(data[,2:606]) > 0
data <- data[k,]

df_spacerpop = as.data.frame(data[, c(2:606)])
row.names(df_spacerpop)=data$Isolate

metadata_df <- as.data.frame(data[, c(609)])
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
plot_PCoA <- ggplot(pcoavectors_spacerpop_corr,aes(x=V1, y=V2,color=Type))+
  geom_point(size = 2)+
  scale_color_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  stat_ellipse(aes(color=Type), type = "t")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position="top") +#c(0.85,1.0)
  xlab("PCoA1 (13.20%)")+ylab("PCoA2 (8.26%)") + #jaccard T
  coord_fixed((ratio=1))
plot_PCoA
ggsave("250724_spacerpop_PCoA_JAC_SIGtype.pdf", plot_PCoA , width=5, height=5)


#Adonis test for significance (permANOVA)

df_spacerpop_meta = merge(df_spacerpop, metadata_df, by = "row.names")
adonis2(df_spacerpop ~ df_spacerpop_meta$Type,
        permutations = 999, method = 'jaccard', binary=T)

