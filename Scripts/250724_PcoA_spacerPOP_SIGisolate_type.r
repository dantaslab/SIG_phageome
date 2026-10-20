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

setwd("~/Google Drive/WashU/Phageome/SIG/Manuscript/Figures/Figure5/Fig5C")

data = read_excel("250722_SIGisolate_spacerPOP_mtx_plus_MSLT.xlsx")
data = filter(data, Isolate_Cas_kind!="None")
k <- rowSums(data[,2:606]) > 0
data <- data[k,]

df_spacerpop = as.data.frame(data[, c(2:606)])
row.names(df_spacerpop)=data$Isolate

metadata_df <- as.data.frame(data[, c(609, 612)])
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

# plot_PCoA2 <- ggMarginal(plot_PCoA, groupColour = TRUE, groupFill = TRUE, type = "boxplot")
# ggsave("MN_PhagePop_PCoA_JAC_Cohort_V2.pdf", plot_PCoA2 , width=7, height=5)

# plot_axis <- ggplot(pcoavectors_phagepop_corr, aes(x=Type, y=V1, fill=Type)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
#   geom_jitter(color="gray", size=2, alpha=1) +
#   #theme_ipsum() +
#   theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12),legend.position="none") +
#   ylab("V1 Value") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   geom_signif(comparisons = list(c("Colonizing", "Diagnostic"), c("Environment", "Colonizing"), c("Environment", "Diagnostic")),
#               map_signif_level = F, textsize=3, step_increase = .06)
# ggsave("MN_PhagePop_box_JAC_Cohort_V1.pdf", plot_axis , width=6, height=7)
# 
# plot_axis <- ggplot(pcoavectors_phagepop_corr, aes(x=Type, y=V2, fill=Type)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
#   geom_jitter(color="gray", size=2, alpha=1) +
#   #theme_ipsum() +
#   theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12),legend.position="none") +
#   ylab("V2 Value") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   geom_signif(comparisons = list(c("Colonizing", "Diagnostic"), c("Environment", "Colonizing"), c("Environment", "Diagnostic")),
#               map_signif_level = F, textsize=3, step_increase = .06)
# ggsave("MN_PhagePop_box_JAC_Cohort_V2.pdf", plot_axis , width=6, height=7)


#Adonis test for significance (permANOVA)

df_spacerpop_meta = merge(df_spacerpop, metadata_df, by = "row.names")
adonis2(df_spacerpop ~ df_spacerpop_meta$Type,
        permutations = 999, method = 'jaccard', binary=T)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df_spacerpop ~ df_spacerpop_meta$Type, permutations = 999, method = "jaccard", binary = T)
# Df SumOfSqs      R2    F Pr(>F)    
# df_spacerpop_meta$Type   2    2.009 0.01893 2.19  0.001 ***
#   Residual               227  104.142 0.98107                
# Total                  229  106.152 1.00000                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(df_spacerpop ~ df_spacerpop_meta$Type + df_spacerpop_meta$MLST,
        by = "margin", permutations = 999, method = 'jaccard', binary=T)

# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df_spacerpop ~ df_spacerpop_meta$Type + df_spacerpop_meta$MLST, permutations = 999, method = "jaccard", by = "margin", binary = T)
# Df SumOfSqs      R2        F Pr(>F)    
# df_spacerpop_meta$Type   2    0.013 0.00012   1.6103  0.008 ** 
# df_spacerpop_meta$MLST 141  103.806 0.97790 188.1688  0.001 ***
# Residual                86    0.336 0.00317                    
# Total                  229  106.152 1.00000                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


