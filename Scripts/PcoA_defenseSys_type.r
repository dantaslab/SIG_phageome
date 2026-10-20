#-----------------------------------------Phage_defense_PCoA-------------------------------------------------------------------

#load vegan for jaccard distance, ape for pcoa, ggplot2 for visualization
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(labdsv)
library(ggpubr)
library(ggExtra)
library(readxl)

data = read_excel("230903_isolate_defense_plus_MLST.xlsx")
row.names(data) = data$Isolate
df_defenseSys = data[, c(3:32)]
metadata_df = data[, c(2,35,36,37,38,39,41)]

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
  scale_color_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environmental"="#B3B6B7")) +
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
# plot_PCoA2 <- ggMarginal(plot_PCoA, groupColour = TRUE, groupFill = TRUE, type = "boxplot")
# ggsave("MN_defenseSys_PCoA_JAC_Cohort_V2.pdf", plot_PCoA2 , width=7, height=5)

# plot_axis <- ggplot(pcoavectors_defenseSys_corr, aes(x=Type, y=V1, fill=Type)) +
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
# ggsave("MN_defenseSys_box_JAC_Cohort_V1.pdf", plot_axis , width=6, height=7)
# 
# plot_axis <- ggplot(pcoavectors_defenseSys_corr, aes(x=Type, y=V2, fill=Type)) +
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
# ggsave("MN_defenseSys_box_JAC_Cohort_V2.pdf", plot_axis , width=6, height=7)


#Adonis test for significance (permANOVA)

df_defenseSys_meta = merge(df_defenseSys,metadata_df, by = "row.names")
adonis2(df_defenseSys ~ df_defenseSys_meta$Type,
        permutations = 999, method = "jaccard")

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df_defenseSys ~ df_defenseSys_meta$Type, permutations = 999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)
# df_defenseSys_meta$Type   2    0.225 0.00368 0.9049  0.524
# Residual                490   61.029 0.99632              
# Total                   492   61.254 1.00000   

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df_defenseSys ~ df_defenseSys_meta$Type, permutations = 999, method = "jaccard")
# Df SumOfSqs      R2      F Pr(>F)
# df_defenseSys_meta$Type   2    0.431 0.00464 1.1418  0.305
# Residual                490   92.375 0.99536              
# Total                   492   92.806 1.00000   

adonis2(df_defenseSys ~ df_defenseSys_meta$Type + df_defenseSys_meta$MLST,
        by = "margin", permutations = 999, method = "jaccard")

# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df_defenseSys ~ df_defenseSys_meta$Type + df_defenseSys_meta$MLST, permutations = 999, method = "jaccard", by = "margin")
# Df SumOfSqs      R2      F Pr(>F)    
# df_defenseSys_meta$Type   2    0.131 0.00141 0.5209  0.924    
# df_defenseSys_meta$MLST 307   69.382 0.74760 1.7987  0.001 ***
#   Residual                183   22.993 0.24776                  
# Total                   492   92.806 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

