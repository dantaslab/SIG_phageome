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

data = read_excel("230831_phage_info.xlsx")
data = filter(data, Group!="Household")
df_phagepop = as.data.frame.matrix(table(data$Isolate, data$Population))

metadata_df <- data[, c(2,4,5,6,7,8)] %>% group_by(Isolate,Type,Source,Cohort,Group,Bodysite) %>%
  summarise(.groups = 'drop') %>%
  as.data.frame()
row.names(metadata_df) = metadata_df$Isolate
head(metadata_df)

df_phagepop_meta = merge(metadata_df, df_phagepop, by = "row.names")
row.names(df_phagepop_meta) = df_phagepop_meta$Row.names
df_phagepop = select(df_phagepop_meta, -c(Row.names,Isolate,Type,Source,Cohort,Group,Bodysite))
metadata_df = select(df_phagepop_meta, c(Type,Source,Cohort,Group,Bodysite))

#Calculate jaccard distance. Set to matrix
k <- rowSums(df_phagepop) > 0
df_phagepop <- df_phagepop[k,]

df_phagepop[df_phagepop!= 0] <- 1

JAC_phagepop<-vegdist(df_phagepop, method='jaccard',binary=TRUE)
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
# write.csv(pcoavectors_phagepop_corr,file = "pcoavectors_phagepop_corr.csv")
#plot first two principal coordinate axes
plot_PCoA <- ggplot(pcoavectors_phagepop_corr,aes(x=V1, y=V2,color=as.factor(Group)))+
  geom_point(size = 2)+
  scale_color_manual(values = c("AI"="#F1948A", "HI"="#CB4335","Household"="#48C9B0")) +
  stat_ellipse(aes(color=as.factor(Group)), type = "t")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position=c(0.85,0.85))+
  # xlab("PCoA1 (11.40%)")+ylab("PCoA2 (8.73%)") + # bray
  xlab("PCoA1 (10.16%)")+ylab("PCoA2 (6.98%)") + # jaccard
  coord_fixed((ratio=1))
plot_PCoA
ggsave("250804_PhagePop_PCoA_JAC_group.pdf", plot_PCoA , width=5, height=5)
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

df_phagepop_meta = merge(df_phagepop,metadata_df, by = "row.names")
adonis2(df_phagepop ~ df_phagepop_meta$Group,
        permutations = 999, method = "jaccard")

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df_phagepop ~ df_phagepop_meta$Group, permutations = 999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# df_phagepop_meta$Group   2     2.16 0.01036 2.4022  0.001 ***
#   Residual               459   206.40 0.98964                  
# Total                  461   208.56 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df_phagepop ~ df_phagepop_meta$Group, permutations = 999, method = "jaccard")
# Df SumOfSqs      R2      F Pr(>F)    
# df_phagepop_meta$Group   2    1.982 0.00924 2.1414  0.001 ***
#   Residual               459  212.391 0.99076                  
# Total                  461  214.373 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
