library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
# library(beeswarm)
# library(ggbeeswarm)
library(reshape2)
library(rstatix)
set.seed(125)

setwd("~/Google Drive/WashU/Phageome/SIG/Manuscript/Figures/Figure1/Fig1B-C")
all_data <- read.csv('250712_Jaccard_phyloseq.csv',
                      header = T,
                      row.names = 1)
all_data$Bodysite_group[all_data$Bodysite_group!="wound"&all_data$Bodysite_group!="nares"&all_data$Bodysite_group!="tissue"&all_data$Bodysite_group!="respiratory"] = "NA"


#Household
household_data = filter(all_data, Household_group!="NA")
# beeswarm(
#   household_data$Jaccard,
#   pwcol = household_data$Household_group,
#   method = "swarm",
#   fast=TRUE)
# 
# ggplot(household_data,aes(x=1, y=Jaccard, colour=Household_group)) +
#   geom_beeswarm()

plot_household <-
ggplot(household_data, aes(x=Jaccard, color=Household_group, fill=Household_group)) +
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.3, color = "NA", bins=50)+
  geom_density(alpha=0.5)+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),legend.position=c(0.3,0.8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
plot_household
ggsave("230913_BC_phagePOP_household_violin.pdf", plot_household, width=7, height=6)

ks_test <- ks.test(household_data$Jaccard[household_data$Household_group=="Within"],
                   household_data$Jaccard[household_data$Household_group=="Between"],
                   alternative = "two.sided")
print(ks_test) # p-value < 2.2e-16
t_test <- t.test(household_data$Jaccard[household_data$Household_group=="Within"],
                 household_data$Jaccard[household_data$Household_group=="Between"],
                 paired = FALSE, alternative = "two.sided")
print(t_test) # p-value < 2.2e-16
wilcox_test <- wilcox.test(household_data$Jaccard[household_data$Household_group=="Within"],
                           household_data$Jaccard[household_data$Household_group=="Between"],
                           alternative = "two.sided") #“two.sided” (default), “greater” or “less”
print(wilcox_test) # p-value < 2.2e-16


sub_data = all_data[,c(3,18:24)]
sub_data.melt = melt(sub_data, id = "Jaccard")
sub_data.melt = filter(sub_data.melt, value!="NA")
plot_jaccard <-
  ggplot(sub_data.melt, aes(x=Jaccard, color=value, fill=value)) +
  # geom_histogram(aes(y=..density..), position="identity", alpha=0.3, color = "NA", bins=50)+
  geom_density(alpha=0.7, fill = "white")+
  # scale_fill_brewer(palette="Set2") +
  # scale_color_brewer(palette="Set2") +
  # scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  # scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),legend.position=c(0.67,0.15),legend.direction = "horizontal") +
  # guides(fill=guide_legend(ncol=2)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~variable, scales = "fixed") +
  guides(color=guide_legend(nrow=5,byrow=TRUE))
  # scale_y_log10() +
  # xlim(0.8,1)
plot_jaccard
ggsave("250712_JAC_phagePOP_density.pdf", plot_jaccard, width=8, height=7)
# npk.aov <- aov(Jaccard ~ Source_group+Type_group+Group_group+Bodysite_group+MLST_group+Household_group+Intervention_group, sub_data)
# summary(npk.aov)
# coefficients(npk.aov)

# kruskal.test(list(sub_data$Jaccard[sub_data$Source_group=="Animal"],
#                   sub_data$Jaccard[sub_data$Source_group=="Human"],
#                   sub_data$Jaccard[sub_data$Source_group=="Environment"]))

kruskal.test(Jaccard ~ Source_group, data = sub_data, subset=Source_group!="NA") # Kruskal-Wallis chi-squared = 298.95, df = 2, p-value < 2.2e-16
dunn1 = dunn_test(Jaccard ~ Source_group, data = sub_data, p.adjust.method = "BH", detailed = FALSE)
kruskal.test(Jaccard ~ Type_group, data = sub_data, subset=Type_group!="NA") # Kruskal-Wallis chi-squared = 293.31, df = 2, p-value < 2.2e-16
dunn2 = dunn_test(Jaccard ~ Type_group, data = sub_data, p.adjust.method = "BH", detailed = FALSE)
kruskal.test(Jaccard ~ Group_group, data = sub_data, subset=Group_group!="NA") # Kruskal-Wallis chi-squared = 285.86, df = 2, p-value < 2.2e-16
dunn3 = dunn_test(Jaccard ~ Group_group, data = sub_data, p.adjust.method = "BH", detailed = FALSE)
kruskal.test(Jaccard ~ Bodysite_group, data = sub_data, subset=Bodysite_group!="NA") # Kruskal-Wallis chi-squared = 21.122, df = 2, p-value = 2.591e-05
dunn4 = dunn_test(Jaccard ~ Bodysite_group, data = sub_data, p.adjust.method = "BH", detailed = FALSE)
kruskal.test(Jaccard ~ Intervention_group, data = sub_data, subset=Intervention_group!="NA") # Kruskal-Wallis chi-squared = 57.027, df = 3, p-value = 2.536e-12
dunn5 = dunn_test(Jaccard ~ Intervention_group, data = sub_data, p.adjust.method = "BH", detailed = FALSE)
kruskal.test(Jaccard ~ MLST_group, data = sub_data, subset=MLST_group!="NA") # Kruskal-Wallis chi-squared = 2217.2, df = 1, p-value < 2.2e-16
dunn6 = dunn_test(Jaccard ~ MLST_group, data = sub_data, p.adjust.method = "BH", detailed = FALSE)
kruskal.test(Jaccard ~ Household_group, data = sub_data, subset=Household_group!="NA") # Kruskal-Wallis chi-squared = 786.72, df = 1, p-value < 2.2e-16
dunn7 = dunn_test(Jaccard ~ Household_group, data = sub_data, p.adjust.method = "BH", detailed = FALSE)

dunn_result = rbind(dunn1, dunn2, dunn3, dunn4, dunn5, dunn6, dunn7)

write.csv(dunn_result, "250712_dunn_result.csv", row.names = FALSE)
