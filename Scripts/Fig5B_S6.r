library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
set.seed(125)

setwd("~/Google Drive/WashU/Phageome/SIG/Manuscript/Figures/Figure3/Figure3C")
df = read.csv("230924_phagebacteria_defense.csv",
                head = TRUE,
                row.names = 1)
df = filter(df, Delphini == FALSE)

##-------------------------------------------------------------------------------------------------------------------------
## Fig. 5B
data1 = filter(df, bAbiJ>=1&pAbiJ==0)
data1$Belongs = replicate(dim(data1)[1], "Bacteria")
data2 = filter(df, pAbiJ>=1&bAbiJ==0)
data2$Belongs = replicate(dim(data2)[1], "Phage")
data2$Phage = data2$Phage-data2$pAbiJ
data3 = filter(df, bAbiJ==0&pAbiJ==0)
data3$Belongs = replicate(dim(data3)[1], "zNone")
data4 = filter(df, bAbiJ>=1&pAbiJ>=1)
data4$Belongs = replicate(dim(data4)[1], "aBoth")
data4$Phage = data4$Phage-data4$pAbiJ
data  = rbind(data1,data2,data3,data4)

#AbiJ vs phage count
plot_Phagecount <- ggplot(data, aes(x=Belongs, y=Phage, fill=Belongs)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Bacteria"="#8ECDDD", "Phage"="#FFBB5C", "zNone" = "white", "aBoth" = "#AE445A")) +
  # geom_jitter(color="gray", size=2, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),legend.position="none") +
  ylab("Count of phage without AbiJ") +
  xlab("AbiJ") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons =  split(t(combn(levels(as.factor(data$Belongs)), 2)), seq(nrow(t(combn(levels(as.factor(data$Belongs)), 2))))), #list(c("Bacteria","Phage")), 
              map_signif_level = F, textsize=4, step_increase = .06)
plot_Phagecount
ggsave("250807_AbiJ_phagecount_boxplot_v1.pdf", plot_Phagecount , width=5, height=5)

##-------------------------------------------------------------------------------------------------------------------------
## Fig. S6A
data1 = filter(df, bAbi2>=1&pAbi2==0)
data1$Belongs = replicate(dim(data1)[1], "Bacteria")
data2 = filter(df, pAbi2>=1&bAbi2==0)
data2$Belongs = replicate(dim(data2)[1], "Phage")
data2$Phage = data2$Phage-data2$pAbi2
data3 = filter(df, bAbi2==0&pAbi2==0)
data3$Belongs = replicate(dim(data3)[1], "zNone")
data4 = filter(df, bAbi2>=1&pAbi2>=1)
data4$Belongs = replicate(dim(data4)[1], "aBoth")
data4$Phage = data4$Phage-data4$pAbi2
data  = rbind(data1,data2,data3,data4)

#Abi2 vs phage count
plot_Phagecount <- ggplot(data, aes(x=Belongs, y=Phage, fill=Belongs)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Bacteria"="#8ECDDD", "Phage"="#FFBB5C", "zNone" = "white", "aBoth" = "#AE445A")) +
  # geom_jitter(color="gray", size=2, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),legend.position="none") +
  ylab("Count of phage without Abi2") +
  xlab("Abi2") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons =  split(t(combn(levels(as.factor(data$Belongs)), 2)), seq(nrow(t(combn(levels(as.factor(data$Belongs)), 2))))), #list(c("Bacteria","Phage")), 
              map_signif_level = F, textsize=4, step_increase = .06)
plot_Phagecount
ggsave("250807_Abi2_phagecount_boxplot_v1.pdf", plot_Phagecount , width=5, height=5)

##-------------------------------------------------------------------------------------------------------------------------
## Fig. S6B
data1 = filter(df, bGabija>=1&pGabija==0)
data1$Belongs = replicate(dim(data1)[1], "Bacteria")
data2 = filter(df, pGabija>=1&bGabija==0)
data2$Belongs = replicate(dim(data2)[1], "Phage")
data2$Phage = data2$Phage-data2$pGabija
data3 = filter(df, bGabija==0&pGabija==0)
data3$Belongs = replicate(dim(data3)[1], "zNone")
data4 = filter(df, bGabija>=1&pGabija>=1)
data4$Belongs = replicate(dim(data4)[1], "aBoth")
data4$Phage = data4$Phage-data4$pGabija
data  = rbind(data1,data2,data3,data4)

#Gabija vs phage count
plot_Phagecount <- ggplot(data, aes(x=Belongs, y=Phage, fill=Belongs)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Bacteria"="#8ECDDD", "Phage"="#FFBB5C", "zNone" = "white", "aBoth" = "#AE445A")) +
  # geom_jitter(color="gray", size=2, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),legend.position="none") +
  ylab("Count of phage without Gabija") +
  xlab("Gabija") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons =  split(t(combn(levels(as.factor(data$Belongs)), 2)), seq(nrow(t(combn(levels(as.factor(data$Belongs)), 2))))), #list(c("Bacteria","Phage")), 
              map_signif_level = F, textsize=4, step_increase = .06)
plot_Phagecount
ggsave("250807_Gabija_phagecount_boxplot_v1.pdf", plot_Phagecount , width=5, height=5)

