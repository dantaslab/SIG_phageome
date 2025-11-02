library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library(ggpmisc)
library('caret')
set.seed(126)

data = read_excel("Data_S5_v2.xlsx", sheet = "Spacer_prophage_match")
data = filter(data, Isolates>=5)

##--------------------------------------------------------------------------------------------------------------------------------------------
## Fig. 3E
plot1 <- ggplot(data, aes(x=`CasIII%`, y=`CasII%`,color = Cas_kind_50)) +
  geom_count(alpha=0.8) +
  # scale_size(breaks = c(1,10,100,200,300)) +
  scale_size_area(breaks = c(1,5,20,40)) +
  scale_color_manual(values = c("CasII"="#618264", "Both"="#FFBB5C","CasIII"="#9A3B3B")) +
  geom_hline(yintercept=90, linetype="dashed", color = "gray", size=0.5) +
  geom_vline(xintercept=90, linetype="dashed", color = "gray", size=0.5) +
  geom_hline(yintercept=50, linetype="dashed", color = "gray", size=0.5) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12),legend.position=c(0.3,0.25)) +
  theme(legend.direction = "vertical", legend.box = "horizontal") +
  xlab("Spacer pop in isolates with CasIII %") +
  ylab("Spacer pop in isolates with CasII %") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed((ratio=1)) + ylim(0,100)
plot1
ggsave("250722_scatter_spacerPOP_Cas_5_kind_v1.pdf", plot1, width=4, height=4)

##--------------------------------------------------------------------------------------------------------------------------------------------
## Fig. 3F
plot_phagePOP_count <- ggplot(data, aes(x=Cas_kind_50, y=PhagePOP_count, fill=Cas_kind_50)) +
  geom_boxplot() +
  scale_fill_manual(values = c("CasII"="#618264", "Both"="#FFBB5C","CasIII"="#9A3B3B")) +
  # geom_jitter(color="gray", size=0.5, alpha=0.5) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 2)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Count of targeted vOTU") +
  xlab("Spacer Cas kind") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Cas_kind_50)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cas_kind_50)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
plot_phagePOP_count
ggsave("250722_Box_CasSpacer_matched_phagePop_count_5.pdf", plot_phagePOP_count, width=4, height=5)


