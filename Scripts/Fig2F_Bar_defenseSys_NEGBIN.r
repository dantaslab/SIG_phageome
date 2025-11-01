library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
set.seed(126)

setwd("~/Google Drive/WashU/Phageome/SIG/Manuscript/Figures/Figure2/Fig2C")
data_group = read.table("230903_DefenseCount_Group_NEGBIN_0.01/significant_results.tsv", header = TRUE)
data_source = read.table("230903_DefenseCount_Source_NEGBIN_0.01/significant_results.tsv", header = TRUE)
data_type = read.table("230903_DefenseCount_Type_NEGBIN_0.01/significant_results.tsv", header = TRUE)
data = rbind(data_group, data_source, data_type)
data$Sig_associations = -log10(data$qval)*sign(data$coef)

#-------------------------------------------------------------------------------------------------------------------------------------------

plot_group <- ggplot(data, aes(x=value, y=coef, fill = value)) +
  geom_bar(stat = "identity", width=0.7, color="black") + coord_flip() +
  scale_fill_manual(values = c("Diagnostic"="#EF8A62", "AI"="#F1948A", "HI"="#CB4335", "Human"="#435334")) +
  #theme_ipsum() +
  geom_hline(yintercept=0) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  xlab("Phage defense system") +
  ylab("Coefficient") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
# ylim(0,12.5)
plot_group 
ggsave("250710_Bar_defenseSys_NEGBIN.pdf", plot_group, width=5, height=5)


plot_group2 <- ggplot(data, aes(x=value, y=Sig_associations, fill = value)) +
  geom_bar(stat = "identity", width=0.7, color="black") + coord_flip() +
  scale_fill_manual(values = c("Diagnostic"="#EF8A62", "AI"="#F1948A", "HI"="#CB4335", "Human"="#435334")) +
  #theme_ipsum() +
  geom_hline(yintercept=0) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  xlab("Phage defense system") +
  ylab("Significant associations to household") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(limits = c('AI','HI', 'Human', 'Diagnostic'))
# ylim(0,12.5)
plot_group2
ggsave("250710_Bar_defenseSys_NEGBIN_v2.pdf", plot_group2, width=5, height=5)

