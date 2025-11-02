library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library(ggpmisc)
library('caret')
set.seed(126)

data = read_excel("250722_phage_spacer_hit.xlsx", sheet = "PhagePop_member")

#--------------------------------------------------------------------------------------------------------------------------------------------
plot2 <- ggplot(data, aes(x=vOTU_Members, y=matched_spacerCluster_count)) +
  # geom_point(aes(size=Group_size, color = Match_with_spacers),alpha=0.35) + 
  geom_count(aes(color = Match_with_spacers), alpha=0.5) +
  scale_size_area(breaks = c(1,5,9)) +
  # scale_size(breaks = c(1,5,9)) +
  stat_poly_line(color = "black") +
  stat_poly_eq(use_label(c("eq","adj.R2", "p")), label.x = 0.1) +
  # scale_color_manual(values = c("No"="#CAEDFF", "Yes"="#ff6666")) +
  scale_color_manual(values = c("No"="#ec8100", "Yes"="#4085c6")) +
  # stat_poly_eq(use_label(c("eq", "adj.R2", "f", "p", "n"))) +
  # scale_color_manual(values = c("#293B5F")) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12),legend.position=c(0.2,0.65)) + #c(0.2,0.6)
  xlab("Phage population size") +
  ylab("Matched spacer cluters") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + scale_x_log10()
plot2
ggsave("250722_scatter_phagePOP_spacerPOP_v2.pdf", plot2, width=5, height=5)
