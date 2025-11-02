library(ggstatsplot)
library(readxl)
library(dplyr)
library(ggplot2)

target_data = read_excel("Data_S5_v2.xlsx", sheet = "CRISPRTarget")
spacer_info = read_excel("Data_S5_v2.xlsx", sheet = "Spacer_prophage_match")

##################### Spacer Target #####################

data = merge(target_data, spacer_info, by.x = "Spacer_cluster", by.y = "Pop")
data$Protospacer_description_refine[data$Protospacer_description_refine=="Virus"]="Unknown"
data$Protospacer_description_refine[data$Protospacer_description_refine=="UncultivatedPhage"] = "Phage"
data$Protospacer_description_refine[data$Protospacer_seq_id=="None"] = "zNo target" 
data$Protospacer_description_refine[data$Protospacer_description_refine=="Plasmid"] = "qPlasmid" 
                                    
##--------------------------------------------------------------------------------------------------------
## Fig. S4A
pie_spacer_target1 <- ggpiestats(
  data = data,
  x = Protospacer_description_refine,
  legend.title ="Spacer target",
  xlab = "Spacer population source",
  ylab = "Percentage",
  package = "awtools",
  palette = "a_palette") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="right")
pie_spacer_target1
ggsave("250722_pie_spacer_pop_target.pdf", pie_spacer_target1, width=5, height=5)

##--------------------------------------------------------------------------------------------------------
## Fig. S4B
pie_spacer_target2 <- ggpiestats(
  data = data,
  x = Protospacer_description_refine,
  y = Match_phages,
  legend.title ="Spacer target",
  xlab = "Spacer population source",
  ylab = "Percentage",
  package = "awtools",
  palette = "a_palette") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="none")
pie_spacer_target2
ggsave("250722_pie_spacer_pop_target_HIT.pdf", pie_spacer_target2, width=8, height=5)

##--------------------------------------------------------------------------------------------------------
## Fig. 3G
# data_cas = data
# data_cas$Cas_kind_50[data$Isolates<5] = "NA"
data_cas = filter(data, Isolates>=5 & Cas_kind_50!="Both")

pie_spacer_target3 <- ggpiestats(
  data = data_cas,
  x = Protospacer_description_refine,
  y = Cas_kind_50,
  legend.title ="Spacer target",
  xlab = "Spacer population source",
  ylab = "Percentage",
  package = "awtools",
  palette = "a_palette") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="none")
pie_spacer_target3
ggsave("250722_pie_spacer_pop_target_Cas.pdf", pie_spacer_target3, width=8, height=5)


##################### Spacer Bacterial Host #####################

data = merge(target_data, spacer_info, by.x = "Spacer_cluster", by.y = "Pop")
data = filter(data,Protospacer_description_refine!="None")
# Refine at genus level
data$Host[data$Host=="Staphylococcaceae"]="Unknown"
data$Host[data$Host=="Leptotrichiaceae"]="Unknown"
data$Host[data$Host!="Unknown" & data$Host!="Multiple" & data$Host!="Hymenobacter" & data$Host!="Staphylococcus" & data$Host!="Klebsiella" & data$Host!="Clostridium"] = "Others"

##--------------------------------------------------------------------------------------------------------
## Fig. S4C
pie_spacer_target_host1 <- ggpiestats(
  data = data,
  x = Host,
  legend.title ="Spacer target_host",
  xlab = "Spacer population source",
  ylab = "Percentage",
  perc.k=1L,
  package = "awtools",
  palette = "a_palette") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="right")+
  scale_fill_brewer(palette="RdBu", direction = -1)
pie_spacer_target_host1
ggsave("250809_pie_spacer_pop_target_host.pdf", pie_spacer_target_host1, width=5, height=5)


##--------------------------------------------------------------------------------------------------------
## Fig. S4D
# data_cas = data
# data_cas$Cas_kind_50[data$Isolates<5] = "NA"
data_cas = filter(data, Isolates>=5 & Cas_kind_50!="Both")

pie_spacer_target_host3 <- ggpiestats(
  data = data_cas,
  x = Host,
  y = Cas_kind_50,
  legend.title ="Spacer target_host",
  xlab = "Spacer population source",
  ylab = "Percentage",
  perc.k=1L,
  package = "awtools",
  palette = "a_palette") +
  labs(caption = NULL) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16),legend.position="none") +
  scale_fill_brewer(palette="RdBu", direction = -1)
pie_spacer_target_host3
ggsave("250809_pie_spacer_pop_target_host_Cas.pdf", pie_spacer_target_host3, width=8, height=5)



