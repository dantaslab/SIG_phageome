library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')
set.seed(125)

data = read_excel("Data_S3_v2.xlsx", sheet = "All_AMG")
colnames(data)  = as.matrix(data[1,])
data = data[-1,]
df_isolateAMG = as.data.frame.matrix(table(data$Isolate, data$query_name))
df_isolateAMG$total = rowSums(df_isolateAMG)

data_isolate = read_excel("230831_isolate_info.xlsx")
metadata_df <- as.data.frame(data_isolate[, c(1:10)])
row.names(metadata_df) = metadata_df$Isolate

data = merge(df_isolateAMG, metadata_df, by = "row.names")

#--------------------------------------------------------------------------------------------------------------------------------------------
# stat_pvalue <- data %>% 
#   rstatix::wilcox_test(total ~ Type) %>%
#   filter(p < 0.05) %>% 
#   rstatix::add_significance("p") %>% 
#   rstatix::add_y_position() %>% 
#   mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

fun_mean <- function(x){
  return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}

fun_median <- function(x){
  return(data.frame(y=median(x),label=median(x,na.rm=T)))}

# data$Type[data$Type=="Colonizing"] = "Household"
# data$Type[data$Type=="Environment"] = "Household"
plot_type <- ggplot(data, aes(x=Type, y=total, fill=Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Colonizing"="#67A9CF", "Diagnostic"="#EF8A62","Environment"="#B3B6B7")) +
  # scale_fill_manual(values = c("Diagnostic"="#EF8A62","Household"="#48C9B0")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Type)), 2)), seq(nrow(t(combn(levels(as.factor(data$Type)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
  # stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  # stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7)
  # ylim(0,12.5)
  #stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_type 
ggsave("250804_Box_AMGcount_type.pdf", plot_type, width=5, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------
data_group = filter(data, Group!="Household")
plot_group <- ggplot(data_group, aes(x=Group, y=total, fill=Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("AI"="#F1948A", "HI"="#CB4335","Household"="#48C9B0")) +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data_group$Group)), 2)), seq(nrow(t(combn(levels(as.factor(data_group$Group)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
# ylim(0,12.5)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_group 
ggsave("250804_Box_AMGcount_group.pdf", plot_group, width=5, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------
data = filter(data, Source == "Human")
data = filter(data, Bodysite == "nares"|Bodysite == "respiratory"|Bodysite == "tissue"|Bodysite == "wound")

plot_bodysite <- ggplot(data, aes(x=Bodysite, y=total, fill=Bodysite)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=1, option = "E") +
  # geom_jitter(color="gray", size=0.8, alpha=1) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size = 4)+
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = split(t(combn(levels(as.factor(data$Bodysite)), 2)), seq(nrow(t(combn(levels(as.factor(data$Bodysite)), 2))))), #stat_pvalue$groups, #list(c("G", "S")), split(t(combn(levels(as.factor(data$Cluster)), 2)), seq(nrow(t(combn(levels(as.factor(data$Cluster)), 2))))),
              map_signif_level = F, textsize=5,step_increase = .1)
  # scale_y_continuous(breaks = seq(0, 11, by = 2))
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_bodysite
ggsave("250804_Box_AMGcount_bodysite.pdf", plot_bodysite, width=5, height=5)

#--------------------------------------------------------------------------------------------------------------------------------------------
data = filter(data, Source == "Human")
data = filter(data, Age_group != "NA")

stat_pvalue <- data %>%
  rstatix::wilcox_test(total ~ Age_group) %>%
  filter(p < 0.05) %>%
  rstatix::add_significance("p") %>%
  rstatix::add_y_position() %>%
  mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

plot_age_group <- ggplot(data, aes(x=Age_group, y=total, fill=Age_group)) +
  geom_boxplot() +
  scale_fill_brewer(palette="RdBu") +
  geom_jitter(color="gray", size=0.8, alpha=1) +
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = stat_pvalue$groups, #split(t(combn(levels(as.factor(data$Age_group)), 2)), seq(nrow(t(combn(levels(as.factor(data$Age_group)), 2))))), #stat_pvalue$groups, 
              map_signif_level = F, textsize=5,step_increase = .1)
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_age_group
ggsave("231005_Box_AMGcount_age.pdf", plot_age_group, width=5, height=5)

library(ggpmisc)
plot_age <- ggplot(data, aes(x=as.numeric(Age), y=total)) +
  geom_smooth(method="loess", color="dark gray",alpha=0.5) +
  geom_point(size=1,alpha=1, aes(color=Age_group)) + #, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_color_brewer(palette="RdBu") +
  # stat_poly_line() +
  # stat_poly_eq(use_label(c("adj.R2", "p"))) +
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position=c(0.15,0.85)) +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_continuous(limits = c(0, 90), breaks = c(0,10,20,30,40,50,60,70,80,90))
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_age
ggsave("231005_Scatter_AMGcount_age.pdf", plot_age, width=5, height=5)

plot_age2 <- ggplot(data, aes(x=as.numeric(Age), y=total,color=Age_group2)) +
  # geom_smooth(method="loess") +
  geom_point(size=1,alpha=1) + #, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_color_brewer(palette="Set1") +
  stat_poly_line() +
  stat_poly_eq(use_label(c("adj.R2", "p"))) +
  #theme_ipsum() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("AMG count") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_continuous(limits = c(0, 90), breaks = c(0,10,20,30,40,50,60,70,80,90))
#stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 1.5)
plot_age2
ggsave("231005_Scatter_AMGcount_age_v2.pdf", plot_age2, width=5, height=5)

