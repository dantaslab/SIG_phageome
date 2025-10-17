library("compositions")
library("readxl")
library("xlsx")
library("dplyr")
library("vegan")
library("Matrix")
library("phyloseq")
library("data.table")
library("reshape2")
library("ggplot2")
library(data.table, warn.conflicts=FALSE)

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
####################### https://rpubs.com/lconteville/714853

set.seed(129)
df <- read.csv("isolate_vOTU_matrix.csv", header = T)
df = filter(df, !grepl("SF_0046|SF_0070|SF_0126|SF_0255|SF_0256|SF_0257|SF_0258|SF_0260", Isolate)) # S. delphini
k <- rowSums(df[,3:166]) > 0
df <- df[k,]
df_all <- select(df, select = -c(Isolate, Cohort))

#SM_df = filter(SM_df, NEC_onset2=="Late")
#df = filter(SM_df, NEC_status=="case")
taxmat <- as.data.frame(colnames(df_all[,1:164]))
colnames(taxmat)[1] <- "Phage_Population"

RowName <- df$Isolate
metadata <- as.data.frame(df$Cohort)
rownames(metadata) = RowName
colnames(metadata)[1] <- "Cohort"
#df <- select(df, -SampleID, -NEC_onset2, -NEC_status, -DOL, -Cohort)
rownames(df) = RowName
t_df = t(df_all)

setdiff(rownames(t_df),taxmat$Phage_Population)

rownames(taxmat) <- taxmat$Phage_Population
taxmat =  as.matrix(taxmat)

OTU = otu_table(t_df, taxa_are_rows = TRUE)
colnames(OTU) = df$Isolate
sampledata = sample_data(metadata)
TAX = tax_table(taxmat)
physeq1 = phyloseq(OTU, TAX, sampledata)

physeq1
head(otu_table(physeq1)[,1:6])

######################

abrel_jaccard <- phyloseq::distance(physeq1, method = "jaccard")
abrel_jaccard <- as.matrix(abrel_jaccard)
head(abrel_jaccard)[,1:6]

case_JAC <- melt(abrel_jaccard)
case_JAC_pair <- case_JAC[case_JAC$Var1!=case_JAC$Var2,]
colnames(case_JAC_pair) = c("Isolate1","Isolate2","Jaccard")

# write.csv(case_JAC_pair,"Jaccard_pairs.csv")

isolate_info <- read_excel("230831_isolate_info.xlsx")
isolate_info_sub = select(isolate_info, Isolate, Source, Type, Group, Bodysite, MLST, Household, Intervention)
isolate_info_sub1 = isolate_info_sub
colnames(isolate_info_sub1) = c("Isolate1", "Source1", "Type1", "Group1", "Bodysite1", "MLST1", "Household1", "Intervention1")
JAC_pair_meta1 = merge(case_JAC_pair, isolate_info_sub1, by = "Isolate1")
isolate_info_sub2 = isolate_info_sub
colnames(isolate_info_sub2) = c("Isolate2", "Source2", "Type2", "Group2", "Bodysite2", "MLST2", "Household2", "Intervention2")
JAC_pair_meta = merge(JAC_pair_meta1, isolate_info_sub2, by = "Isolate2")

JAC_pair_meta$Source_group <- ifelse(JAC_pair_meta$Source1==JAC_pair_meta$Source2, JAC_pair_meta$Source1, 'NA')
JAC_pair_meta$Type_group <- ifelse(JAC_pair_meta$Type1==JAC_pair_meta$Type2, JAC_pair_meta$Type1, 'NA')
JAC_pair_meta$Group_group <- ifelse(JAC_pair_meta$Group1==JAC_pair_meta$Group2, JAC_pair_meta$Group1, 'NA')
JAC_pair_meta$Bodysite_group <- ifelse(JAC_pair_meta$Bodysite1==JAC_pair_meta$Bodysite2, JAC_pair_meta$Bodysite1, 'NA')
JAC_pair_meta$MLST_group <- ifelse(JAC_pair_meta$MLST1==JAC_pair_meta$MLST2, 'Within', 'Between')
JAC_pair_meta$Household_group <- ifelse((JAC_pair_meta$Household1=='NA' | JAC_pair_meta$Household2=='NA'), 'NA',  ifelse(JAC_pair_meta$Household1==JAC_pair_meta$Household2, 'Within', 'Between'))
JAC_pair_meta$Intervention_group <- ifelse((JAC_pair_meta$Intervention1=='NA' | JAC_pair_meta$Intervention2=='NA'), 'NA',  ifelse(JAC_pair_meta$Intervention1==JAC_pair_meta$Intervention2, JAC_pair_meta$Intervention1, 'Between'))

write.csv(JAC_pair_meta,"250712_Jaccard_phyloseq.csv")
