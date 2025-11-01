##Recursive feature elimination using caret
library(ggplot2)
library(caret)
library(vegan)
library(dendextend)
library(gplots)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library("compositions")
library(readxl)
set.seed(125)

gene_data_all= read_excel("Data_S4_v1.xlsx", sheet = "defense_sys_sum")
# gene_data_all$type = factor(gene_data_all$type, levels = c('Abi2', 'AbiD', 'AbiJ', 'AbiL', 'AbiN', 'AbiP2', 'AbiR', 
#                                                            'AbiU', 'Pycsar', 'SEFIR', 'Stk2', 'Thoeris',
#                                                           'Cas', 'Gabija', 'Nhi', 'RloC', 'RM', 'CBASS', 'Lamassu-Fam', 
#                                                           'Retron', 'Dodola', 'DRT', 'Eleos',
#                                                           'Kiwa', 'PD-Lambda-1', 'PD-T4-9', 'PD-T7-2', 'RosmerTA', 'ShosTA', 'SoFIC'))
gene_data = as.matrix(table(gene_data_all$isolate, gene_data_all$type))

#calculate correlation matrix
correlationMatrix <- cor(gene_data)

# #summarize the correlation matrix
# highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0) # 0.75
# highlyCorrelated = sort(highlyCorrelated)
# 
# #Remove highly correlated variables from dataset. Left with 153 uncorrelated features
# reduced_data     = gene_data[,-c(highlyCorrelated)]
# correlate_data   = gene_data[,c(highlyCorrelated)]
# correlate_Matrix = correlationMatrix[c(highlyCorrelated),c(highlyCorrelated)]
# #write.csv(as.data.frame(as.table(correlate_Matrix)),"high_correlate_data_scoary.csv")

require(corrplot)
correlationMatrix[correlationMatrix==1]=0
# r = rbind(c('Abi2','AbiD','AbiJ','AbiK','AbiL','AbiN','AbiP2','AbiR','AbiU','Pycsar','SEFIR','Stk2','Thoeris'),
#           c('Cas','Gabija','Nhi','RloC','RM'),
#           c('CBASS','Lamassu-Fam','Retron'),
#           c('Dodola','DRT','Eleos','Kiwa','PD-Lambda-1','PD-T4-9','PD-T7-2','RosmerTA','ShosTA','SoFIC'))
testRes = cor.mtest(gene_data, conf.level = 0.95)
corrplot(correlationMatrix, tl.cex = 0.9, order = 'original', type = 'upper', method = 'color', diag = FALSE, is.corr = FALSE,
         p.mat = testRes$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', col.lim = c(-0.35,0.55)) %>%
  corrRect(name = c('Abi2', 'Cas', 'CBASS', 'Dodola', 'SoFIC'))

min(correlationMatrix)
max(correlationMatrix)
