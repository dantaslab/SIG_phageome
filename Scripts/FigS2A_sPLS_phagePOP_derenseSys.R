## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '90%')


## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(5249) # for reproducibility, remove for normal use

require(doParallel)
cl = makeCluster(6)
registerDoParallel(cl)

## -------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library("compositions")

setwd('~/Google Drive/WashU/Phageome/SIG/Manuscript/Figures/Figure2/Fig2F')

# data_all = read_excel("230831_phage_info.xlsx")
# data = as.data.frame.matrix(table(data_all$Isolate, data_all$Population))
# write.csv(data, file = "230904_isolate_phagePOP_matrix.csv")

X_all   = read_excel("230904_isolate_defense.xlsx")
X_range = X_all[,c(1:31)] #filter(X_all, Days_before_NEC_onset >= 0) #filter(X_all, Days_before_NEC_onset == 0) #filter(X_all,NEC_onset2 == "Late")
X_data  = subset(X_range, select = -c(IsolateID))

Y_all   = read_excel("230904_isolate_phagePOP_matrix.xlsx")
Y_range = Y_all
#Y_range = filter(Y_all, Days_before_NEC_onset >= 0) #filter(Y_all, Days_before_NEC_onset == 0) #filter(Y_all,NEC_onset2 == "Late")
Y_data <- subset(Y_range, select = -c(IsolateID))

# low.count.removal <- function(
#     data, # OTU count df of size n (sample) x p (OTU)
#     percent=0 # cutoff chosen
# ) 
# {
#   keep.otu = which(colSums(data!=0)*100/dim(data)[1] > percent)
#   data.filter = data[,keep.otu]
#   return(list(data.filter = data.filter, keep.otu = keep.otu))
# }
# X_data.filter <- low.count.removal(X_data, percent=1)
# X.filter <- X_data.filter$data.filter
X <- logratio.transfo(as.matrix(X_data), logratio = 'CLR', offset = 1) 
# Y_data.filter <- low.count.removal(Y_data, percent=1)
# Y.filter <- Y_data.filter$data.filter
Y <- logratio.transfo(as.matrix(Y_data), logratio = 'CLR', offset = 1) 

# X = scale(X)
# Y = scale(Y)

## ---- out.width = "49%", fig.show = "hold", fig.cap = "FIGURE 1: Barplot of the variance each principal component explains of the liver toxicity data."----
# pca.gene <- pca(X, ncomp = 5, center = TRUE, scale = FALSE)
# pca.clinic <- pca(Y, ncomp = 5, center = TRUE, scale = FALSE)
# 
# plot(pca.gene)
# plot(pca.clinic)
# 
# 
# ## ---- out.width = "49%", fig.show = "hold", fig.cap = "FIGURE 2: Preliminary (unsupervised) analysis with PCA on the liver toxicity data"----
# plotIndiv(pca.gene, comp = c(1, 2),
#           group = X_range$type,
#           ind.names = X_range$type,
#           legend = TRUE, title = 'Liver gene, PCA comp 1 - 2')
# 
# plotIndiv(pca.clinic, comp = c(1, 2),
#           group =  Y_range$Group,
#           ind.names = Y_range$DOL,
#           legend = TRUE, title = 'Liver clinic, PCA comp 1 - 2')


## -------------------------------------------------------------------------------------------------------------------
spls.liver <- spls(X = X, Y = Y, ncomp = 3, mode = 'regression', scale = FALSE) # form basic model


## ---- fig.cap = "FIGURE 3: Tuning the number of components in PLS on the liver toxicity data. For each component, the repeated cross-validation (5 × 10−fold CV) $Q^2$ score is shown. Horizontal line depicts $Q^2$ = 0.0975. The bars represent the variation of these values across the repeated folds."----
perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
                         folds = 10, nrepeat = 5) # repeated CV tuning of component count
 
plot(perf.spls.liver, criterion = 'Q2.total')


## ---- fig.cap = "FIGURE 4: Tuning plot for sPLS2."------------------------------------------------------------------
set.seed(222)
list.keepX <- c(seq(1, 25, 5))  # set range of test values for number of variables to use from X dataframe
# c(seq(2, 26, 3))
list.keepY <- c(seq(1, 47, 5)) ##1:c(seq(3, 45, 5))  5:c(seq(1:11)) # set range of test values for number of variables to use from Y dataframe

tune.spls.liver <- tune.spls(X, Y, ncomp = 2,
                              test.keepX = list.keepX,
                              test.keepY = list.keepY,
                              nrepeat = 1, folds = 10, # use 10 folds
                              mode = 'regression', measure = 'cor') # use the correlation measure for tuning
plot(tune.spls.liver)


## -------------------------------------------------------------------------------------------------------------------
tune.spls.liver$choice.keepX
tune.spls.liver$choice.keepY


## -------------------------------------------------------------------------------------------------------------------
optimal.keepX <- tune.spls.liver$choice.keepX # extract optimal number of variables for X dataframe
optimal.keepY <- tune.spls.liver$choice.keepY # extract optimal number of variables for Y dataframe
optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components



## -------------------------------------------------------------------------------------------------------------------
final.spls.liver <- spls(X, Y, ncomp = optimal.ncomp, # use all tuned values from above
                    keepX = optimal.keepX,
                    keepY = optimal.keepY,
                    mode = "regression") # explanitory approach being used, hence use regression mode


## ---- out.width = "49%", fig.show = "hold", fig.cap = "FIGURE 5: Sample plot for sPLS2 performed on the liver.toxicity data. Samples are projected into the space spanned by the components associated to each data set (or block)."----
# plotIndiv(final.spls.liver, ind.names = FALSE, 
#          rep.space = "X-variate", # plot in X-variate subspace
#          group = X_range$Group, # colour by time group
#          #pch = as.factor(X_range$DOL), # select symbol by dose group
#          col.per.group = color.mixo(1:3), 
#          legend = TRUE, legend.title = 'NEC status', legend.title.pch = 'DOL')
# 
# plotIndiv(final.spls.liver, ind.names = FALSE,
#          rep.space = "Y-variate", # plot in Y-variate subspace
#          group = Y_range$Group, # colour by time group
#          #pch = as.factor(liver.toxicity$treatment$Dose.Group), # select symbol by dose group
#          col.per.group = color.mixo(1:3), 
#          legend = TRUE, legend.title = 'NEC status', legend.title.pch = 'DOL')
# 
# 
# ## ---- fig.cap = "FIGURE 6: Sample plot for sPLS2 performed on the liver.toxicity data. Samples are projected into the space spanned by the averaged components of both datasets."----
# plotIndiv(final.spls.liver, ind.names = FALSE, 
#          rep.space = "XY-variate", # plot in averaged subspace
#          group = X_range$Group, # colour by time group
#          #pch = as.factor(liver.toxicity$treatment$Dose.Group), # select symbol by dose group
#          col.per.group = color.mixo(1:3), 
#          legend = TRUE, legend.title = 'NEC status', legend.title.pch = 'DOL')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
# col.tox <- color.mixo(as.numeric(as.factor(X_range$Group))) # create set of colours
# plotIndiv(final.spls.liver, ind.names = FALSE,
#           rep.space = "XY-variate", # plot in averaged subspace
#           axes.box = "both", col = col.tox, style = '3d')


## ---- fig.cap = "FIGURE 8:  Arrow plot from the sPLS2 performed on the liver.toxicity data. The start of the arrow indicates the location of a given sample in the space spanned by the components associated to the gene data set, and the tip of the arrow the location of that same sample in the space spanned by the components associated to the clinical data set."----
# plotArrow(final.spls.liver, ind.names = FALSE,
#           group = X_range$Group, # colour by time group
#           col.per.group = color.mixo(1:4),
#           legend.title = 'Time.Group')


## ---- fig.cap = "FIGURE 9:  Stability of variable selection from the sPLS on the Liver Toxicity gene expression data. The barplot represents the frequency of selection across repeated CV folds for each selected gene for component 1 (a) and 2 (b)."----
# form new perf() object which utilises the final model
perf.spls.liver <- perf(final.spls.liver, 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)

# plot the stability of each feature for the first two components, 'h' type refers to histogram
par(mfrow=c(1,2)) 
plot(perf.spls.liver$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls.liver$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))


## ---- fig.cap = "FIGURE 10:  Correlation circle plot from the sPLS2 performed on the liver.toxicity data. This plot should be interpreted in relation to Figure 5 to better understand how the expression levels of these molecules may characterise specific sample groups."----
plotVar(final.spls.liver, cex = c(3,4), var.names = c(TRUE, TRUE))


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
color.edge <- color.GreenRed(30)  # set the colours of the connecting lines

# X11() # To open a new window for Rstudio
network(final.spls.liver, comp = 1:2,
        cutoff = 0.05, # only show connections with a correlation above 0.7
        #shape.node = c("rectangle", "circle"),
        color.node = c("orange","lightblue"), # "orange","lightblue" ;"cyan", "pink"
        color.edge = color.edge,
        lwd.edge = 3,
        save = 'pdf', # save as a png to the current working directory
        name.save = '250710_sPLS_phagePOP_defense_Network_0_cut0.05_v2')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
#cim(final.spls.liver, comp = 1:2, xlab = "Bactiera", ylab = "Phage host", row.cex = 2)
cim(final.spls.liver, comp = 1:2, xlab = "Phage population", ylab = "Phage defense system", 
    margins = c(6, 6), transpose = T,
    save = 'pdf', # save as a png to the current working directory
    name.save = '250710_sPLS_phagePOP_defense_heatmap_0_v2')

A = cim(final.spls.liver, comp = 1:2, xlab = "Phage population", ylab = "Phage defense system", 
        margins = c(6, 6), transpose = T)
min(A$mat)
max(A$mat)

B = A$mat
stopCluster(cl)
registerDoSEQ()

