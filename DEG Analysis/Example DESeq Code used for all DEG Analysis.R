## Ligand/Receptor Project 
## Rounds 3, 3.5 and 4 Bulk RNA-Seq Analysis 

## 7-10-21 

## Anna Voss | Sloan Lab 

##########################################################
#DESeq Ligands vs. Control d60-90 1363.1
########################################################## 
#Adapted from: 
#https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html#example-3-two-conditions-two-genotypes-with-an-interaction-term
#Build SummarizedExperiment from count matrix

setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis")
cts <- read.csv("R3_and_4_feature_counts.csv",header=TRUE)
rownames(cts) <- make.names(cts[,1], unique=TRUE)
gene_names<-cts[,1]
cts <- cts[,-1]
cts <- cts[,-1]
View(cts)

round3_5_M6_10 <- cts[,c(8:11,15:17)]
View(round3_5_M6_10)
gene_names <- gene_names

#Set Up Col Data 
condition <- c("control", "control", "control","control","M6-10","M6-10", "M6-10")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)

##Setup DESeq2
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = round3_5_M6_10,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)
dds$condition = relevel(dds$condition, "control") #Set hypoxia as the default
dds$condition
resultsNames(dds)

#Run Comparison
res <- results(dds, contrast=c("condition","M6-10", "control"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "R3_5_M6-10.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

# Make a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = gene_names,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8),
                pCutoff = 0.001,
                title = "Control (left) vs. M6-10 (right) d60-90 (R3.5)")

##########################################################
#DESeq R4 Ligands vs. Control d90-120 1363.1
########################################################## 
#Adapted from: 
#https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html#example-3-two-conditions-two-genotypes-with-an-interaction-term
#Build SummarizedExperiment from count matrix

cts <- read.csv("R3_and_4_feature_counts.csv",header=TRUE)
rownames(cts) <- make.names(cts[,1], unique=TRUE)
gene_names<-cts[,1]
cts <- cts[,-1]
cts <- cts[,-1]

round4_M6_10 <- cts[,c(21:24,28:30)]
View(round4_M6_10)
gene_names <- gene_names

#Set Up Col Data 
condition <- c("control", "control", "control","control","M6-10","M6-10", "M6-10")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)

##Setup DESeq2
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = round4_M6_10,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)
dds$condition = relevel(dds$condition, "control") #Set hypoxia as the default
dds$condition
resultsNames(dds)

#Run Comparison
res <- results(dds, contrast=c("condition","M6-10", "control"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "R4_M6-10.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

# Make a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = gene_names,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8),
                pCutoff = 0.001,
                title = "Control (left) vs. M6-10 (right) d90-120 (R4)")


##########################################################
##Unsupervised Clustering Ligands d60-90 1363.1 
##########################################################
round3_5 <- cts[,c(21:33)]

#Variance 
V <- apply(round3_5, 1, var, na.rm = TRUE)
View(cts)

#Sort variance in decreasing order and select top 100, 1000, and 5000 genes 
selectedGenes100 <- names(V[order(V, decreasing = TRUE)][1:100])
selectedGenes1k <- names(V[order(V, decreasing = TRUE)][1:1000])
selectedGenes5k <- names(V[order(V, decreasing = TRUE)][1:5000])

#Plot that sucker 
library(pheatmap)
genes100 <- pheatmap(round3_5[selectedGenes100,], scale = 'row', show_rownames = FALSE)
genes1k <- pheatmap(round3_5[selectedGenes1k,], scale = 'row', show_rownames = FALSE)
genes5k <- pheatmap(round3_5[selectedGenes5k,], scale = 'row', show_rownames = FALSE)

##########################################################
##PCA d60-90 1363.1 
########################################################## 
library(stats)
library(ggplot2)
library(ggfortify) #needed to let ggplot2 know about PCA data structure.

#transpose the matrix 
transpose100 <- t(round3_5[selectedGenes100,])
transpose1k <- t(round3_5[selectedGenes1k,])
transpose5k <- t(round3_5[selectedGenes5k,])

# transform the counts to log2 scale 
transpose100 <- log2(transpose100 + 1)
transpose1k <- log2(transpose1k + 1)
transpose5k <- log2(transpose5k + 1)

#compute PCA 
pcaResults100 <- prcomp(transpose100)
pcaResults1k <- prcomp(transpose1k)
pcaResults5k <- prcomp(transpose5k)

#Add Metadata
condition <- c("control", "control", "control", "control", "M1-5", "M1-5","M1-5",
            "M6-10","M6-10","M6-10","M1-10","M1-10","M1-10")

transpose100 <- cbind(transpose100, condition)
transpose1k <- cbind(transpose1k, condition)
transpose5k <- cbind(transpose5k, condition)

#plot PCA results making use of ggplot2's autoplot function
autoplot(pcaResults100, data = transpose100, colour = "condition")
autoplot(pcaResults1k, data = transpose1k, colour = "condition")
autoplot(pcaResults5k, data = transpose5k, colour = "condition")


##########################################################
##PCA d90-120 1363.1 
########################################################## 
library(stats)
library(ggplot2)
library(ggfortify) #needed to let ggplot2 know about PCA data structure.

#transpose the matrix 
transpose100 <- t(round4[selectedGenes100,])
transpose1k <- t(round4[selectedGenes1k,])
transpose5k <- t(round4[selectedGenes5k,])

# transform the counts to log2 scale 
transpose100 <- log2(transpose100 + 1)
transpose1k <- log2(transpose1k + 1)
transpose5k <- log2(transpose5k + 1)

#compute PCA 
pcaResults100 <- prcomp(transpose100)
pcaResults1k <- prcomp(transpose1k)
pcaResults5k <- prcomp(transpose5k)

#Add Metadata
condition <- c("control", "control", "control", "control", "M1-5", "M1-5","M1-5",
               "M6-10","M6-10","M6-10","M1-10","M1-10","M1-10")

transpose100 <- cbind(transpose100, condition)
transpose1k <- cbind(transpose1k, condition)
transpose5k <- cbind(transpose5k, condition)

#plot PCA results making use of ggplot2's autoplot function
autoplot(pcaResults100, data = transpose100, colour = "condition")
autoplot(pcaResults1k, data = transpose1k, colour = "condition")
autoplot(pcaResults5k, data = transpose5k, colour = "condition")


#######################################################
#Supervised Clustering: By TPMs
#With edited gene list-made more stringent
####################################################### 
library(pheatmap)
setwd("~/Desktop/Round 3.5 and 4 RNA Seq Analysis")
cts <- read.csv("R3_and_4_feature_counts.csv",header=TRUE)
rownames(cts) <- make.names(cts[,1], unique=TRUE)
length <- as.numeric(as.vector(cts[,2]))
cts <- cts[,-1]
cts <- cts[,-1]
cts_nao <- na.omit(cts)

##Function that actually works 
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

##Make the matrix 
tpms <- tpm3(cts_nao, length)
View(tpms)

genes <- c("PCNA", "MKI67", "MCM2", "TOP2A", "RIIAD1", "PAX6", "HES5", "SOX2", "THBS4",
           "PTPRZ1", "HOPX", "FAM107A", "LIFR", "SPARC", "EGFR", "SLC4A4", "SOX9", "FGFR3",
           "EDNRB", "BMPR1B", "ATP1A2", "TTYH1", "CRYAB", "SLC1A3", "SLC1A2", "SPARCL1",
           "AQP4", "ALDH1L1", "WIF1", "SPON1", "MLC1", "DIO2", "APOE", "GJA1", "CLU", "AGT", "MLC1",
           "DCX", "NEUROD2", "HMGN2", "CCND2", "NEUROD1",
           "SOX11", "TUBB3", "SNAP25", "THY1", "SYT1", "STMN2", "ENO2", "SYN1", "SLC17A6",
           "GRIN1", "GRIN2B", "RELN", "BCL11B", "NEUROD6", "CNTNAP2", "RBFOX3", "SLC17A7",
           "SLC17A6", "SLC32A1")
length(genes)
tpm_genes <- tpms[genes,]
View(tpm_genes)
#Subsetted <- tpm_genes[rowSums(tpm_genes) > 33,]

heatmap <- pheatmap(tpm_genes, cluster_rows = FALSE, scale = 'row', show_rownames = TRUE)

##########################################################
##Unsupervised Clustering: By Variance 
##########################################################
round4 <- cts[,c(21:33)]

#Variance 
V <- apply(round4, 1, var, na.rm = TRUE)
View(V)

#Sort variance in decreasing order and select top 100, 1000, and 5000 genes 
selectedGenes100 <- names(V[order(V, decreasing = TRUE)][1:100])
selectedGenes1k <- names(V[order(V, decreasing = TRUE)][1:1000])
selectedGenes5k <- names(V[order(V, decreasing = TRUE)][1:5000])

#Plot that sucker 
library(pheatmap)
genes100 <- pheatmap(round4[selectedGenes100,], scale = 'row', show_rownames = FALSE)
genes1k <- pheatmap(round4[selectedGenes1k,], scale = 'row', show_rownames = FALSE)
genes5k <- pheatmap(round4[selectedGenes5k,], scale = 'row', show_rownames = FALSE)

##########################################################
##Unsupervised Clustering: By TPMs
##########################################################

setwd("~/Desktop/Round 3.5 and 4 RNA Seq Analysis")
cts <- read.csv("R3_and_4_feature_counts.csv",header=TRUE)
rownames(cts) <- make.names(cts[,1], unique=TRUE)
length <- as.numeric(as.vector(cts[,2]))
cts <- cts[,-1]
cts <- cts[,-1]
cts_nao <- na.omit(cts)
cts_nao <- cts[,c(8:20)]

#https://www.biostars.org/p/335187/ 
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

##Make the matrix 
tpms <- tpm3(cts_nao, length)
View(tpms_filtered)

tpms_filtered <- tpms[rowSums(tpms[])>0,]
heatmap <- pheatmap(tpms_filtered, scale = 'row', show_rownames = FALSE) 
