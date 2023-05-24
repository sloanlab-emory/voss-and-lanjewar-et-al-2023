## Ligand/Receptor Project 
## Revision Experiments Bulk RNA-Seq Analysis 

## 2-15-23

## Anna Voss | Sloan Lab 

##########################################################
#Setup
########################################################## 

setwd("~/Desktop/Ligand Project/CD15 RNA Seq/")
cts <- read.csv("cd15_seq_featureconts_reordered.csv",header=TRUE)
rownames(cts) <- make.names(cts[,1], unique=TRUE)
gene_names<-cts[,1]
cts <- cts[,-1]
cts <- cts[,-1]
View(cts)

##########################################################
#Count Normalization 
########################################################## 

##CPM Normalization
library(edgeR)
cpm <- cpm(cts, log = FALSE, normalized.lib.sizes=TRUE)
#write.csv(cpm, "ligands_cd15_seq_cpm_values.csv")

##DESeq Normalization

#Set Up Col Data 
condition <- c("control", "control", "control", "5ligs", "5ligs", "5ligs",
               "BT", "BT", "BT", "DTN", "DTN", "DTN", "LIF", "LIF", "LIF",
               
               "control", "control", "control", "5ligs", "5ligs", "5ligs", "BT", 
               "BT", "BT", "DTN", "DTN", "DTN", "LIF", "LIF", "LIF",
               
               "control", "control", "control", "BT", "BT", "BT", "DTN", "DTN",
               "DTN", "LIF", "LIF", "LIF", 
               
               "control", "control", "5ligs", "BT", "BT", "BT", "DTN", "LIF",
               "LIF", "LIF", 
               
               "control", "control", "control", "5ligs", "5ligs", "5ligs", "BT",
               "BT", "BT", "DTN", "DTN", "DTN")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)
coldata

##Setup DESeq2
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#write.csv(normalized_counts, "ligands_cd15_seq_deseq_normalized_counts.csv")

##########################################################
#DESeq
########################################################## 

#############################
#5 Ligands vs. Control
#############################

ligands_v_control <- cts[,c(1:6,16:21, 31:33, 43:45, 53:58)]
View(ligands_v_control)
gene_names <- rownames(ligands_v_control)

#Set Up Col Data 
condition <- c("control", "control", "control", "5ligs", "5ligs", "5ligs",
               
               "control", "control", "control", "5ligs", "5ligs", "5ligs", 
               
               "control", "control", "control",
               
               "control", "control", "5ligs", 
               
               "control", "control", "control", "5ligs", "5ligs", "5ligs")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)
coldata

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = ligands_v_control,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)

dds$condition = relevel(dds$condition, "5ligs") 
dds$condition
resultsNames(dds)

res <- results(dds, contrast=c("condition","5ligs", "control"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "cd15 5ligs vs control.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

# Make a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = gene_names,
                x = 'log2FoldChange',
                #selectLab = c("CRYAB"),
                y = 'pvalue',
                xlim = c(-9, 9),
                ylim = c(0, 60),
                pCutoff = 0.001,
                title = "Control (left) vs. 5 ligands (right)")

#############################
#LIF/CNTF vs. Control
#############################

LIF_v_control <- cts[,c(1:3, 13:18, 28:33, 40:44, 50:55)]
View(LIF_v_control)
gene_names <- rownames(LIF_v_control)

#Set Up Col Data 
condition <- c("control", "control", "control", "LIF", "LIF", "LIF",
               
               "control", "control", "control", "LIF", "LIF", "LIF",
               
               "control", "control", "control", "LIF", "LIF", "LIF", 
               
               "control", "control", "LIF", "LIF", "LIF", 
               
               "control", "control", "control")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)
coldata

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = LIF_v_control,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)

dds$condition = relevel(dds$condition, "control") 
dds$condition
resultsNames(dds)

res <- results(dds, contrast=c("condition","LIF", "control"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "cd15 LIF_CNTF vs control.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

# Make a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = gene_names,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                ylim = c(0, 100),
                pCutoff = 0.001,
                title = "Control (left) vs. LIF/CNTF (right)")

#############################
#BMP4/TGFB2 vs. Control
#############################

BT_v_control <- cts[,c(1:3, 7:9, 16:18, 22:24, 31:36, 43:44, 46:48, 53:55, 59:61)]
View(BT_v_control)
gene_names <- rownames(BT_v_control)

#Set Up Col Data 
condition <- c("control", "control", "control", "BT", "BT", "BT",
               
               "control", "control", "control", "BT", "BT", "BT",
               
               "control", "control", "control", "BT", "BT", "BT", 
               
               "control", "control", "BT", "BT", "BT", 
               
               "control", "control", "control", "BT", "BT", "BT")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)
coldata

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = BT_v_control,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)

dds$condition = relevel(dds$condition, "BT") 
dds$condition
resultsNames(dds)

res <- results(dds, contrast=c("condition","BT", "control"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "cd15 BMP4_TGFB2 vs control.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

# Make a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = gene_names,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                ylim = c(0, 100),
                pCutoff = 0.001,
                title = "Control (left) vs. BMP4/TGFB2 (right)")

#############################
#DKK1/TSLP/NLGN1 vs. Control
#############################

DNT_v_control <- cts[,c(1:3, 10:12, 16:18, 25:27, 31:33, 37:39, 43:44, 49, 53:55, 62:64)]
View(DNT_v_control)
gene_names <- rownames(DNT_v_control)

#Set Up Col Data 
condition <- c("control", "control", "control", "DTN", "DTN", "DTN",
               
               "control", "control", "control", "DTN", "DTN", "DTN",
               
               "control", "control", "control", "DTN", "DTN", "DTN",
               
               "control", "control", "DTN",
               
               "control", "control", "control", "DTN", "DTN", "DTN")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)
coldata

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = DNT_v_control,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)

dds$condition = relevel(dds$condition, "control") 
dds$condition
resultsNames(dds)

res <- results(dds, contrast=c("condition","control", "DTN"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "cd15 DKK1_NLGN1_TSLP vs control.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

# Make a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = gene_names,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                ylim = c(0, 100),
                pCutoff = 0.001,
                title = "Control (left) vs. DKK1/NLGN1/TSLP (right)")

##########################################################
#DESeq Figures for Paper: Vocanoes unlabeled/bolded
########################################################## 

#############################
#5 Ligands vs. Control
#############################

ligands_v_control <- cts[,c(1:6,16:21, 31:33, 43:45, 53:58)]
View(ligands_v_control)
gene_names <- rownames(ligands_v_control)

#Set Up Col Data 
condition <- c("control", "control", "control", "5ligs", "5ligs", "5ligs",
               
               "control", "control", "control", "5ligs", "5ligs", "5ligs", 
               
               "control", "control", "control",
               
               "control", "control", "5ligs", 
               
               "control", "control", "control", "5ligs", "5ligs", "5ligs")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)
coldata

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = ligands_v_control,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)

dds$condition = relevel(dds$condition, "5ligs") 
dds$condition
resultsNames(dds)

res <- results(dds, contrast=c("condition","5ligs", "control"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "cd15 5ligs vs control.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

#Highlight only neuronal and astro genes of interest in volcanos
supptable5_neuronal_genes <- c("BCL11B", "CALB1", "CBLN2", "CDH8", "DCX", "EOMES",
                               "GJB7", "GRIA4", "GRIK1", "GRIN2A", "GRM1", "KCNJ3",
                               "LHX6", "NDNF", "NES", "NEUROD4", "NEUROG2", "PLD5",
                               "PVALB", "SP8", "STMN2", "SYT10", "TBR1", "TRH", "TUBB3")

supptable5_astro_genes <- c("CRYAB", "IL33", "GFAP", "ITGA7", "TNC", "ALDH1L1", 
                            "SLC1A3", "VIM", "NCAN", "CDON", "ID4", "SOX9", "IGFBPL1",
                            "PDLIM3", "ASPM", "SLC1A2", "SFRP2", "AGT", "OTX1", 
                            "GJA1", "FABP7", "PAX6", "AQP4", "FGFR3", "DIO2")
library(dplyr)

keyvals <- ifelse(
  rownames(res) %in% supptable5_neuronal_genes, 'blue',
  ifelse(rownames(res) %in% supptable5_astro_genes, 'red', 'grey'))

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'blue'] <- 'neurons'
names(keyvals)[keyvals == 'red'] <- 'astrocytes'
names(keyvals)[keyvals == 'grey'] <- 'other genes'

# Make a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                colCustom = keyvals,
                xlim = c(-9, 9),
                ylim = c(0, 60),
                pCutoff = 0.001,
                title = "Control (left) vs. 5 ligands (right)")

#############################
#5 Ligands vs. LIF/CNTF 
#############################

allligs_vs_lif <- cts[,c(4:6, 13:15, 19:21, 28:30, 40:42, 45, 50:52, 56:58)]
View(allligs_vs_lif)
gene_names <- rownames(allligs_vs_lif)

#Set Up Col Data 
condition <- c("5ligs", "5ligs", "5ligs", "LIF", "LIF", "LIF",
               
               "5ligs", "5ligs", "5ligs", "LIF", "LIF", "LIF",
               
               "LIF", "LIF", "LIF",
               
               "5ligs", "LIF", "LIF", "LIF", 
               
               "5ligs", "5ligs", "5ligs")

coldata<-cbind(condition)
coldata <- as.data.frame(coldata)
coldata$condition <- factor(coldata$condition)
coldata

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = allligs_vs_lif,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds, betaPrior=FALSE)

dds$condition = relevel(dds$condition, "5ligs") 
dds$condition
resultsNames(dds)

res <- results(dds, contrast=c("condition","5ligs", "LIF"))
ix = which.min(res$padj) # most significant
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
rownames(res)<-gene_names
write.csv(res, "cd15 all5_vs_LIF.csv")
# Shows a barplot of the top DE gene to confirm pattern of expression is what you think

# Make a volcano plot 

supptable5_neuronal_genes <- c("BCL11B", "CALB1", "CBLN2", "CDH8", "DCX", "EOMES",
                               "GJB7", "GRIA4", "GRIK1", "GRIN2A", "GRM1", "KCNJ3",
                               "LHX6", "NDNF", "NES", "NEUROD4", "NEUROG2", "PLD5",
                               "PVALB", "SP8", "STMN2", "SYT10", "TBR1", "TRH", "TUBB3")

supptable5_astro_genes <- c("CRYAB", "IL33", "GFAP", "ITGA7", "TNC", "ALDH1L1", 
                            "SLC1A3", "VIM", "NCAN", "CDON", "ID4", "SOX9", "IGFBPL1",
                            "PDLIM3", "ASPM", "SLC1A2", "SFRP2", "AGT", "OTX1", 
                            "GJA1", "FABP7", "PAX6", "AQP4", "FGFR3", "DIO2")
library(dplyr)

keyvals <- ifelse(
  rownames(res) %in% supptable5_neuronal_genes, 'blue',
  ifelse(rownames(res) %in% supptable5_astro_genes, 'red', 'grey'))

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'blue'] <- 'neurons'
names(keyvals)[keyvals == 'red'] <- 'astrocytes'
names(keyvals)[keyvals == 'grey'] <- 'other genes'

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                colCustom = keyvals,
                xlim = c(-9, 9),
                ylim = c(0, 60),
                pCutoff = 0.001,
                title = "LIF/CNTF (left) vs. All 5 ligands (right)")






