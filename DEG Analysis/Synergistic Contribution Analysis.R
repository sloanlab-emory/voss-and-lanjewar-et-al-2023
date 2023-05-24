##Synergistic Contribution Analysis 

##Anna Voss | Sloan Lab 

##1/27/22 

##################################################################################
#0.001 P-value
##################################################################################

##########################################################
#day 60-90 1363.1 Organoids (R3.5)
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis/DESeq CSVs + Volcanos")
r3_5 <- read.csv("R3_5_M6-10.csv")

#Get Up Regulated and Significant Genes: 147 
R3_5_up_DEGs <- r3_5[r3_5$pvalue <= 0.001 & r3_5$log2FoldChange >=1,] 
R3_5_up_nao <- na.omit(R3_5_up_DEGs) 
View(R3_5_up_nao)
colnames(R3_5_up_nao)[1] <- "Var1"

###############
#TGFB2 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
tgfb2 <- read.csv("R3_TGFB2.csv")

#Get Up Regulated and Significant Genes: 61
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
NLGN1 <- read.csv("R3_NLGN1.csv")

#Get Up Regulated and Significant Genes: 34
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.001 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
TSLP <- read.csv("R3_TSLP.csv")

#Get Up Regulated and Significant Genes: 65
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.001 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
DKK1 <- read.csv("R3_DKK1.csv")

#Get Up Regulated and Significant Genes: 79
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.001 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
BMP4 <- read.csv("R3_BMP4.csv")

#Get Up Regulated and Significant Genes: 61
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.001 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #300 DEGs Total (see validation below) 

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #129
two_occurances <- frequency_table[frequency_table$Freq == 2,] #25
three_occurances <- frequency_table[frequency_table$Freq == 3,] #19
four_occurances <- frequency_table[frequency_table$Freq == 4,] #16
five_occurances <- frequency_table[frequency_table$Freq == 5,] #0

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

###############
#Validation
###############

#Values in each occurance df add up to 300!! 
library(tidyverse)

#Using add_row() function to add observation to data frame  
one_occurance_1 <- one_occurance %>% add_row(Var1 = "BAMBI", Freq = 1)
one_occurance_2 <- one_occurance_1 %>% add_row(Var1 = "ID3", Freq = 1)
one_occurance_3 <- one_occurance_2 %>% add_row(Var1 = "WNT1", Freq = 1)
#These three genes appear in the line of code below: validates intersect script 

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 3 (CRYAB, FBLN2, and GPC3)
#p-val = 0.01 -> 6 (CRYAB, GJD2, MT1E, FBLN2, PENK, GPC3)
#p-val = 0.05 -> 8 (AGT, CRYAB, GJD2, MT1E, FBLN2, ELN, PENK, GPC3)
common_DEGs_1 <- inner_join(R3_5_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R3 1 I occurance and M6-10 001 sig.csv")

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 2 (RXRG and COL9A3)
#p-val = 0.01 -> 0
#p-val = 0.05 -> 1 (COL2A1)
common_DEGs_2 <- inner_join(R3_5_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R3 2 I occurances and M6-10 001 sig.csv")

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 0 
#p-val = 0.01 -> 1 (COL9A3)
#p-val = 0.05 -> 1 (COL9A3)
common_DEGs_3 <- inner_join(R3_5_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R3 3 I occurances and M6-10 001 sig.csv")

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 0
#p-val = 0.01 -> 1 (RXRG) 
#p-val = 0.05 -> 1 (RXRG) 
common_DEGs_4 <- inner_join(R3_5_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R3 4 I occurances and M6-10 001 sig.csv") 

##########################################################
#day 90-120 1363.1 organoids (R4)
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis/DESeq CSVs + Volcanos")
r4 <- read.csv("R4_M6-10.csv")

#Get Up Regulated and Significant Genes: 276
R4_up_DEGs <- r4[r4$pvalue <= 0.001 & r4$log2FoldChange >=1,] 
R4_up_nao <- na.omit(R4_up_DEGs) 
View(R4_up_nao)
colnames(R4_up_nao)[1] <- "Var1"

###############
#TGFB2 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
tgfb2 <- read.csv("R4_TGFB2.csv")

#Get Up Regulated and Significant Genes: 33
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
NLGN1 <- read.csv("R4_NLGN1.csv")

#Get Up Regulated and Significant Genes: 9
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.001 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
TSLP <- read.csv("R4_TSLP.csv")

#Get Up Regulated and Significant Genes: 28
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.001 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
DKK1 <- read.csv("R4_DKK1.csv")

#Get Up Regulated and Significant Genes: 79
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.001 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
BMP4 <- read.csv("R4_BMP4.csv")

#Get Up Regulated and Significant Genes: 58
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.001 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #207 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #164
two_occurances <- frequency_table[frequency_table$Freq == 2,] #15
three_occurances <- frequency_table[frequency_table$Freq == 3,] #3
four_occurances <- frequency_table[frequency_table$Freq == 4,] #1
five_occurances <- frequency_table[frequency_table$Freq == 5,] #0

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 23
common_DEGs_1 <- inner_join(R4_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R4 1 I occurance and M6-10 001 sig.csv")

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 5
#p-val = 0.01 -> 
#p-val = 0.05 -> 
common_DEGs_2 <- inner_join(R4_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R4 2 I occurances and M6-10 001 sig.csv")

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 0
#p-val = 0.01 -> 
#p-val = 0.05 -> 
common_DEGs_3 <- inner_join(R4_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R4 3 I occurances and M6-10 001 sig.csv")

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 1
#p-val = 0.01 -> 
#p-val = 0.05 -> 
common_DEGs_4 <- inner_join(R4_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R4 4 I occurances and M6-10 001 sig.csv") 

dim(common_DEGs_1)
dim(common_DEGs_2)
dim(common_DEGs_3)
dim(common_DEGs_4)

##########################################################
#day 60-90 8858.3 organoids (R5)
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
r5 <- read.csv("R5_M6-10.csv")

#Get Up Regulated and Significant Genes: 377
R5_up_DEGs <- r5[r5$pvalue <= 0.001 & r5$log2FoldChange >=1,] 
R5_up_nao <- na.omit(R5_up_DEGs) 
View(R5_up_nao)
colnames(R5_up_nao)[1] <- "Var1"

###############
#TGFB2 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
tgfb2 <- read.csv("R5_TGFB2.csv")

#Get Up Regulated and Significant Genes: 41
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
NLGN1 <- read.csv("R5_NLGN1.csv")

#Get Up Regulated and Significant Genes: 18
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.001 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
TSLP <- read.csv("R5_TSLP.csv")

#Get Up Regulated and Significant Genes: 20
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.001 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
DKK1 <- read.csv("R5_DKK1.csv")

#Get Up Regulated and Significant Genes: 18
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.001 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
BMP4 <- read.csv("R5_BMP4.csv")

#Get Up Regulated and Significant Genes: 421
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.001 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #518 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #446
two_occurances <- frequency_table[frequency_table$Freq == 2,] #10
three_occurances <- frequency_table[frequency_table$Freq == 3,] #1
four_occurances <- frequency_table[frequency_table$Freq == 4,] #6
five_occurances <- frequency_table[frequency_table$Freq == 5,] #5

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 321
common_DEGs_1 <- inner_join(R5_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R5 1 I occurance and M6-10 001 sig.csv")
dim(common_DEGs_1)

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 3
common_DEGs_2 <- inner_join(R5_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R5 2 I occurances and M6-10 001 sig.csv")
dim(common_DEGs_2)

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 0
common_DEGs_3 <- inner_join(R5_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R5 3 I occurances and M6-10 001 sig.csv")
dim(common_DEGs_3)

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 4
common_DEGs_4 <- inner_join(R5_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R5 4 I occurances and M6-10 001 sig.csv") 
dim(common_DEGs_4)

#Common DEGs: Five occurances + 6-10: 
#p-val = 0.001 -> 5
common_DEGs_5 <- inner_join(R5_up_nao, five_occurances, by = "Var1")
write.csv(common_DEGs_5, "Common DEGs R5 5 I occurances and M6-10 001 sig.csv") 
dim(common_DEGs_5)

##########################################################
#day 90-120 8858.3 organoids (R6)
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
r6 <- read.csv("R6_M6-10.csv")

#Get Up Regulated and Significant Genes: 200
R6_up_DEGs <- r6[r6$pvalue <= 0.001 & r6$log2FoldChange >=1,] 
R6_up_nao <- na.omit(R6_up_DEGs) 
View(R6_up_nao)
colnames(R6_up_nao)[1] <- "Var1"

###############
#TGFB2 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
tgfb2 <- read.csv("R6_TGFB2.csv")

#Get Up Regulated and Significant Genes: 14
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
NLGN1 <- read.csv("R6_NLGN1.csv")

#Get Up Regulated and Significant Genes: 44
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.001 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
TSLP <- read.csv("R6_TSLP.csv")

#Get Up Regulated and Significant Genes: 56
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.001 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
DKK1 <- read.csv("R6_DKK1.csv")

#Get Up Regulated and Significant Genes: 32
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.001 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
BMP4 <- read.csv("R6_BMP4.csv")

#Get Up Regulated and Significant Genes: 180
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.001 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #326 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #219
two_occurances <- frequency_table[frequency_table$Freq == 2,] #22
three_occurances <- frequency_table[frequency_table$Freq == 3,] #4
four_occurances <- frequency_table[frequency_table$Freq == 4,] #4
five_occurances <- frequency_table[frequency_table$Freq == 5,] #7

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 130
common_DEGs_1 <- inner_join(R6_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R6 1 I occurance and M6-10 001 sig.csv")
dim(common_DEGs_1)

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 3
common_DEGs_2 <- inner_join(R6_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R6 2 I occurances and M6-10 001 sig.csv")
dim(common_DEGs_2)

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 0
common_DEGs_3 <- inner_join(R6_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R6 3 I occurances and M6-10 001 sig.csv")
dim(common_DEGs_3)

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 3
common_DEGs_4 <- inner_join(R6_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R6 4 I occurances and M6-10 001 sig.csv") 
dim(common_DEGs_4)

#Common DEGs: Five occurances + 6-10: 
#p-val = 0.001 -> 5
common_DEGs_5 <- inner_join(R6_up_nao, five_occurances, by = "Var1")
write.csv(common_DEGs_5, "Common DEGs R6 5 I occurances and M6-10 001 sig.csv") 
dim(common_DEGs_5)

####################################################################################################
#Venn Diagram Comparisons: Overlap of TGFB2 with M6-10, R3 (Up Regualted Only)
####################################################################################################
#TGFB2 I
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
library(dplyr)

#R3 TGFB2: 61 DEGs 
R3_TGFB2 <- read.csv("R3_TGFB2.csv")
R3_TGFB2_DEGs <- R3_TGFB2[R3_TGFB2$pvalue <= 0.001 & R3_TGFB2$log2FoldChange >=1,]
R3_TGFB2_nao <- na.omit(R3_TGFB2_DEGs)

##########################
#M6-10
##########################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/R3 Mega Wells DESeq")

R3_M6_10 <- read.csv("R3_M6-10.csv")
R3_M6_10_DEGs <- R3_M6_10[R3_M6_10$pvalue <= 0.001 & R3_M6_10$log2FoldChange >=1,]
R3_M6_10_nao <- na.omit(R3_M6_10_DEGs)

#Common DEGs R3 TGFB2 and M6-10: 12 for R3 
common_DEGs <- inner_join(R3_TGFB2_nao, R3_M6_10_nao, by = "X")

#Determining only TGFB2 DEGs 
#Same TGFB2 and other ligands overlap so only filter away other ligands
TGFB2_only <- anti_join(R3_TGFB2_nao, R3_M6_10_nao, by = "X")

#Determining only M6-10 DEGs 
R3_M6_10_only <- anti_join(R3_M6_10_nao, R3_TGFB2_nao, by = "X")

#Write CSV's
setwd("~/Desktop/Venn Diagram Out2")
write.csv(common_DEGs, "R3 Common DEGs TGFB2 + M6-10.csv")
write.csv(TGFB2_only, "R3 M6-10 TGFB2_only_DEGs.csv")
write.csv(R3_M6_10_only, "R3 M6-10 Only DEGs.csv")

#R3.5
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis/DESeq CSVs + Volcanos")
R3_5_M6_10 <- read.csv("R3_5_M6-10.csv")
R3_5_M6_10_DEGs <- R3_5_M6_10[R3_5_M6_10$pvalue <= 0.001 & R3_5_M6_10$log2FoldChange >=1,]
R3_5_M6_10_nao <- na.omit(R3_5_M6_10_DEGs)

common_DEGs <- inner_join(R3_TGFB2_nao, R3_5_M6_10_nao, by = "X")
View(common_DEGs)

TGFB2_only <- anti_join(R3_TGFB2_nao, R3_5_M6_10_nao, by = "X")

#Determining only M6-10 DEGs 
R3_M6_10_only <- anti_join(R3_5_M6_10_nao, R3_TGFB2_nao, by = "X")

##################################################################################
#0.01 P-value
##################################################################################


#########################
#R3.5
#########################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis/DESeq CSVs + Volcanos")
r3_5 <- read.csv("R3_5_M6-10.csv")

#Get Up Regulated and Significant Genes: 192 
R3_5_up_DEGs <- r3_5[r3_5$pvalue <= 0.01 & r3_5$log2FoldChange >=1,] 
R3_5_up_nao <- na.omit(R3_5_up_DEGs) 
View(R3_5_up_nao)
colnames(R3_5_up_nao)[1] <- "Var1"

###############
#TGFB2 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
tgfb2 <- read.csv("R3_TGFB2.csv")

#Get Up Regulated and Significant Genes: 73
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.01 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
NLGN1 <- read.csv("R3_NLGN1.csv")

#Get Up Regulated and Significant Genes: 53
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.01 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
TSLP <- read.csv("R3_TSLP.csv")

#Get Up Regulated and Significant Genes: 106
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.01 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
DKK1 <- read.csv("R3_DKK1.csv")

#Get Up Regulated and Significant Genes: 118
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.01 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
BMP4 <- read.csv("R3_BMP4.csv")

#Get Up Regulated and Significant Genes: 61
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.01 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

##########################################################
#Frequency of I ligand DEGs 
##########################################################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #467 DEGs Total    

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #230
two_occurances <- frequency_table[frequency_table$Freq == 2,] #39
three_occurances <- frequency_table[frequency_table$Freq == 3,] #25
four_occurances <- frequency_table[frequency_table$Freq == 4,] #21
five_occurances <- frequency_table[frequency_table$Freq == 5,] #0

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

##########################################################
#Intersect of I frequencies with 6-10 DEGs  
##########################################################
library(dplyr)
setwd("~/Desktop/")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 3 (CRYAB, FBLN2, and GPC3)
#p-val = 0.01 -> 6 (CRYAB, GJD2, MT1E, FBLN2, PENK, GPC3)
#p-val = 0.05 -> 8 (AGT, CRYAB, GJD2, MT1E, FBLN2, ELN, PENK, GPC3)
common_DEGs_1 <- inner_join(R3_5_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R3 1 I occurance and M6-10 01 sig.csv")

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 2 (RXRG and COL9A3)
#p-val = 0.01 -> 0
#p-val = 0.05 -> 1 (COL2A1)
common_DEGs_2 <- inner_join(R3_5_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R3 2 I occurances and M6-10 01 sig.csv")

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 0 
#p-val = 0.01 -> 1 (COL9A3)
#p-val = 0.05 -> 1 (COL9A3)
common_DEGs_3 <- inner_join(R3_5_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R3 3 I occurances and M6-10 01 sig.csv")

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 0
#p-val = 0.01 -> 1 (RXRG) 
#p-val = 0.05 -> 1 (RXRG) 
common_DEGs_4 <- inner_join(R3_5_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R3 4 I occurances and M6-10 01 sig.csv") 

dim(common_DEGs_1)
dim(common_DEGs_2)
dim(common_DEGs_3)
dim(common_DEGs_4)

#########################
#R4
#########################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis/DESeq CSVs + Volcanos")
r4 <- read.csv("R4_M6-10.csv")

#Get Up Regulated and Significant Genes: 326
R4_up_DEGs <- r4[r4$pvalue <= 0.01 & r4$log2FoldChange >=1,] 
R4_up_nao <- na.omit(R4_up_DEGs) 
View(R4_up_nao)
colnames(R4_up_nao)[1] <- "Var1"

###############
#TGFB2 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
tgfb2 <- read.csv("R4_TGFB2.csv")

#Get Up Regulated and Significant Genes: 63
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.01 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
NLGN1 <- read.csv("R4_NLGN1.csv")

#Get Up Regulated and Significant Genes: 19
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.01 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
TSLP <- read.csv("R4_TSLP.csv")

#Get Up Regulated and Significant Genes: 47
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.01 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
DKK1 <- read.csv("R4_DKK1.csv")

#Get Up Regulated and Significant Genes: 159
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.01 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
BMP4 <- read.csv("R4_BMP4.csv")

#Get Up Regulated and Significant Genes: 107
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.01 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

##########################################################
#Frequency of I ligand DEGs 
##########################################################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #395 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #296
two_occurances <- frequency_table[frequency_table$Freq == 2,] #37
three_occurances <- frequency_table[frequency_table$Freq == 3,] #7
four_occurances <- frequency_table[frequency_table$Freq == 4,] #1
five_occurances <- frequency_table[frequency_table$Freq == 5,] #0

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

##########################################################
#Intersect of I frequencies with 6-10 DEGs  
##########################################################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.01 -> 37
common_DEGs_1 <- inner_join(R4_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R4 1 I occurance and M6-10 01 sig.csv")

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.01 -> 8
common_DEGs_2 <- inner_join(R4_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R4 2 I occurances and M6-10 01 sig.csv")

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.01 -> 3
common_DEGs_3 <- inner_join(R4_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R4 3 I occurances and M6-10 01 sig.csv")

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.01 -> 1
common_DEGs_4 <- inner_join(R4_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R4 4 I occurances and M6-10 01 sig.csv") 

dim(common_DEGs_1)
dim(common_DEGs_2)
dim(common_DEGs_3)
dim(common_DEGs_4)

#########################
#R5
#########################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
r5 <- read.csv("R5_M6-10.csv")

#Get Up Regulated and Significant Genes: 432
R5_up_DEGs <- r5[r5$pvalue <= 0.01 & r5$log2FoldChange >=1,] 
R5_up_nao <- na.omit(R5_up_DEGs) 
View(R5_up_nao)
colnames(R5_up_nao)[1] <- "Var1"

###############
#TGFB2 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
tgfb2 <- read.csv("R5_TGFB2.csv")

#Get Up Regulated and Significant Genes: 73
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.01 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
NLGN1 <- read.csv("R5_NLGN1.csv")

#Get Up Regulated and Significant Genes: 26
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.01 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
TSLP <- read.csv("R5_TSLP.csv")

#Get Up Regulated and Significant Genes: 44
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.01 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
DKK1 <- read.csv("R5_DKK1.csv")

#Get Up Regulated and Significant Genes: 63
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.01 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
BMP4 <- read.csv("R5_BMP4.csv")

#Get Up Regulated and Significant Genes: 553
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.01 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

##########################################################
#Frequency of I ligand DEGs 
##########################################################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #642 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #530
two_occurances <- frequency_table[frequency_table$Freq == 2,] #25
three_occurances <- frequency_table[frequency_table$Freq == 3,] #3
four_occurances <- frequency_table[frequency_table$Freq == 4,] #7
five_occurances <- frequency_table[frequency_table$Freq == 5,] #5

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

##########################################################
#Intersect of I frequencies with 6-10 DEGs  
##########################################################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.01 -> 351
common_DEGs_1 <- inner_join(R5_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R5 1 I occurance and M6-10 01 sig.csv")
dim(common_DEGs_1)

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.01 -> 13
common_DEGs_2 <- inner_join(R5_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R5 2 I occurances and M6-10 01 sig.csv")
dim(common_DEGs_2)

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.01 -> 1
common_DEGs_3 <- inner_join(R5_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R5 3 I occurances and M6-10 01 sig.csv")
dim(common_DEGs_3)

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.01 -> 5
common_DEGs_4 <- inner_join(R5_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R5 4 I occurances and M6-10 01 sig.csv") 
dim(common_DEGs_4)

#Common DEGs: Five occurances + 6-10: 
#p-val = 0.01 -> 
common_DEGs_5 <- inner_join(R5_up_nao, five_occurances, by = "Var1")
write.csv(common_DEGs_5, "Common DEGs R5 5 I occurances and M6-10 01 sig.csv") 
dim(common_DEGs_5)

#########################
#R6
#########################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
r6 <- read.csv("R6_M6-10.csv")

#Get Up Regulated and Significant Genes: 278
R6_up_DEGs <- r6[r6$pvalue <= 0.01 & r6$log2FoldChange >=1,] 
R6_up_nao <- na.omit(R6_up_DEGs) 
View(R6_up_nao)
colnames(R6_up_nao)[1] <- "Var1"

###############
#TGFB2 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
tgfb2 <- read.csv("R6_TGFB2.csv")

#Get Up Regulated and Significant Genes: 23
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.01 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
NLGN1 <- read.csv("R6_NLGN1.csv")

#Get Up Regulated and Significant Genes: 70
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.01 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
TSLP <- read.csv("R6_TSLP.csv")

#Get Up Regulated and Significant Genes: 103
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.01 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
DKK1 <- read.csv("R6_DKK1.csv")

#Get Up Regulated and Significant Genes: 49
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.01 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
BMP4 <- read.csv("R6_BMP4.csv")

#Get Up Regulated and Significant Genes: 212
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.01 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

##########################################################
#Frequency of I ligand DEGs 
##########################################################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #457 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #261
two_occurances <- frequency_table[frequency_table$Freq == 2,] #35
three_occurances <- frequency_table[frequency_table$Freq == 3,] #22
four_occurances <- frequency_table[frequency_table$Freq == 4,] #5
five_occurances <- frequency_table[frequency_table$Freq == 5,] #8

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

##########################################################
#Intersect of I frequencies with 6-10 DEGs  
##########################################################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.01 -> 145
common_DEGs_1 <- inner_join(R6_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R6 1 I occurance and M6-10 01 sig.csv")
dim(common_DEGs_1)

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.01 -> 4
common_DEGs_2 <- inner_join(R6_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R6 2 I occurances and M6-10 01 sig.csv")
dim(common_DEGs_2)

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.01 -> 2
common_DEGs_3 <- inner_join(R6_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R6 3 I occurances and M6-10 01 sig.csv")
dim(common_DEGs_3)

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.01 -> 1
common_DEGs_4 <- inner_join(R6_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R6 4 I occurances and M6-10 01 sig.csv") 
dim(common_DEGs_4)

#Common DEGs: Five occurances + 6-10: 
#p-val = 0.01 -> 6
common_DEGs_5 <- inner_join(R6_up_nao, five_occurances, by = "Var1")
write.csv(common_DEGs_5, "Common DEGs R6 5 I occurances and M6-10 01 sig.csv") 
dim(common_DEGs_5)


##################################################################################
#0.05 P-value
##################################################################################

##########################################################
#R3.5
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis/DESeq CSVs + Volcanos")
r3_5 <- read.csv("R3_5_M6-10.csv")

#Get Up Regulated and Significant Genes: 238 
R3_5_up_DEGs <- r3_5[r3_5$pvalue <= 0.05 & r3_5$log2FoldChange >=1,] 
R3_5_up_nao <- na.omit(R3_5_up_DEGs) 
View(R3_5_up_nao)
colnames(R3_5_up_nao)[1] <- "Var1"

###############
#TGFB2 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
tgfb2 <- read.csv("R3_TGFB2.csv")

#Get Up Regulated and Significant Genes: 81
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
NLGN1 <- read.csv("R3_NLGN1.csv")

#Get Up Regulated and Significant Genes: 76
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.05 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
TSLP <- read.csv("R3_TSLP.csv")

#Get Up Regulated and Significant Genes: 145
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.05 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
DKK1 <- read.csv("R3_DKK1.csv")

#Get Up Regulated and Significant Genes: 156
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.05 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R3
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
BMP4 <- read.csv("R3_BMP4.csv")

#Get Up Regulated and Significant Genes: 195
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.05 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #653 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #346
two_occurances <- frequency_table[frequency_table$Freq == 2,] #64
three_occurances <- frequency_table[frequency_table$Freq == 3,] #29
four_occurances <- frequency_table[frequency_table$Freq == 4,] #23
five_occurances <- frequency_table[frequency_table$Freq == 5,] #0

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 3 (CRYAB, FBLN2, and GPC3)
#p-val = 0.01 -> 6 (CRYAB, GJD2, MT1E, FBLN2, PENK, GPC3)
#p-val = 0.05 -> 8 (AGT, CRYAB, GJD2, MT1E, FBLN2, ELN, PENK, GPC3)
common_DEGs_1 <- inner_join(R3_5_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R3 1 I occurance and M6-10 05 sig.csv")

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 2 (RXRG and COL9A3)
#p-val = 0.01 -> 0
#p-val = 0.05 -> 1 (COL2A1)
common_DEGs_2 <- inner_join(R3_5_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R3 2 I occurances and M6-10 05 sig.csv")

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 0 
#p-val = 0.01 -> 1 (COL9A3)
#p-val = 0.05 -> 1 (COL9A3)
common_DEGs_3 <- inner_join(R3_5_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R3 3 I occurances and M6-10 05 sig.csv")

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 0
#p-val = 0.01 -> 1 (RXRG) 
#p-val = 0.05 -> 1 (RXRG) 
common_DEGs_4 <- inner_join(R3_5_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R3 4 I occurances and M6-10 05 sig.csv") 

dim(common_DEGs_1)
dim(common_DEGs_2)
dim(common_DEGs_3)
dim(common_DEGs_4)

##########################################################
#R4
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Round 3.5 and 4 RNA Seq Analysis/DESeq CSVs + Volcanos")
r4 <- read.csv("R4_M6-10.csv")

#Get Up Regulated and Significant Genes: 375
R4_up_DEGs <- r4[r4$pvalue <= 0.05 & r4$log2FoldChange >=1,] 
R4_up_nao <- na.omit(R4_up_DEGs) 
View(R4_up_nao)
colnames(R4_up_nao)[1] <- "Var1"

###############
#TGFB2 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
tgfb2 <- read.csv("R4_TGFB2.csv")

#Get Up Regulated and Significant Genes: 77
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
NLGN1 <- read.csv("R4_NLGN1.csv")

#Get Up Regulated and Significant Genes: 26
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.05 & NLGN1$log2FoldChange >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
TSLP <- read.csv("R4_TSLP.csv")

#Get Up Regulated and Significant Genes: 63
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.05 & TSLP$log2FoldChange >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
DKK1 <- read.csv("R4_DKK1.csv")

#Get Up Regulated and Significant Genes: 265
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.05 & DKK1$log2FoldChange >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
BMP4 <- read.csv("R4_BMP4.csv")

#Get Up Regulated and Significant Genes: 172
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.05 & BMP4$log2FoldChange >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #603 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #436
two_occurances <- frequency_table[frequency_table$Freq == 2,] #59
three_occurances <- frequency_table[frequency_table$Freq == 3,] #11
four_occurances <- frequency_table[frequency_table$Freq == 4,] #4
five_occurances <- frequency_table[frequency_table$Freq == 5,] #0

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 23
common_DEGs_1 <- inner_join(R4_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R4 1 I occurance and M6-10 05 sig.csv")

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 5
#p-val = 0.01 -> 
#p-val = 0.05 -> 
common_DEGs_2 <- inner_join(R4_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R4 2 I occurances and M6-10 05 sig.csv")

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 0
#p-val = 0.01 -> 
#p-val = 0.05 -> 
common_DEGs_3 <- inner_join(R4_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R4 3 I occurances and M6-10 05 sig.csv")

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 1
#p-val = 0.01 -> 
#p-val = 0.05 -> 
common_DEGs_4 <- inner_join(R4_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R4 4 I occurances and M6-10 05 sig.csv") 

dim(common_DEGs_1)
dim(common_DEGs_2)
dim(common_DEGs_3)
dim(common_DEGs_4)

##########################################################
#R5
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
r5 <- read.csv("R5_M6-10.csv")

#Get Up Regulated and Significant Genes: 485
R5_up_DEGs <- r5[r5$pvalue <= 0.05 & abs(r5$log2FoldChange) >=1,] 
R5_up_nao <- na.omit(R5_up_DEGs) 
View(R5_up_nao)
colnames(R5_up_nao)[1] <- "Var1"

###############
#TGFB2 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
tgfb2 <- read.csv("R5_TGFB2.csv")

#Get Up Regulated and Significant Genes: 94
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
NLGN1 <- read.csv("R5_NLGN1.csv")

#Get Up Regulated and Significant Genes: 42
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.05 & abs(NLGN1$log2FoldChange) >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
TSLP <- read.csv("R5_TSLP.csv")

#Get Up Regulated and Significant Genes: 44
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.05 & abs(TSLP$log2FoldChange) >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
DKK1 <- read.csv("R5_DKK1.csv")

#Get Up Regulated and Significant Genes: 63
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.05 & abs(DKK1$log2FoldChange) >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
BMP4 <- read.csv("R5_BMP4.csv")

#Get Up Regulated and Significant Genes: 421
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.05 & abs(BMP4$log2FoldChange) >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #796 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #446
two_occurances <- frequency_table[frequency_table$Freq == 2,] #10
three_occurances <- frequency_table[frequency_table$Freq == 3,] #1
four_occurances <- frequency_table[frequency_table$Freq == 4,] #6
five_occurances <- frequency_table[frequency_table$Freq == 5,] #5

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 389
common_DEGs_1 <- inner_join(R5_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R5 1 I occurance and M6-10 05 sig up and down.csv")
dim(common_DEGs_1)

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 19
common_DEGs_2 <- inner_join(R5_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R5 2 I occurances and M6-10 05 sig up and down.csv")
dim(common_DEGs_2)

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 7
common_DEGs_3 <- inner_join(R5_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R5 3 I occurances and M6-10 05 sig up and down.csv")
dim(common_DEGs_3)

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 5
common_DEGs_4 <- inner_join(R5_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R5 4 I occurances and M6-10 05 up and down sig.csv") 
dim(common_DEGs_4)

#Common DEGs: Five occurances + 6-10: 
#p-val = 0.001 -> 5
common_DEGs_5 <- inner_join(R5_up_nao, five_occurances, by = "Var1")
write.csv(common_DEGs_5, "Common DEGs R5 5 I occurances and M6-10 05 up and down sig.csv") 
dim(common_DEGs_5)

##########################################################
#R6
##########################################################
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
r6 <- read.csv("R6_M6-10.csv")

#Get Up Regulated and Significant Genes: 278
R6_up_DEGs <- r6[r6$pvalue <= 0.05 & abs(r6$log2FoldChange) >=1,] 
R6_up_nao <- na.omit(R6_up_DEGs) 
View(R6_up_nao)
colnames(R6_up_nao)[1] <- "Var1"

###############
#TGFB2 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
tgfb2 <- read.csv("R6_TGFB2.csv")

#Get Up Regulated and Significant Genes: 42
TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
View(TGFB2_up_nao)

###############
#NLGN1 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
NLGN1 <- read.csv("R6_NLGN1.csv")

#Get Up Regulated and Significant Genes: 96
NLGN1_up_DEGs <- NLGN1[NLGN1$pvalue <= 0.05 & abs(NLGN1$log2FoldChange) >=1,] 
NLGN1_up_nao <- na.omit(NLGN1_up_DEGs) 
View(NLGN1_up_nao)

###############
#TSLP R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
TSLP <- read.csv("R6_TSLP.csv")

#Get Up Regulated and Significant Genes: 165
TSLP_up_DEGs <- TSLP[TSLP$pvalue <= 0.05 & abs(TSLP$log2FoldChange) >=1,] 
TSLP_up_nao <- na.omit(TSLP_up_DEGs) 
View(TSLP_up_nao)

###############
#DKK1 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
DKK1 <- read.csv("R6_DKK1.csv")

#Get Up Regulated and Significant Genes: 65 
DKK1_up_DEGs <- DKK1[DKK1$pvalue <= 0.05 & abs(DKK1$log2FoldChange) >=1,] 
DKK1_up_nao <- na.omit(DKK1_up_DEGs) 
View(DKK1_up_nao)

###############
#BMP4 R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos")
BMP4 <- read.csv("R6_BMP4.csv")

#Get Up Regulated and Significant Genes: 180
BMP4_up_DEGs <- BMP4[BMP4$pvalue <= 0.05 & abs(BMP4$log2FoldChange) >=1,] 
BMP4_up_nao <- na.omit(BMP4_up_DEGs) 
View(BMP4_up_nao)

############################
#Frequency of I ligand DEGs 
############################

#Create data frame with all DEGs from all I ligands 
I_combined <- rbind(TGFB2_up_nao, NLGN1_up_nao, TSLP_up_nao, DKK1_up_nao, BMP4_up_nao)
View(I_combined) #628 DEGs Total

#Get frequency of each gene value in I_combined
frequency_table <- table(I_combined$X)
frequency_table <- as.data.frame(frequency_table)

#Identifiy genes present, 1, 2, 3, and 4 times in I_combined
one_occurance <- frequency_table[frequency_table$Freq == 1,] #315
two_occurances <- frequency_table[frequency_table$Freq == 2,] #51
three_occurances <- frequency_table[frequency_table$Freq == 3,] #31
four_occurances <- frequency_table[frequency_table$Freq == 4,] #17
five_occurances <- frequency_table[frequency_table$Freq == 5,] #10

dim(one_occurance)
dim(two_occurances)
dim(three_occurances)
dim(four_occurances)
dim(five_occurances)

############################
#Intersect of I frequencies with 6-10 DEGs  
############################
library(dplyr)
setwd("~/Desktop")

#Common DEGs: One occurance + 6-10: 
#p-val = 0.001 -> 160
common_DEGs_1 <- inner_join(R6_up_nao, one_occurance, by = "Var1")
write.csv(common_DEGs_1, "Common DEGs R6 1 I occurance and M6-10 05 sig up and down.csv")
dim(common_DEGs_1)

#Common DEGs: Two occurances + 6-10: 
#p-val = 0.001 -> 9
common_DEGs_2 <- inner_join(R6_up_nao, two_occurances, by = "Var1")
write.csv(common_DEGs_2, "Common DEGs R6 2 I occurances and M6-10 05 sig up and down.csv")
dim(common_DEGs_2)

#Common DEGs: Three occurances + 6-10: 
#p-val = 0.001 -> 6
common_DEGs_3 <- inner_join(R6_up_nao, three_occurances, by = "Var1")
write.csv(common_DEGs_3, "Common DEGs R6 3 I occurances and M6-10 05 sig up and down.csv")
dim(common_DEGs_3)

#Common DEGs: Four occurances + 6-10: 
#p-val = 0.001 -> 6
common_DEGs_4 <- inner_join(R6_up_nao, four_occurances, by = "Var1")
write.csv(common_DEGs_4, "Common DEGs R6 4 I occurances and M6-10 05 sig up and down.csv") 
dim(common_DEGs_4)

#Common DEGs: Five occurances + 6-10: 
#p-val = 0.001 -> 6
common_DEGs_5 <- inner_join(R6_up_nao, five_occurances, by = "Var1")
write.csv(common_DEGs_5, "Common DEGs R6 5 I occurances and M6-10 05 sig up and down.csv") 
dim(common_DEGs_5)

##################################################################################
#What Ligand(s) are the one-occurance overlaps from: 
##################################################################################
library(dplyr)

###############
#R3.5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R3.5 001 sig ") 
r3_5_one_occurance_001 <- read.csv("Common DEGs R3 1 I occurance and M6-10 001 sig.csv")
View(r3_5_one_occurance_001)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
tgfb2 <- read.csv("R3_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: FBLN2
tgfb2_overlap_r3_5 <- inner_join(r3_5_one_occurance_001, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r3_5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
nlgn1 <- read.csv("R3_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.001 & nlgn1$log2FoldChange >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 0
nlgn1_overlap_r3_5 <- inner_join(r3_5_one_occurance_001, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r3_5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
tslp <- read.csv("R3_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.001 & tslp$log2FoldChange >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 0
tslp_overlap_r3_5 <- inner_join(r3_5_one_occurance_001, tslp_up_nao, by = "Var1")
View(tslp_overlap_r3_5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
dkk1 <- read.csv("R3_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.001 & dkk1$log2FoldChange >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: CRYAB, GPC3
dkk1_overlap_r3_5 <- inner_join(r3_5_one_occurance_001, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r3_5)


###############
#R4
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R4 001 sig") 
r4_one_occurance_001 <- read.csv("Common DEGs R4 1 I occurance and M6-10 001 sig.csv")
View(r4_one_occurance_001)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
tgfb2 <- read.csv("R4_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 15
#CYR61, ASPM, EMP1, COL2A1, E2F7, NUSAP1, MFAP4, TOP2A, UBE2C, RASL11B, ADAMTS16, LINC01021, KIF20A, MSX2, ELN
tgfb2_overlap_r4 <- inner_join(r4_one_occurance_001, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r4)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
nlgn1 <- read.csv("R4_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.001 & nlgn1$log2FoldChange >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 0
nlgn1_overlap_r4 <- inner_join(r4_one_occurance_001, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r4)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
tslp <- read.csv("R4_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.001 & tslp$log2FoldChange >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 0
tslp_overlap_r4 <- inner_join(r4_one_occurance_001, tslp_up_nao, by = "Var1")
View(tslp_overlap_r4)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 4")
dkk1 <- read.csv("R4_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.001 & dkk1$log2FoldChange >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 8
#MFAP, MRLN, CDON, CDH15, CYBRD1, INSM1, TBX1, FBLN1 
dkk1_overlap_r4 <- inner_join(r4_one_occurance_001, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r4)

###############
#R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R5 001 sig") 
r5_one_occurance_001 <- read.csv("Common DEGs R5 1 I occurance and M6-10 001 sig.csv")
View(r5_one_occurance_001)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R5_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 1
#CDC25C
tgfb2_overlap_r5 <- inner_join(r5_one_occurance_001, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R5_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.001 & nlgn1$log2FoldChange >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 0
nlgn1_overlap_r5 <- inner_join(r5_one_occurance_001, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R5_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.001 & tslp$log2FoldChange >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 0
tslp_overlap_r5 <- inner_join(r5_one_occurance_001, tslp_up_nao, by = "Var1")
View(tslp_overlap_r5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R5_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.001 & dkk1$log2FoldChange >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 0
dkk1_overlap_r5 <- inner_join(r5_one_occurance_001, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r5)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R5_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.001 & bmp4$log2FoldChange >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 0
bmp4_overlap_r5 <- inner_join(r5_one_occurance_001, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r5)

###############
#R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R6 001 sig") 
r6_one_occurance_001 <- read.csv("Common DEGs R6 1 I occurance and M6-10 001 sig.csv")
View(r6_one_occurance_001)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R6_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.001 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 0
tgfb2_overlap_r6 <- inner_join(r6_one_occurance_001, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r6)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R6_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.001 & nlgn1$log2FoldChange >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 0
nlgn1_overlap_r6 <- inner_join(r6_one_occurance_001, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r6)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R6_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.001 & tslp$log2FoldChange >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 0
tslp_overlap_r6 <- inner_join(r6_one_occurance_001, tslp_up_nao, by = "Var1")
View(tslp_overlap_r6)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R6_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.001 & dkk1$log2FoldChange >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 1
#MIRLET7BHG
dkk1_overlap_r6 <- inner_join(r6_one_occurance_001, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r6)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R6_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.001 & bmp4$log2FoldChange >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 0
bmp4_overlap_r6 <- inner_join(r6_one_occurance_001, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r6) 

##################################################################################
#What Ligand(s) are the one-occurance overlaps from: 0.05 Sig 
##################################################################################
library(dplyr)

###############
#R3.5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R3.5 05 sig/") 
r3_5_one_occurance_05 <- read.csv("Common DEGs R3 1 I occurance and M6-10 05 sig.csv")
View(r3_5_one_occurance_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
tgfb2 <- read.csv("R3_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & tgfb2$log2FoldChange >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: FBLN2
tgfb2_overlap_r3_5 <- inner_join(r3_5_one_occurance_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r3_5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
nlgn1 <- read.csv("R3_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & nlgn1$log2FoldChange >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 2 (AGT, ELN)
nlgn1_overlap_r3_5 <- inner_join(r3_5_one_occurance_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r3_5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
tslp <- read.csv("R3_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & tslp$log2FoldChange >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 0
tslp_overlap_r3_5 <- inner_join(r3_5_one_occurance_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r3_5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 3 and 4 Individual Bulk /DESeq Round 3")
dkk1 <- read.csv("R3_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & dkk1$log2FoldChange >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: CRYAB, GPC3
dkk1_overlap_r3_5 <- inner_join(r3_5_one_occurance_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r3_5)

###############
#R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R5 05 sig/") 
r5_one_occurance_05 <- read.csv("Common DEGs R5 1 I occurance and M6-10 05 sig up and down.csv")
View(r5_one_occurance_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R5_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 0
tgfb2_overlap_r5 <- inner_join(r5_one_occurance_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R5_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 0
nlgn1_overlap_r5 <- inner_join(r5_one_occurance_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R5_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 0
tslp_overlap_r5 <- inner_join(r5_one_occurance_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R5_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 0
dkk1_overlap_r5 <- inner_join(r5_one_occurance_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r5)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R5_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 0
bmp4_overlap_r5 <- inner_join(r5_one_occurance_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r5)

###############
#R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R6 05 sig/") 
r6_one_occurance_05 <- read.csv("Common DEGs R6 1 I occurance and M6-10 05 sig up and down.csv")
View(r6_one_occurance_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R6_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 0
tgfb2_overlap_r6 <- inner_join(r6_one_occurance_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r6)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R6_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 0
nlgn1_overlap_r6 <- inner_join(r6_one_occurance_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r6)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R6_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 0
tslp_overlap_r6 <- inner_join(r6_one_occurance_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r6)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R6_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 1
#MIRLET7BHG
dkk1_overlap_r6 <- inner_join(r6_one_occurance_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r6)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R6_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 0
bmp4_overlap_r6 <- inner_join(r6_one_occurance_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r6) 

##################################################################################
#What Ligand(s) are the two-occurance overlaps from: 0.05 Sig 
##################################################################################
library(dplyr)

###############
#R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R5 05 sig/") 
r5_two_occurances_05 <- read.csv("Common DEGs R5 2 I occurances and M6-10 05 sig up and down.csv")
View(r5_two_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R5_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap:
tgfb2_overlap_r5 <- inner_join(r5_two_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R5_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r5 <- inner_join(r5_two_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R5_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r5 <- inner_join(r5_two_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R5_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 
dkk1_overlap_r5 <- inner_join(r5_two_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r5)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R5_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 
bmp4_overlap_r5 <- inner_join(r5_two_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r5)

###############
#R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R6 05 sig/") 
r6_two_occurances_05 <- read.csv("Common DEGs R6 2 I occurances and M6-10 05 sig up and down.csv")
View(r6_two_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R6_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 
tgfb2_overlap_r6 <- inner_join(r6_two_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r6)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R6_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r6 <- inner_join(r6_two_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r6)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R6_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r6 <- inner_join(r6_two_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r6)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R6_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap:
dkk1_overlap_r6 <- inner_join(r6_two_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r6)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R6_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap:
bmp4_overlap_r6 <- inner_join(r6_two_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r6) 

##################################################################################
#What Ligand(s) are the three-occurance overlaps from: 0.05 Sig 
##################################################################################
library(dplyr)

###############
#R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R5 05 sig/") 
r5_three_occurances_05 <- read.csv("Common DEGs R5 3 I occurances and M6-10 05 sig up and down.csv")
View(r5_three_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R5_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap:
tgfb2_overlap_r5 <- inner_join(r5_three_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R5_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r5 <- inner_join(r5_three_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R5_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r5 <- inner_join(r5_three_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R5_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 
dkk1_overlap_r5 <- inner_join(r5_three_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r5)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R5_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 
bmp4_overlap_r5 <- inner_join(r5_three_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r5)

###############
#R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R6 05 sig/") 
r6_three_occurances_05 <- read.csv("Common DEGs R6 3 I occurances and M6-10 05 sig up and down.csv")
View(r6_three_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R6_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 
tgfb2_overlap_r6 <- inner_join(r6_three_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r6)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R6_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r6 <- inner_join(r6_three_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r6)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R6_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r6 <- inner_join(r6_three_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r6)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R6_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap:
dkk1_overlap_r6 <- inner_join(r6_three_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r6)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R6_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap:
bmp4_overlap_r6 <- inner_join(r6_three_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r6)

##################################################################################
#What Ligand(s) are the four-occurance overlaps from: 0.05 Sig 
##################################################################################
library(dplyr)

###############
#R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R5 05 sig/") 
r5_four_occurances_05 <- read.csv("Common DEGs R5 4 I occurances and M6-10 05 up and down sig.csv")
View(r5_four_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R5_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap:
tgfb2_overlap_r5 <- inner_join(r5_four_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R5_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r5 <- inner_join(r5_four_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R5_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r5 <- inner_join(r5_four_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R5_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 
dkk1_overlap_r5 <- inner_join(r5_four_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r5)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R5_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 
bmp4_overlap_r5 <- inner_join(r5_four_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r5)

###############
#R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R6 05 sig/") 
r6_four_occurances_05 <- read.csv("Common DEGs R6 4 I occurances and M6-10 05 sig up and down.csv")
View(r6_four_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R6_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 
tgfb2_overlap_r6 <- inner_join(r6_four_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r6)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R6_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r6 <- inner_join(r6_four_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r6)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R6_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r6 <- inner_join(r6_four_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r6)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R6_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap:
dkk1_overlap_r6 <- inner_join(r6_four_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r6)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R6_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & abs(bmp4$log2FoldChange) >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap:
bmp4_overlap_r6 <- inner_join(r6_four_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r6)

##################################################################################
#What Ligand(s) are the five-occurance overlaps from: 0.05 Sig 
##################################################################################
library(dplyr)

###############
#R5
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R5 05 sig/") 
r5_five_occurances_05 <- read.csv("Common DEGs R5 5 I occurances and M6-10 05 up and down sig.csv")
View(r5_five_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R5_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap:
tgfb2_overlap_r5 <- inner_join(r5_five_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r5)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R5_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & abs(nlgn1$log2FoldChange) >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r5 <- inner_join(r5_five_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r5)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R5_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & abs(tslp$log2FoldChange) >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r5 <- inner_join(r5_five_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r5)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R5_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & abs(dkk1$log2FoldChange) >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap: 
dkk1_overlap_r5 <- inner_join(r5_five_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r5)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R5_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & bmp4$log2FoldChange >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap: 
bmp4_overlap_r5 <- inner_join(r5_five_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r5)

###############
#R6
###############
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Synergistic Contribution Analysis/Synergistic Analysis Outputs/R6 05 sig/") 
r6_five_occurances_05 <- read.csv("Common DEGs R6 5 I occurances and M6-10 05 sig up and down.csv")
View(r6_four_occurances_05)

#Overlap with TGFBB2 
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tgfb2 <- read.csv("R6_TGFB2.csv")

TGFB2_up_DEGs <- tgfb2[tgfb2$pvalue <= 0.05 & abs(tgfb2$log2FoldChange) >=1,] 
TGFB2_up_nao <- na.omit(TGFB2_up_DEGs) 
colnames(TGFB2_up_nao)[1] <- "Var1"

#Overlap: 
tgfb2_overlap_r6 <- inner_join(r6_five_occurances_05, TGFB2_up_nao, by = "Var1")
View(tgfb2_overlap_r6)

#Overlap with NLGN1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
nlgn1 <- read.csv("R6_NLGN1.csv")

nlgn1_up_DEGs <- nlgn1[nlgn1$pvalue <= 0.05 & nlgn1$log2FoldChange >=1,] 
nlgn1_up_nao <- na.omit(nlgn1_up_DEGs) 
colnames(nlgn1_up_nao)[1] <- "Var1"

#Overlap: 
nlgn1_overlap_r6 <- inner_join(r6_five_occurances_05, nlgn1_up_nao, by = "Var1")
View(nlgn1_overlap_r6)

#Overlap with TSLP
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
tslp <- read.csv("R6_TSLP.csv")

tslp_up_DEGs <- tslp[tslp$pvalue <= 0.05 & tslp$log2FoldChange >=1,] 
tslp_up_nao <- na.omit(tslp_up_DEGs) 
colnames(tslp_up_nao)[1] <- "Var1"

#Overlap: 
tslp_overlap_r6 <- inner_join(r6_five_occurances_05, tslp_up_nao, by = "Var1")
View(tslp_overlap_r6)

#Overlap with DKK1
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
dkk1 <- read.csv("R6_DKK1.csv")

dkk1_up_DEGs <- dkk1[dkk1$pvalue <= 0.05 & dkk1$log2FoldChange >=1,] 
dkk1_up_nao <- na.omit(dkk1_up_DEGs) 
colnames(dkk1_up_nao)[1] <- "Var1"

#Overlap:
dkk1_overlap_r6 <- inner_join(r6_five_occurances_05, dkk1_up_nao, by = "Var1")
View(dkk1_overlap_r6)

#Overlap with BMP4
setwd("~/Desktop/Sloan Lab/Ligand Receptor Project/Rounds 5 and 6 RNA Seq/DESeq Analyses and Volcanos/")
bmp4 <- read.csv("R6_BMP4.csv")

bmp4_up_DEGs <- bmp4[bmp4$pvalue <= 0.05 & bmp4$log2FoldChange >=1,] 
bmp4_up_nao <- na.omit(bmp4_up_DEGs) 
colnames(bmp4_up_nao)[1] <- "Var1"

#Overlap:
bmp4_overlap_r6 <- inner_join(r6_five_occurances_05, bmp4_up_nao, by = "Var1")
View(bmp4_overlap_r6)


