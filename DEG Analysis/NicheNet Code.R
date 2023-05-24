##Anna Voss | Sloan Lab 

##Adapted from NicheNet Vinette: https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md


##Step 0 
setwd("~/aligned_files/bam_files/featurecounts")
library(nichenetr)
library(tidyverse)
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds")) 
x <- read.delim("feature_counts2_counts_only.txt")
x2 <- x[,-1]
rownames(x2) <- make.names(x[,1], unique=TRUE)

##Step 1 Adapted from https://stackoverflow.com/questions/37758738/r-remove-rows-with-fewer-than-certain-threshold-non-zero-values
#Receiver Cell
astros <- x2[,1:2]
drop2 <- NULL

for (i in 1:nrow(astros)){
  sums2 <- sum(astros[i,])
  drop2 <- c(drop2, sums2 <= 500)
}
astros_drop_500 <- astros[!drop2==TRUE,]

nrow(astros_drop_500)
astro_names <- rownames(astros_drop_500)
astro_names_upper <- toupper(astro_names)
head(astro_names_upper)

#Repeat for Sender Cell
#Neuron
neuron <- x2[,3:4]
drop <- NULL

for (i in 1:nrow(neuron)){
  sums <- sum(neuron[i,])
  drop <- c(drop, sums <= 500)
}
neuron_drop_500 <- neuron[!drop==TRUE,]

nrow(neuron_drop_500)
neuron_names <- rownames(neuron_drop_500)
neuron_names_upper <- toupper(neuron_names)

#OPC
opc <- x2[,5:6]
drop <- NULL

for (i in 1:nrow(opc)){
  sums <- sum(opc[i,])
  drop <- c(drop, sums <= 500)
}
opc_drop_500 <- opc[!drop==TRUE,]

nrow(opc_drop_500)
opc_names <- rownames(opc_drop_500)
opc_names_upper <- toupper(opc_names)

#ONF
onf <- x2[,7:8]
drop <- NULL

for (i in 1:nrow(onf)){
  sums <- sum(onf[i,])
  drop <- c(drop, sums <= 500)
}
onf_drop_500 <- onf[!drop==TRUE,]

nrow(onf_drop_500)
onf_names <- rownames(onf_drop_500)
onf_names_upper <- toupper(onf_names)

#OM
om <- x2[,9:10]
drop <- NULL

for (i in 1:nrow(om)){
  sums <- sum(om[i,])
  drop <- c(drop, sums <= 500)
}
om_drop_500 <- om[!drop==TRUE,]

nrow(om_drop_500)
om_names <- rownames(om_drop_500)
om_names_upper <- toupper(om_names)

#Microglia
micro <- x2[,11:12]
drop <- NULL

for (i in 1:nrow(micro)){
  sums <- sum(micro[i,])
  drop <- c(drop, sums <= 500)
}
micro_drop_500 <- micro[!drop==TRUE,]

nrow(micro_drop_500)
micro_names <- rownames(micro_drop_500)
micro_names_upper <- toupper(micro_names)

#Endo
endo <- x2[,13:14]
drop <- NULL

for (i in 1:nrow(endo)){
  sums <- sum(endo[i,])
  drop <- c(drop, sums <= 500)
}
endo_drop_500 <- endo[!drop==TRUE,]

nrow(endo_drop_500)
endo_names <- rownames(endo_drop_500)
endo_names_upper <- toupper(endo_names)

#Germany Receiving 
osvz <- x2[,13:18]
drop <- NULL

for (i in 1:nrow(osvz)){
  sums <- sum(osvz[i,])
  drop <- c(drop, sums <= 500)
}
osvz_drop_500 <- osvz[!drop==TRUE,]

nrow(osvz_drop_500)
osvz_names <- rownames(osvz_drop_500)
osvz_names_upper <- toupper(osvz_names)

#Germany Sending
vz <- x2[,19:24]
drop <- NULL

for (i in 1:nrow(vz)){
  sums <- sum(vz[i,])
  drop <- c(drop, sums <= 500)
}
vz_drop_500 <- vz[!drop==TRUE,]

nrow(vz_drop_500)
vz_names <- rownames(vz_drop_500)
vz_names_upper <- toupper(vz_names)

##Step 2 
geneset_oi <- read.csv("immature_astro_genes.csv", header = FALSE) #First Excel Sheet (Immature Only) 
geneset_oi <- read.csv("geneset_full.csv", header = FALSE) #Both Excel Sheets (Immature and Mature)  
geneset_oi <- read.csv("mature_astro_genes.csv", header = FALSE) #Second Excel Sheet (Mature Only) 

rownames(geneset_oi) <- make.names(geneset_oi[,1], unique = TRUE)
geneset_oi_names <- rownames(geneset_oi)

background_expressed_genes = osvz_names_upper %>% .[. %in% rownames(ligand_target_matrix)]


##Step 3
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands, endo_names_upper)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors, astro_names_upper)

lr_network_expressed = lr_network %>% filter(from  %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

##Step 4
ligand_activities = predict_ligand_activities(geneset = geneset_oi_names, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities %>% arrange(-pearson)

best_upstream_ligands = ligand_activities %>% top_n(50, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange", binwidth = 0.01)  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

##Step 5
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi_names, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)


order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
vis_ligand_target = active_ligand_target_links[,order_ligands] %>% t() #different than the documentation- the length of order targets is 160 instead of 159 like in active_ligand_target_links

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Ventricular Zone Ligands","Mature Astrocyte Genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

##Follow up Analysis 1
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized Ventricular Zone Ligands","Receptors expressed in the oSVZ", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network
