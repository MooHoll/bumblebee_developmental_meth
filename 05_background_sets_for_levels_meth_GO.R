#---------------------------------------------------
# Background gene sets for levels of meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/weighted_meth_genes")
library(readr)

weighted_meth_by_group <- read_delim("weighted_meth_by_group_genes_only.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
colnames(weighted_meth_by_group) <- c("geneID","group","weighted_meth","meth_category")

Bumble_bee_ensemble_GO_terms <- read_delim("../GO_analysis/Revisions/Bombus_terrestris_HGD_go_annotation.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
Bumble_bee_ensemble_GO_terms <- Bumble_bee_ensemble_GO_terms[,-2]
colnames(Bumble_bee_ensemble_GO_terms) <- c("geneID", "GOIds")

#---------------------------------------------------
# Look at genes with zero methylation compared to all genes?
#---------------------------------------------------

# all genes background set
all_genes_A <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_brain"]))
colnames(all_genes_A) <- "geneID"
all_genes_B <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_ovaries"]))
colnames(all_genes_B) <- "geneID"
all_genes_C <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="larvae"]))
colnames(all_genes_C) <- "geneID"
all_genes_D <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="pupae"]))
colnames(all_genes_D) <- "geneID"
all_genes_E <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="male_brain"]))
colnames(all_genes_E) <- "geneID"
all_genes_F <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="sperm"]))
colnames(all_genes_F) <- "geneID"
# A 10563
# B 10562
# C 10545
# D 10567
# E 10564
# F 10043

all_genes_A_GO <- merge(all_genes_A, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(all_genes_A_GO$geneID))
all_genes_B_GO <- merge(all_genes_B, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(all_genes_B_GO$geneID))
all_genes_C_GO <- merge(all_genes_C, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(all_genes_C_GO$geneID))
all_genes_D_GO <- merge(all_genes_D, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(all_genes_D_GO$geneID))
all_genes_E_GO <- merge(all_genes_E, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(all_genes_E_GO$geneID))
all_genes_F_GO <- merge(all_genes_F, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(all_genes_F_GO$geneID))

# A 9080/10563
# B 9077/10562
# C 9071/10545
# D 9081/10567
# E 9078/10564
# F 8645/10043

write.table(all_genes_A_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/all_genes/A_allgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(all_genes_B_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/all_genes/B_allgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(all_genes_C_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/all_genes/C_allgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(all_genes_D_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/all_genes/D_allgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(all_genes_E_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/all_genes/E_allgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(all_genes_F_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/all_genes/F_allgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Now get the none meth gene lists
none_A <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_brain"
                                                                  & weighted_meth_by_group$meth_category=="none"]))
colnames(none_A) <- "geneID"
none_B <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_ovaries"
                                                                  & weighted_meth_by_group$meth_category=="none"]))
colnames(none_B) <- "geneID"
none_C <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="larvae"
                                                                  & weighted_meth_by_group$meth_category=="none"]))
colnames(none_C) <- "geneID"
none_D <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="pupae"
                                                                  & weighted_meth_by_group$meth_category=="none"]))
colnames(none_D) <- "geneID"
none_E <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="male_brain"
                                                                  & weighted_meth_by_group$meth_category=="none"]))
colnames(none_E) <- "geneID"
none_F <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="sperm"
                                                                  & weighted_meth_by_group$meth_category=="none"]))
colnames(none_F) <- "geneID"
# A 92
# B 46
# C 232
# D 18
# E 87
# F 1169

write.table(none_A, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/A_no_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_B, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/B_no_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_C, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/C_no_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_D, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/D_no_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_E, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/E_no_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_F, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/F_no_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)


#---------------------------------------------------
# Look at highly methylated genes compared to all methylated genes?
#---------------------------------------------------

# all meth gene lists 
meth_A <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_brain"
                                                             & weighted_meth_by_group$weighted_meth> 0.05]))
colnames(meth_A) <- "geneID"
meth_B <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_ovaries"
                                                             & weighted_meth_by_group$weighted_meth> 0.05]))
colnames(meth_B) <- "geneID"
meth_C <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="larvae"
                                                             & weighted_meth_by_group$weighted_meth> 0.05]))
colnames(meth_C) <- "geneID"
meth_D <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="pupae"
                                                             & weighted_meth_by_group$weighted_meth> 0.05]))
colnames(meth_D) <- "geneID"
meth_E <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="male_brain"
                                                             & weighted_meth_by_group$weighted_meth> 0.05]))
colnames(meth_E) <- "geneID"
meth_F <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="sperm"
                                                             & weighted_meth_by_group$weighted_meth> 0.05]))
colnames(meth_F) <- "geneID"

# A 4706
# B 3893
# C 4466
# D 4710
# E 4773
# F 3705

meth_A_GO <- merge(meth_A, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(meth_A_GO$geneID))
meth_B_GO <- merge(meth_B, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(meth_B_GO$geneID))
meth_C_GO <- merge(meth_C, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(meth_C_GO$geneID))
meth_D_GO <- merge(meth_D, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(meth_D_GO$geneID))
meth_E_GO <- merge(meth_E, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(meth_E_GO$geneID))
meth_F_GO <- merge(meth_F, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(meth_F_GO$geneID))

# A 4458/4706
# B 3693/3893
# C 4231/4466
# D 4458/4710
# E 4522/4773
# F 3524/3705

write.table(meth_A_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/A_methgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(meth_B_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/B_methgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(meth_C_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/C_methgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(meth_D_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/D_methgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(meth_E_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/E_methgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(meth_F_GO, file="../GO_analysis/Revisions/background_levels_GO_lists/meth_genes/F_methgenes_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)


# Highly meth gene lists
high_A <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_brain"
                                                             & weighted_meth_by_group$meth_category=="high"]))
colnames(high_A) <- "geneID"
high_B <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_ovaries"
                                                             & weighted_meth_by_group$meth_category=="high"]))
colnames(high_B) <- "geneID"
high_C <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="larvae"
                                                             & weighted_meth_by_group$meth_category=="high"]))
colnames(high_C) <- "geneID"
high_D <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="pupae"
                                                             & weighted_meth_by_group$meth_category=="high"]))
colnames(high_D) <- "geneID"
high_E <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="male_brain"
                                                             & weighted_meth_by_group$meth_category=="high"]))
colnames(high_E) <- "geneID"
high_F <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="sperm"
                                                             & weighted_meth_by_group$meth_category=="high"]))
colnames(high_F) <- "geneID"
# A 35
# B 1
# C 42
# D 9
# E 45
# F 450

write.table(high_A, file="../GO_analysis/A_high_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_B, file="../GO_analysis/B_high_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_C, file="../GO_analysis/C_high_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_D, file="../GO_analysis/D_high_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_E, file="../GO_analysis/E_high_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_F, file="../GO_analysis/F_high_methylation_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)
