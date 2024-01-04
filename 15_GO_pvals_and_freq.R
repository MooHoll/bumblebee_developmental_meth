# Merging pvalues and frequency information 

# Diff meth genes
setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/Diff_meth_FDR0.05")
library(readr)

p_vals <- read_delim("AvsB_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("AvsB_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="AvsB_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("AvsE_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("AvsE_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="AvsE_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("BvsC_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("BvsC_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="BvsC_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("BvsF_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("BvsF_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="BvsF_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("CvsD_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("CvsD_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="CvsD_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("DvsE_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("DvsE_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="DvsE_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("EvsF_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("EvsF_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="EvsF_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)


# JSD index
setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/JSD_index")

p_vals <- read_delim("outliers_vs_all_genes_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("outlier_vs_all_genes_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="outlier_vs_all_genes_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("weird_chr_vs_all_outliers_enriched_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("weird_chr_vs_all_outliers_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="weird_chr_vs_all_outliers_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)


# High meth genes
setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/Levels_FDR0.05/High_meth")

p_vals <- read_delim("A_high_meth_against_all_meth_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("A_high_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="A_high_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("C_high_meth_against_all_meth_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("C_high_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="C_high_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("D_high_meth_against_all_meth_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("D_high_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="D_high_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("E_high_meth_against_all_meth_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("E_high_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="E_high_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("F_high_meth_against_all_meth_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("F_high_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="F_high_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)


# No meth genes
setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/Levels_FDR0.05/No_meth")

p_vals <- read_delim("A_no_meth_against_all_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("A_no_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="A_no_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("B_no_meth_against_all_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("B_no_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="B_no_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("C_no_meth_against_all_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("C_no_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="C_no_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("D_no_meth_against_all_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("D_no_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="D_no_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("E_no_meth_against_all_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("E_no_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="E_no_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)

p_vals <- read_delim("F_no_meth_against_all_genes_GOs.txt",  delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
freq_info <- read_delim("F_no_Revigo_BP_Table.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"
both <- merge(p_vals, freq_info, by = "GOBPID")
write.table(both, file="F_no_info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)


