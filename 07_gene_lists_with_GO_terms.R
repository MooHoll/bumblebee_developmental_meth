#---------------------------------------------------
# Gene lists with GO terms
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/GO_analysis")
library(readr)

# Main files
Bumble_bee_ensemble_GO_terms <- read_delim("../Bumble_bee_ensemble_GO_terms.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
colnames(Bumble_bee_ensemble_GO_terms) <- c("geneID", "GOIds")

weighted_meth_by_group <- read_delim("../../weighted_meth_genes/weighted_meth_by_group_genes_only.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
colnames(weighted_meth_by_group)<-c("geneID","group","weighted_methylation","meth_category")

# Get gene lists for no methylation
none_A <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_brain" & weighted_meth_by_group$meth_category=="none"]))
none_B <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_ovaries"  & weighted_meth_by_group$meth_category=="none"]))
none_C <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="larvae" & weighted_meth_by_group$meth_category=="none"]))
none_D <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="pupae" & weighted_meth_by_group$meth_category=="none"]))
none_E <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="male_brain" & weighted_meth_by_group$meth_category=="none"]))
none_F <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="sperm"  & weighted_meth_by_group$meth_category=="none"]))

colnames(none_A) <- "geneID"
colnames(none_B) <- "geneID"
colnames(none_C) <- "geneID"
colnames(none_D) <- "geneID"
colnames(none_E) <- "geneID"
colnames(none_F) <- "geneID"

none_A_GO <- merge(none_A, Bumble_bee_ensemble_GO_terms, by = "geneID")
none_B_GO <- merge(none_B, Bumble_bee_ensemble_GO_terms, by = "geneID")
none_C_GO <- merge(none_C, Bumble_bee_ensemble_GO_terms, by = "geneID")
none_D_GO <- merge(none_D, Bumble_bee_ensemble_GO_terms, by = "geneID")
none_E_GO <- merge(none_E, Bumble_bee_ensemble_GO_terms, by = "geneID")
none_F_GO <- merge(none_F, Bumble_bee_ensemble_GO_terms, by = "geneID")

write.table(none_A_GO, file="A_no_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_B_GO, file="B_no_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_C_GO, file="C_no_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_D_GO, file="D_no_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_E_GO, file="E_no_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(none_F_GO, file="F_no_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)


# Get gene lists for high methylation
high_A <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_brain" & weighted_meth_by_group$meth_category=="high"]))
high_B <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="repro_ovaries"  & weighted_meth_by_group$meth_category=="high"]))
high_C <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="larvae" & weighted_meth_by_group$meth_category=="high"]))
high_D <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="pupae" & weighted_meth_by_group$meth_category=="high"]))
high_E <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="male_brain" & weighted_meth_by_group$meth_category=="high"]))
high_F <- as.data.frame(unique(weighted_meth_by_group$geneID[weighted_meth_by_group$group=="sperm"  & weighted_meth_by_group$meth_category=="high"]))

colnames(high_A) <- "geneID"
colnames(high_B) <- "geneID"
colnames(high_C) <- "geneID"
colnames(high_D) <- "geneID"
colnames(high_E) <- "geneID"
colnames(high_F) <- "geneID"

high_A_GO <- merge(high_A, Bumble_bee_ensemble_GO_terms, by = "geneID")
high_B_GO <- merge(high_B, Bumble_bee_ensemble_GO_terms, by = "geneID")
high_C_GO <- merge(high_C, Bumble_bee_ensemble_GO_terms, by = "geneID")
high_D_GO <- merge(high_D, Bumble_bee_ensemble_GO_terms, by = "geneID")
high_E_GO <- merge(high_E, Bumble_bee_ensemble_GO_terms, by = "geneID")
high_F_GO <- merge(high_F, Bumble_bee_ensemble_GO_terms, by = "geneID")

write.table(high_A_GO, file="A_high_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_B_GO, file="B_high_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_C_GO, file="C_high_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_D_GO, file="D_high_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_E_GO, file="E_high_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(high_F_GO, file="F_high_meth_genes_with_GOs.txt", sep="\t",quote = F, col.names = T, row.names = F)


# Get gene lists for diff meth
setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/GO_analysis/diff_meth_genes")
file.list <- list.files(pattern="*diff_meth_genes_final.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("AvsB","AvsE","BvsC","BvsF","CvsD","DvsE","EvsF")
names(samples) <- sample_names

for(i in seq_along(samples)){
  df <- samples[[i]]
  final_file <- merge(df, Bumble_bee_ensemble_GO_terms, by = "gene_id")
  myfile <- file.path("./", paste0(names(samples[i]),"_","with_GOs.txt"))
  write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F, col.names = T)
}
