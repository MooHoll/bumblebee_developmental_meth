#---------------------------------------------------
# Background gene sets for differential meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth")
library(readr)
library(data.table)

file.list <- list.files(pattern="*weighted_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("A1","A2","A3","A4","B1","B2","B3","B4","C1","C2","C3", "C4",
                     "D1","D2","D3","D4","E1","E2","F1","F2")
names(samples) <- sample_names

# Pull out lists of methylated genes, i.e. genes >=0.005 
for(i in seq_along(samples)){
  df <- samples[[i]]
  df <- df[,c(1,7)]
  df <- df[df$weightedMeth >= 0.005,] # the weighted meth level of the lambda genome
  myfile <- file.path("./", paste0(names(samples[i]),"_","methylated_genes.txt"))
  write.table(df, file=myfile, quote=F, sep="\t", row.names=F)
}


file.list <- list.files(pattern="*methylated_genes.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("A1","A2","A3","A4","B1","B2","B3","B4","C1","C2","C3", "C4",
                     "D1","D2","D3","D4","E1","E2","F1","F2")
names(samples) <- sample_names

for(i in seq_along(samples)){
  samples[[i]] <- samples[[i]][,1]
}

AvsB <- samples[1:8]
BvsC <- samples[5:12]
CvsD <- samples[9:16]
DvsE <- samples[13:18]
EvsF <- samples[17:20]
BvsF <- samples[c(5:8,19:20)]
AvsE <- samples[c(1:4,17:18)]

AvsB_genes <- Reduce(function(...) merge(..., by = "geneID", all = TRUE),AvsB)
BvsC_genes <- Reduce(function(...) merge(..., by = "geneID", all = TRUE),BvsC)
CvsD_genes <- Reduce(function(...) merge(..., by = "geneID", all = TRUE),CvsD)
DvsE_genes <- Reduce(function(...) merge(..., by = "geneID", all = TRUE),DvsE)
EvsF_genes <- Reduce(function(...) merge(..., by = "geneID", all = TRUE),EvsF)
BvsF_genes <- Reduce(function(...) merge(..., by = "geneID", all = TRUE),BvsF)
AvsE_genes <- Reduce(function(...) merge(..., by = "geneID", all = TRUE),AvsE)

# AvsB 10125
# BvsC 10053
# CvsD 10025
# DvsE 9534
# EvsF 9100
# BvsF 9843
# AvsE 9696

# Now add the annotated GO terms (file made by Alun)
Bumble_bee_ensemble_GO_terms <- read_delim("../GO_analysis/Bumble_bee_ensemble_GO_terms.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
colnames(Bumble_bee_ensemble_GO_terms) <- c("geneID", "GOIds")

AvsB_GO <- merge(AvsB_genes, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(AvsB_GO$geneID))
BvsC_GO <- merge(BvsC_genes, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(BvsC_GO$geneID))
CvsD_GO <- merge(CvsD_genes, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(CvsD_GO$geneID))
DvsE_GO <- merge(DvsE_genes, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(DvsE_GO$geneID))
EvsF_GO <- merge(EvsF_genes, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(EvsF_GO$geneID))
BvsF_GO <- merge(BvsF_genes, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(BvsF_GO$geneID))
AvsE_GO <- merge(AvsE_genes, Bumble_bee_ensemble_GO_terms, by = "geneID")
length(unique(AvsE_GO$geneID))

# AvsB 8204/10125
# BvsC 8155/10053
# CvsD 8125/10025
# DvsE 7747/9534
# EvsF 7529/9100
# BvsF 8021/9843
# AvsE 7898/9696

write.table(AvsB_GO, file="../GO_analysis/AvsB_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(BvsC_GO, file="../GO_analysis/BvsC_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(CvsD_GO, file="../GO_analysis/CvsD_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(DvsE_GO, file="../GO_analysis/DvsE_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(EvsF_GO, file="../GO_analysis/EvsF_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(BvsF_GO, file="../GO_analysis/BvsF_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(AvsE_GO, file="../GO_analysis/AvsE_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Prepare the diff meth gene lists as well from Ben's excel sheet
setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/GO_analysis")

file.list <- list.files(pattern="*diff_meth_genes.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("AvsB","AvsE","BvsC","BvsF","CvsD","DvsE","EvsF")
names(samples) <- sample_names

for(i in seq_along(samples)){
  colnames(samples[[i]]) <- "geneID"
  samples[[i]] <- samples[[i]][!duplicated(samples[[i]]$geneID),]
  samples[[i]] <- samples[[i]][!samples[[i]]$geneID %like% '#',]
  final_file <- samples[[i]]
  print(nrow(final_file))
  myfile <- file.path("./", paste0(names(samples[i]),"_","diff_meth_genes_final.txt"))
  write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}

# Unique diff meth genes
# AvsB 1925
# AvsE 223
# BvsC 585
# BvsF 244
# CvsD 87
# DvsE 631
# EvsF 200


