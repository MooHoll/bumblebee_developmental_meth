#---------------------------------------------------
# Background gene sets for differential meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_genes")
library(readr)
library(data.table)

file.list <- list.files(pattern="*just_genes.txt")

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
  df <- df[,c(2,8)]
  df <- df[df$weightedMeth >= 0.05,] # the weighted meth level of the lambda genome
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

AvsB_genes <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE),AvsB)
BvsC_genes <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE),BvsC)
CvsD_genes <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE),CvsD)
DvsE_genes <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE),DvsE)
EvsF_genes <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE),EvsF)
BvsF_genes <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE),BvsF)
AvsE_genes <- Reduce(function(...) merge(..., by = "gene_id", all = TRUE),AvsE)

# AvsB 5083
# BvsC 5052
# CvsD 5130
# DvsE 5069
# EvsF 5142
# BvsF 4983
# AvsE 5112

# Now add the annotated GO terms (file made by Alun)
Bumble_bee_ensemble_GO_terms <- read_delim("../GO_analysis/Bumble_bee_ensemble_GO_terms.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
colnames(Bumble_bee_ensemble_GO_terms) <- c("gene_id", "GOIds")

AvsB_GO <- merge(AvsB_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(AvsB_GO$gene_id))
BvsC_GO <- merge(BvsC_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(BvsC_GO$gene_id))
CvsD_GO <- merge(CvsD_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(CvsD_GO$gene_id))
DvsE_GO <- merge(DvsE_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(DvsE_GO$gene_id))
EvsF_GO <- merge(EvsF_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(EvsF_GO$gene_id))
BvsF_GO <- merge(BvsF_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(BvsF_GO$gene_id))
AvsE_GO <- merge(AvsE_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(AvsE_GO$gene_id))

# AvsB 4715/5083
# BvsC 4685/5052
# CvsD 4758/5130
# DvsE 4716/5069
# EvsF 4774/5142
# BvsF 4621/4983
# AvsE 4752/5112

write.table(AvsB_GO, file="../GO_analysis/background_GO_lists/AvsB_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(BvsC_GO, file="../GO_analysis/background_GO_lists/BvsC_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(CvsD_GO, file="../GO_analysis/background_GO_lists/CvsD_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(DvsE_GO, file="../GO_analysis/background_GO_lists/DvsE_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(EvsF_GO, file="../GO_analysis/background_GO_lists/EvsF_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(BvsF_GO, file="../GO_analysis/background_GO_lists/BvsF_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)
write.table(AvsE_GO, file="../GO_analysis/background_GO_lists/AvsE_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)

# Prepare the diff meth gene lists as well from Ben's excel sheet
setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists")

file.list <- list.files(pattern="*diff_meth_genes.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("AvsB","AvsE","BvsC","BvsF","CvsD","DvsE","EvsF")
names(samples) <- sample_names

for(i in seq_along(samples)){
  samples[[i]] <- samples[[i]][!duplicated(samples[[i]]$gene_id),]
  samples[[i]] <- samples[[i]][!samples[[i]]$gene_id %like% '#',]
  final_file <- samples[[i]][,1]
  print(nrow(final_file))
  myfile <- file.path("./", paste0(names(samples[i]),"_","diff_meth_genes_final.txt"))
  write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}

# Unique diff meth genes
# AvsB 1222
# AvsE 37
# BvsC 286
# BvsF 134
# CvsD 5
# DvsE 156
# EvsF 49


