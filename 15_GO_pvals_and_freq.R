# Merging pvalues and frequency information 

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/GO_analysis/Levels_FDR0.05/No_meth")
library(readr)

p_vals <- read_delim("./F_no_meth_against_all_genes_GOs.txt", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE)

freq_info <- read_delim("./F_no_meth_Revigo_BP_Table.tsv", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"

both <- merge(p_vals, freq_info, by = "GOBPID")

write.table(both, file="info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)
