## -------------------------------------------------------------------------
# Take average weighted methylation level of feature across bio replicates
## -------------------------------------------------------------------------
setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_genes")
library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(reshape2)

# Make one file covering all samples
file.list = list.files(("./"),pattern="*genes.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)


# Make one dataframe for each population
A <- samples[1:4]
B <- samples[5:8]
C <- samples[9:12]
D <- samples[13:16]
E <- samples[17:18]
f <- samples[19:20]

for(i in seq_along(A)){
  A[[i]]$stage <- "repro_brain"
}
A_all <- as.data.frame(bind_rows(A))
A_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                         stage, data = A_all, FUN=mean)

for(i in seq_along(B)){
  B[[i]]$stage <- "repro_ovaries"
}
B_all <- as.data.frame(bind_rows(B))
B_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                         stage, data = B_all, FUN=mean)

for(i in seq_along(C)){
  C[[i]]$stage <- "larvae"
}
C_all <- as.data.frame(bind_rows(C))
C_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                         stage, data = C_all, FUN=mean)

for(i in seq_along(D)){
  D[[i]]$stage <- "pupae"
}
D_all <- as.data.frame(bind_rows(D))
D_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                         stage, data = D_all, FUN=mean)

for(i in seq_along(E)){
  E[[i]]$stage <- "male_brain"
}
E_all <- as.data.frame(bind_rows(E))
E_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                         stage, data = E_all, FUN=mean)

for(i in seq_along(f)){
  f[[i]]$stage <- "sperm"
}
f_all <- as.data.frame(bind_rows(f))
f_merged <- summaryBy(weightedMeth ~ feature + gene_id + start + end +
                        stage, data = f_all, FUN=mean)


all_data <- rbind(A_merged, B_merged, C_merged, D_merged, E_merged, f_merged)

all_data2 <- dcast(all_data, feature + gene_id + start + end  ~ stage, value.var = "weightedMeth.mean")

write.table(all_data2, file="weighted_meth_annotation_by_stage_genes_only.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")



