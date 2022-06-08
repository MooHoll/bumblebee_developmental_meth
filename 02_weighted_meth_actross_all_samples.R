#---------------------------------------------------
# Make a file with weighted meth for all genes for all samples for plotting
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth")
library(readr)
library(dplyr)
library(reshape2)
library(data.table)

file.list <- list.files(pattern="*weighted_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("A1","A2","A3","A4","B1","B2","B3","B4","C1","C2","C3", "C4",
                     "D1","D2","D3","D4","E1","E2","F1","F2")
names(samples) <- sample_names


# Add column with methylation ranges for later global methylation plotting
for(i in seq_along(samples)){
  samples[[i]] <- samples[[i]][,c(1,7)]
  colnames(samples[[i]])[2] <- paste0(names(samples[i]))
}

all_data <-  Reduce(function(...) merge(..., by = "geneID", all.x = TRUE),
                samples)
all_data_long <- melt(all_data)
colnames(all_data_long) <- c("geneID","sample","weighted_meth")

all_data_long$meth_category[(all_data_long$weighted_meth >= 0 & all_data_long$weighted_meth <= 0.005)] <- "none"
all_data_long$meth_category[(all_data_long$weighted_meth > 0.005 & all_data_long$weighted_meth <= 0.3)] <- "low"
all_data_long$meth_category[(all_data_long$weighted_meth > 0.3 & all_data_long$weighted_meth <= 0.7)] <- "medium"
all_data_long$meth_category[(all_data_long$weighted_meth > 0.7 & all_data_long$weighted_meth <= 1)] <- "high"

all_data_long$group[all_data_long$sample  %like% "A"] <- "repro_worker_brain"
all_data_long$group[all_data_long$sample  %like% "B"] <- "repro_worker_ovary"
all_data_long$group[all_data_long$sample  %like% "C"] <- "larvae_head"
all_data_long$group[all_data_long$sample  %like% "D"] <- "pupae_head"
all_data_long$group[all_data_long$sample  %like% "E"] <- "male_brain"
all_data_long$group[all_data_long$sample  %like% "F"] <- "male_sperm"

write.table(all_data_long, file="weighted_meth_all_samples.txt", sep="\t",
           quote = F, col.names = T, row.names = F)
