# Make weighted meth levels for bew TE annotations

# Making the file which has the total CpGs per gene information
setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/TEs")
library(sqldf)
library(readr)
library(doBy)

cpgs <- read.delim("total_cpgs_in_genome.txt", header=F)
colnames(cpgs) <- c("scaffold", "cpg_position")
cpgs$cpg_position <- as.numeric(cpgs$cpg_position)

genes <- read.delim("te_annotation.txt", header=F)
colnames(genes)<- c("scaffold","feature","gene_id","start","end")
genes$start <- as.numeric(genes$start)
genes$end <- as.numeric(genes$end)

output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.feature,
                genes.start,
                genes.end,
                genes.gene_id
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$gene_id),]

output$cpg_counter <- 1
final <- summaryBy(cpg_counter ~ gene_id+scaffold+start+end+feature, data = output, FUN=sum) 
write.table(final, file="TEs_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)


# Calculate meth level for each gene for each sample (worth submitting as job - R/4.3.1)
library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*final_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("A1","A2","A3","A4","B1","B2","B3","B4","C1","C2","C3", "C4",
                     "D1","D2","D3","D4","E1","E2","F1","F2")
names(samples) <- sample_names

# Read in gene with start/end and total CpGs per gene
annotation_with_total_cpgs <- read_table2("TEs_with_total_cpgs.txt")
colnames(annotation_with_total_cpgs)[6] <- "cpg_count"
colnames(annotation_with_total_cpgs)[2] <- "chr"

registerDoParallel(cores = 8)

# Calculate weighted meth for each gene for each sample
foreach(i = seq_along(samples)) %dopar% {
  df <- samples[[i]]
  df <- subset(df, total_coverage > 5)
  output <- sqldf("SELECT sample.chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.chr,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.feature,
                    annot.cpg_count
                    FROM df AS sample
                    LEFT JOIN annotation_with_total_cpgs AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$gene_id),]
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ chr + feature + gene_id + start + end + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$count_c.sum)/(check$total_coverage.sum)
  myfile <- file.path("./", paste0(names(samples[i]),"_","weighted_meth_TEs.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}
