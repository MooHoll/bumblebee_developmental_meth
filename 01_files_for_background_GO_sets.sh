###-------------------------------------------------
# Make background sets for GO enrichment
###-------------------------------------------------

# --- Bash ---
# Use one genome-wide cytosine report from methylation_extraction from bismark
# to create a file which consists of just the scaffold name and the CpG position
# Only take + strand coordinate otherwise get each CpG counted twice
grep "+" A1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt > plus_only.txt
cut -f1,2 plus_only.txt > total_cpgs_in_genome.txt

# Gene coordinate file
/scratch/monoallelic/bh204/static/bter_genome/gene_with_start_end.txt

# --- R ---
module load R/3.6.1

# Making the file which has the total CpGs per gene information
library(sqldf)
library(readr)
library(doBy)

cpgs <- read.delim("total_cpgs_in_genome.txt", header=F)
colnames(cpgs) <- c("scaffold", "cpg_position")
cpgs$cpg_position <- as.numeric(cpgs$cpg_position)

genes <- read.delim("gene_with_start_end.txt", header=T)
colnames(genes)[1] <- "scaffold"
genes$start <- as.numeric(genes$start)
genes$end <- as.numeric(genes$end)

output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.geneID
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$geneID),]

output$cpg_counter <- 1
final <- summaryBy(cpg_counter ~ geneID+scaffold+start+end, data = output, FUN=sum) #11,023 genes
write.table(final, file="Bter1.0_genes_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)


# --- Bash ---
# Edit coverage files to calculate meth per gene levels
for f in /data/monoallelic/bh204/dmd_final/data/extract/*/*cov.gz 
do 
    cp $f /scratch/monoallelic/hjm32/ben_development
done

gunzip *cov.gz

for file in $(ls *cov)
do
    base=$(basename ${file} "_val_1_bismark_bt2_pe.deduplicated.bismark.cov")
    cut -f1,2,5,6 ${file} > ${base}_coverage.txt
done


# --- R ---
module load R/3.6.1

library(readr)

file.list = list.files(("./"),pattern="*_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", col_names=F)
}

samples <- lapply(file.list, read_file1)

sample_names <- list("A1","A2","A3","A4","B1","B2","B3","B4","C1","C2","C3", "C4",
                     "D1","D2","D3","D4","E1","E2","F1","F2")
names(samples) <- sample_names

for(i in seq_along(samples)){
    colnames(samples[[i]]) <- c("chr", "cpg", "count_c", "count_t")
    samples[[i]]$total_coverage <- samples[[i]]$count_c + samples[[i]]$count_t
    samples[[i]] <- samples[[i]][,-4]
    final_file <- samples[[i]]
    myfile <- file.path("./", paste0(names(samples[i]),"_","final_coverage.txt"))
    write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}


# --- R ---
module load R/3.6.1

# Calculate meth level for each gene for each sample (worth submitting as job)
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
annotation_with_total_cpgs <- read_table2("../Bter1.0_genes_with_total_cpgs.txt")
colnames(annotation_with_total_cpgs)[5] <- "cpg_count"

registerDoParallel(cores = 10)

# Calculate weighted meth for each gene for each sample
foreach(i = seq_along(samples)) %dopar% {
  df <- samples[[i]]
  df <- subset(df, total_coverage > 5)
  output <- sqldf("SELECT sample.chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.scaffold,
                    annot.start,
                    annot.end,
                    annot.geneID,
                    annot.cpg_count
                    FROM df AS sample
                    LEFT JOIN annotation_with_total_cpgs AS annot
                    ON sample.chr = annot.scaffold
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$geneID),]
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ chr + feature + geneID + start + end + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$cpg_count*check$count_c.sum)/(check$cpg_count*check$total_coverage.sum)
  myfile <- file.path("./", paste0(names(samples[i]),"_","weighted_meth.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}
