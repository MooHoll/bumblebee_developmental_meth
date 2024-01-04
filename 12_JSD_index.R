# ----------------------------------------------------------------
# Calculating the Jensen-Shannon diversity index for all genes
# ----------------------------------------------------------------

setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/weighted_meth_genes")

# Useful links
# https://github.com/genomaths/MethylIT
# https://rdrr.io/github/genomaths/MethylIT.utils/man/jensenSDiv.html

#install.packages("devtools")
#devtools::install_git("https://github.com/genomaths/MethylIT.utils.git")

library(MethylIT)
library(MethylIT.utils)
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)

# ----------------------------------------------------------------

# jensenSDiv(p, q, Pi = 0.5, logbase = 2)
# p and q are probability vectors based on counts/sum

# Once have a JSD value for each comparison per gene can average them and 
# then see which genes are the most variable across all of development

# ----------------------------------------------------------------

# I think I was right the first time to use the average for each dev stage, it gives a JSD index per pair of datapoints
weighted_meth_annotation_by_stage <- read_delim("weighted_meth_annotation_by_stage_genes_only.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)
weighted_meth_annotation_by_stage <- weighted_meth_annotation_by_stage[weighted_meth_annotation_by_stage$feature=="gene",]
head(weighted_meth_annotation_by_stage)
weighted_meth_annotation_by_stage <- weighted_meth_annotation_by_stage[,c(2,5:10)]
weighted_meth_annotation_by_stage <- weighted_meth_annotation_by_stage[complete.cases(weighted_meth_annotation_by_stage),]

# Only keep genes classed as methylated in at least one developmental stage
weighted_meth_annotation_by_stage <- weighted_meth_annotation_by_stage[weighted_meth_annotation_by_stage$larvae > 0.005 |
                                                                         weighted_meth_annotation_by_stage$male_brain > 0.005 |
                                                                         weighted_meth_annotation_by_stage$pupae > 0.005 |
                                                                         weighted_meth_annotation_by_stage$repro_brain > 0.005 |
                                                                         weighted_meth_annotation_by_stage$repro_ovaries > 0.005 |
                                                                         weighted_meth_annotation_by_stage$sperm > 0.005,]

weighted_meth_annotation_by_stage$JSD_AvsB <- jensenSDiv(weighted_meth_annotation_by_stage$repro_brain, weighted_meth_annotation_by_stage$repro_ovaries, Pi = 0.5, logbase = 2)
weighted_meth_annotation_by_stage$JSD_BvsC <- jensenSDiv(weighted_meth_annotation_by_stage$repro_ovaries, weighted_meth_annotation_by_stage$larvae, Pi = 0.5, logbase = 2)
weighted_meth_annotation_by_stage$JSD_CvsD <- jensenSDiv(weighted_meth_annotation_by_stage$larvae, weighted_meth_annotation_by_stage$pupae, Pi = 0.5, logbase = 2)
weighted_meth_annotation_by_stage$JSD_DvsE <- jensenSDiv(weighted_meth_annotation_by_stage$pupae, weighted_meth_annotation_by_stage$male_brain, Pi = 0.5, logbase = 2)
weighted_meth_annotation_by_stage$JSD_EvsF <- jensenSDiv(weighted_meth_annotation_by_stage$male_brain, weighted_meth_annotation_by_stage$sperm, Pi = 0.5, logbase = 2)

weighted_meth_annotation_by_stage$mean_JSD <- rowMeans(weighted_meth_annotation_by_stage[,8:12], na.rm=TRUE)

for_table <- weighted_meth_annotation_by_stage[,c(1,13)] # 8182 genes

# Let's see what's going on then
plot(for_table$mean_JSD)

# Wow so most values are below 0.001
hist(for_table$mean_JSD[for_table$mean_JSD < 0.001], breaks=50)
look<- boxplot(for_table$mean_JSD)
options(scipen = 999) # remove e values numbers
look[["stats"]]
# 0.037 is the upper whisker, so ok everything above is an outlier
# Remember this is > 1.5X the interquartile range

ggplot(for_table, aes(y=mean_JSD))+
  geom_boxplot()+
  theme_bw()+
  ylab("Mean JSD index")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20))

for_table$outlier <- ifelse(for_table$mean_JSD >= 0.037, "yes", "no") #827 yes, 7355 no
write.table(for_table, file = "genes_with_JSDindex_new_annotation.txt", sep="\t", quote = F, col.names = T, row.names = F)

# ----------------------------------------------------------------
# Hmmm add in chromosome numbers?
genes_new_annotation <- read_delim("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/genes_new_annotation.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   col_names = FALSE, trim_ws = TRUE)
chroms <- genes_new_annotation[,c(1,6)]
colnames(chroms)<-c("chr","gene_id")
chroms$chr[chroms$chr %like% "NW"] <- "unplaced_scaffolds"

JSD_with_chrs <- merge(for_table, chroms, by="gene_id")

ggplot(JSD_with_chrs, aes(x=chr, y= mean_JSD, fill=chr))+
  geom_boxplot()+
  theme_bw()+
  xlab("Chromosome")+
  ylab("Mean JSD index")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title=element_text(size=20),
        legend.position = "none")

# Are there certain chromosomes with more variable meth genes?
# Need to normalise by gene number
number_genes <- as.data.frame(table(chroms$chr))
colnames(number_genes) <- c("chr","num_genes")

outliers <- JSD_with_chrs[JSD_with_chrs$outlier=="yes",] #827
outliers_count <- as.data.frame(table(outliers$chr))
colnames(outliers_count) <- c("chr","outliers")

outliers_for_plot <- merge(outliers_count, number_genes, by="chr")
outliers_for_plot$proportion <- outliers_for_plot$num_genes / outliers_for_plot$outliers

ggplot(outliers_for_plot, aes(x=chr, y=proportion, fill=chr))+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("Chromosome")+
  ylab("Proportion of genes with >0.024 JSD index")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title=element_text(size=20),
        legend.position = "none")

# ----------------------------------------------------------------
# Pull out gene lists for GO enrichment
all_genes <- as.data.frame(for_table$gene_id) #8182
colnames(all_genes) <- "gene_id"
write.table(all_genes, file="all_genes_JSD_GO_background.txt", sep="\t", quote = F, col.names = T, row.names = F)

outliers_for_go <- as.data.frame(outliers$gene_id) #827
colnames(outliers_for_go) <- "gene_id"
write.table(outliers_for_go, file="outliers_JSD_for_GO.txt", sep="\t", quote = F, col.names = T, row.names = F)

weird_chr_genes <- as.data.frame(outliers$gene_id[outliers$chr=="NC_015774.1"]) #12
colnames(weird_chr_genes) <- "gene_id"
write.table(weird_chr_genes, file="weird_chr_JSD_for_GO.txt", sep="\t", quote = F, col.names = T, row.names = F)

# ----------------------------------------------------------------
# Make GO background files

Bumble_bee_ensemble_GO_terms <- read_delim("../GO_analysis/Revisions/Bombus_terrestris_HGD_go_annotation.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
Bumble_bee_ensemble_GO_terms <- Bumble_bee_ensemble_GO_terms[,-2]
colnames(Bumble_bee_ensemble_GO_terms) <- c("gene_id", "GOIds")

library(readr)
setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/JSD_index")

all_genes <- read_csv("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/JSD_index/all_genes_JSD_GO_background.txt")
all_genes_GO <- merge(all_genes, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(all_genes_GO$gene_id)) # 7094/8182
write.table(all_genes_GO, file="all_genes_JSD_with_GOterms.txt", sep="\t", quote = F, col.names = T, row.names = F)

outliers_for_go <- read_csv("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/JSD_index/outliers_JSD_for_GO.txt")
outliers_GO <- merge(outliers_for_go, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(outliers_GO$gene_id)) # 777/827
write.table(outliers_GO, file="outlier_JSD_with_GOterms.txt", sep="\t", quote = F, col.names = T, row.names = F)


# ----------------------------------------------------------------
## UNUSED

# Read in weighted meth levels for all developmental stages and feature
file.list = list.files(("./"),pattern="*features.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

# Add the sample name as the dataframe name
base <- tools::file_path_sans_ext(basename(file.list))
base <- gsub("_weighted_meth_all_features","",base)
names(samples) <- base

# Look at only genes, ensure min cov 10, remove unnecessary cols, add sample name 
for(i in seq_along(samples)){
  samples[[i]] <- samples[[i]][samples[[i]]$feature=="gene",]
  samples[[i]] <- samples[[i]][samples[[i]]$total_coverage.sum > 10,]
  samples[[i]] <- samples[[i]][,c(2,9)]
  samples[[i]]$name <- paste0(names(samples[i]))
}

# Put all the data together
all_data <- as.data.frame(bind_rows(samples))

# Convert to wide format
all_data_wide <- dcast(all_data, gene_id ~ name, value.var="weightedMeth")

# TEST
test <- head(all_data_wide[,c(1:9)], n=3)
test$JSD <- NA

for(i in 1:nrow(test)) {
  p <- test[i,c(2:5)]
  q <- test[i,c(6:9)]
  test$JSD <- jensenSDiv(p, q, Pi = 0.5, logbase = 2)
}

p = c(0.013748854,0.003787879,0.027722772,0.020373514)
q = c(0.007587253,0.014618974,0.023444545,0.020491803)

jensenSDiv(p, q, Pi = 0.5, logbase = 2)