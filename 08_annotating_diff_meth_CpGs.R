## -------------------------------------------------------------------------
# Filtering diff meth CpGs from methylkit with weighted meth of features
## -------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists")
library(sqldf)
library(readr)
library(doBy)
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)

## -------------------------------------------------------------------------
# Get files in order
## -------------------------------------------------------------------------

# Annotation file
annotation <- read_delim("../GCF_000214255.1_Bter_1.0_genomic_numbered_exons.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

# Remove stuff that isn't really informative
annotation <- annotation[!annotation$feature == "CDS",]
annotation <- annotation[!annotation$feature == "mRNA",]
annotation <- annotation[!annotation$feature == "transcript",]
annotation <- annotation[!annotation$feature == "RNA",]

# Remove gene as this is already accounted for by exon/intron etc
annotation <- annotation[!annotation$feature == "gene",]

# All of the diff meth CpG lists
file.list <- list.files(pattern="*diff_CpGs.txt")
read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("AvsB","AvsE","BvsC","BvsF","CvsD","DvsE","EvsF")
names(samples) <- sample_names

# Add sample names and colnames
for(i in seq_along(samples)){
  samples[[i]]$comparison <- paste0(names(samples[i]))
  colnames(samples[[i]]) <- c("chr","cpg_position","methylation_diff","comparison")
}

# Add column for hypermeth sex
samples[[1]]$hypermethylated <- "A"
samples[[1]]$hypermethylated[samples[[1]]$methylation_diff > 0] <- "B"
samples[[2]]$hypermethylated <- "A"
samples[[2]]$hypermethylated[samples[[2]]$methylation_diff > 0] <- "E"
samples[[3]]$hypermethylated <- "B"
samples[[3]]$hypermethylated[samples[[3]]$methylation_diff > 0] <- "C"
samples[[4]]$hypermethylated <- "B"
samples[[4]]$hypermethylated[samples[[4]]$methylation_diff > 0] <- "F"
samples[[5]]$hypermethylated <- "C"
samples[[5]]$hypermethylated[samples[[5]]$methylation_diff > 0] <- "D"
samples[[6]]$hypermethylated <- "D"
samples[[6]]$hypermethylated[samples[[6]]$methylation_diff > 0] <- "E"
samples[[7]]$hypermethylated <- "E"
samples[[7]]$hypermethylated[samples[[7]]$methylation_diff > 0] <- "F"

# Make one dataframe from list of diff meth CpGs
diff_meth_sites <- bind_rows(samples)

## -------------------------------------------------------------------------
# How many hypermeth in each sex?
## -------------------------------------------------------------------------

nrow(samples[[1]][samples[[1]]$hypermethylated=="A",]) #3819
nrow(samples[[1]][samples[[1]]$hypermethylated=="B",]) #91

nrow(samples[[2]][samples[[2]]$hypermethylated=="A",]) #163
nrow(samples[[2]][samples[[2]]$hypermethylated=="E",]) #150

nrow(samples[[3]][samples[[3]]$hypermethylated=="B",]) #47
nrow(samples[[3]][samples[[3]]$hypermethylated=="C",]) #945

nrow(samples[[4]][samples[[4]]$hypermethylated=="B",]) #20
nrow(samples[[4]][samples[[4]]$hypermethylated=="F",]) #381

nrow(samples[[5]][samples[[5]]$hypermethylated=="C",]) #51
nrow(samples[[5]][samples[[5]]$hypermethylated=="D",]) #48

nrow(samples[[6]][samples[[6]]$hypermethylated=="D",]) #147
nrow(samples[[6]][samples[[6]]$hypermethylated=="E",]) #654

nrow(samples[[7]][samples[[7]]$hypermethylated=="E",]) #52
nrow(samples[[7]][samples[[7]]$hypermethylated=="F",]) #575

# Goodness of fit
observed = c(575, 52)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions
chisq.test(x = observed,
           p = expected)

# AvsB X-squared = 3554.5, df = 1, p-value < 2.2e-16
# AvsE X-squared = 0.53994, df = 1, p-value = 0.4625
# BvsC X-squared = 812.91, df = 1, p-value < 2.2e-16
# BvsF X-squared = 623.28, df = 1, p-value < 2.2e-16
# CvsD X-squared = 0.090909, df = 1, p-value = 0.763
# DvsE X-squared = 320.91, df = 1, p-value < 2.2e-16
# EvsF X-squared = 436.25, df = 1, p-value < 2.2e-16

## -------------------------------------------------------------------------
# Annotate the differentially methylated CpGs with genomic features
## -------------------------------------------------------------------------

output <- sqldf("SELECT sample.chr,
                    sample.cpg_position,
                    sample.methylation_diff,
                    sample.comparison,
                    sample.hypermethylated,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.number
                    FROM diff_meth_sites AS sample
                    LEFT JOIN annotation AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg_position >= annot.start AND sample.cpg_position <= annot.end)")

output <- output[,-1]
output$feature <- as.factor(output$feature)

## -------------------------------------------------------------------------
# Where are these CpGs?
## -------------------------------------------------------------------------

not_in_feature <- output[is.na(output$gene_id),] #514 total
table(not_in_feature$comparison)
#Sites not in a feature
# AvsB = 90
# AvsE = 15
# BvsF = 64
# BvsC = 162
# CvsD = 6
# DvsE = 7
# EvsF = 170

# Note: 9516 annotations from 7443 CpGs as some fall over multiple annotations

plot_data <- output[output$comparison=="AvsB",]
plot_data <- output[output$comparison=="AvsE",]
plot_data <- output[output$comparison=="BvsC",]
plot_data <- output[output$comparison=="BvsF",]
plot_data <- output[output$comparison=="CvsD",]
plot_data <- output[output$comparison=="DvsE",]
plot_data <- output[output$comparison=="EvsF",]

ggplot(plot_data, aes(x=feature, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Genomic Feature")+
  ylab("Number of Significant CpGs")+
  ggtitle("Male brain vs sperm")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))+
  scale_fill_manual("Hypermethylated",
                    breaks = c("E","F"),
                    labels=c("Male Brain","Sperm"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","exon","intron","lnc_RNA","intergenic"),
                   labels = c("Putative Promoter","5' UTR","3' UTR","Exon","Intron","lnc RNA","Intergenic"),
                   limits =c("exon","intron","promoter","intergenic","three_prime_UTR","five_prime_UTR","lnc_RNA"))

## -------------------------------------------------------------------------
# Number diff CpGs by feature
## -------------------------------------------------------------------------

#Add column to account for different exons/introns so not averaged across all exons/introns in a gene
output$id <- paste0(output$gene_id, output$number)

# Number of CpGs per feature
output_count <- summaryBy(cpg_position ~ feature + id + comparison, data=output, FUN=length)

# Remove intergenic as not really informative
output_count <- output_count[!output_count$feature=="intergenic",]

plot_data <- output_count[output_count$comparison=="AvsB",]
plot_data <- output_count[output_count$comparison=="AvsE",]
plot_data <- output_count[output_count$comparison=="BvsC",]
plot_data <- output_count[output_count$comparison=="BvsF",]
plot_data <- output_count[output_count$comparison=="CvsD",]
plot_data <- output_count[output_count$comparison=="DvsE",]
plot_data <- output_count[output_count$comparison=="EvsF",]

ggplot(plot_data, aes(x=feature, y=cpg_position.length))+
  geom_boxplot()+
  xlab("Feature")+
  ylab("Number of Significant CpGs")+
  ggtitle("Male brain vs sperm")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_y_continuous(breaks = scales::pretty_breaks(n=14))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","exon","intron","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Exon","Intron","lnc RNA"),
                   limits =c("intron","exon","three_prime_UTR","lnc_RNA","promoter","five_prime_UTR"))

## -------------------------------------------------------------------------
# Have a look at exons and introns by position
## -------------------------------------------------------------------------

# Number of CpGs per exon/intron number
output_count2 <- summaryBy(cpg_position ~ feature + id + comparison + number, data=output, FUN=length)

# For exons
exons <- output_count2[output_count2$feature=="exon",]
exons <- exons[!is.na(exons$number),]

ggplot(exons, aes(x=number, y=cpg_position.length))+
  geom_boxplot()+
  xlab("Exon Number")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_y_continuous(labels = label_number(accuracy = 1))+
  facet_grid(rows = vars(comparison),scales = "free")


# For introns
introns <- output_count2[output_count2$feature=="intron",]
introns <- introns[!is.na(introns$number),]

ggplot(introns, aes(x=number, y=cpg_position.length))+
  geom_boxplot()+
  xlab("Intron Number")+
  ylab("Number of Significant CpGs")+
  ggtitle("Male Brains vs Sperm")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_y_continuous(labels = label_number(accuracy = 1))+
  facet_grid(rows = vars(comparison),scales = "free")


