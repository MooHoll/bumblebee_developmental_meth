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
library(plyr)
library(tidyverse)

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

#-------------------------------------------------------------------------------------------

# Keep only genes which have a diff meth CpG in exon and a weighted meth difference of 15% for that exon
head(output)
exons_with_a_cpg <- output[output$feature=="exon",]
exons_with_a_cpg <- exons_with_a_cpg[,c(3,4,7,8,9,10)]

weighted_meth <- read_delim("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_new_annotation/weighted_meth_annotation_by_stage_new_annotation.txt", 
                                                               delim = "\t", escape_double = FALSE, 
                                                               trim_ws = TRUE)
weighted_meth_exons <- weighted_meth[weighted_meth$feature=="exon",]
weighted_meth_exons <- weighted_meth_exons[,-1]

both <- merge(exons_with_a_cpg, weighted_meth_exons, by=c("gene_id","start","end"))
table(exons_with_a_cpg$comparison)
table(both$comparison) # they match, this works


# Pull out a dataframe for each comparison now
head(both)
AvsB <- both[both$comparison=="AvsB",]
AvsB <- AvsB[,c(1,5,11,12)]
AvsB$overall_diff <- AvsB$repro_brain - AvsB$repro_ovaries
AvsB$diff_15 <- "no"
AvsB$diff_15[AvsB$overall_diff > 0.15] <- "A"
AvsB$diff_15[AvsB$overall_diff < -0.15] <- "B"

# Remove rows which don't match the same direction of the cpg, only a couple it seems
table(AvsB$diff_15)
counts <- ddply(AvsB, .(AvsB$diff_15, AvsB$hypermethylated), nrow)
AvsB <- AvsB[AvsB$hypermethylated == AvsB$diff_15,]
write.table(AvsB, file="AvsB_filtered_diff_meth_genes.txt", sep="\t", quote=F, col.names = T, row.names = F)

# BvsC
BvsC <- both[both$comparison=="BvsC",]
BvsC <- BvsC[,c(1,5,12,8)]
BvsC$overall_diff <- BvsC$repro_ovaries - BvsC$larvae
BvsC$diff_15 <- "no"
BvsC$diff_15[BvsC$overall_diff > 0.15] <- "B"
BvsC$diff_15[BvsC$overall_diff < -0.15] <- "C"
table(BvsC$diff_15)
ddply(BvsC, .(BvsC$diff_15, BvsC$hypermethylated), nrow)
BvsC <- BvsC[BvsC$hypermethylated == BvsC$diff_15,]
write.table(BvsC, file="BvsC_filtered_diff_meth_genes.txt", sep="\t", quote=F, col.names = T, row.names = F)

# CvsD
CvsD <- both[both$comparison=="CvsD",]
CvsD <- CvsD[,c(1,5,8,10)]
CvsD$overall_diff <- CvsD$larvae - CvsD$pupae
CvsD$diff_15 <- "no"
CvsD$diff_15[CvsD$overall_diff > 0.15] <- "C"
CvsD$diff_15[CvsD$overall_diff < -0.15] <- "D"
table(CvsD$diff_15)
ddply(CvsD, .(CvsD$diff_15, CvsD$hypermethylated), nrow)
CvsD <- CvsD[CvsD$hypermethylated == CvsD$diff_15,]
write.table(CvsD, file="CvsD_filtered_diff_meth_genes.txt", sep="\t", quote=F, col.names = T, row.names = F)

# DvsE
DvsE <- both[both$comparison=="DvsE",]
DvsE <- DvsE[,c(1,5,10,9)]
DvsE$overall_diff <- DvsE$pupae - DvsE$male_brain
DvsE$diff_15 <- "no"
DvsE$diff_15[DvsE$overall_diff > 0.15] <- "D"
DvsE$diff_15[DvsE$overall_diff < -0.15] <- "E"
table(DvsE$diff_15)
ddply(DvsE, .(DvsE$diff_15, DvsE$hypermethylated), nrow)
DvsE <- DvsE[DvsE$hypermethylated == DvsE$diff_15,]
write.table(DvsE, file="DvsE_filtered_diff_meth_genes.txt", sep="\t", quote=F, col.names = T, row.names = F)

# EvsF
EvsF <- both[both$comparison=="EvsF",]
EvsF <- EvsF[,c(1,5,9,13)]
EvsF$overall_diff <- EvsF$male_brain - EvsF$sperm
EvsF$diff_15 <- "no"
EvsF$diff_15[EvsF$overall_diff > 0.15] <- "E"
EvsF$diff_15[EvsF$overall_diff < -0.15] <- "F"
table(EvsF$diff_15)
ddply(EvsF, .(EvsF$diff_15, EvsF$hypermethylated), nrow)
EvsF <- EvsF[EvsF$hypermethylated == EvsF$diff_15,]
write.table(EvsF, file="EvsF_filtered_diff_meth_genes.txt", sep="\t", quote=F, col.names = T, row.names = F)

# AvsE
AvsE <- both[both$comparison=="AvsE",]
AvsE <- AvsE[,c(1,5,11,9)]
AvsE$overall_diff <- AvsE$repro_brain - AvsE$male_brain
AvsE$diff_15 <- "no"
AvsE$diff_15[AvsE$overall_diff > 0.15] <- "A"
AvsE$diff_15[AvsE$overall_diff < -0.15] <- "E"
table(AvsE$diff_15)
ddply(AvsE, .(AvsE$diff_15, AvsE$hypermethylated), nrow)
AvsE <- AvsE[AvsE$hypermethylated == AvsE$diff_15,]
write.table(AvsE, file="AvsE_filtered_diff_meth_genes.txt", sep="\t", quote=F, col.names = T, row.names = F)

# BvsF
BvsF <- both[both$comparison=="BvsF",]
BvsF <- BvsF[,c(1,5,12,13)]
BvsF$overall_diff <- BvsF$repro_ovaries - BvsF$sperm
BvsF$diff_15 <- "no"
BvsF$diff_15[BvsF$overall_diff > 0.15] <- "B"
BvsF$diff_15[BvsF$overall_diff < -0.15] <- "F"
table(BvsF$diff_15)
ddply(BvsF, .(BvsF$diff_15, BvsF$hypermethylated), nrow)
BvsF <- BvsF[BvsF$hypermethylated == BvsF$diff_15,]
write.table(BvsF, file="BvsF_filtered_diff_meth_genes.txt", sep="\t", quote=F, col.names = T, row.names = F)

# Lets make an upset of the overlaps
library(UpSetR)
library(grid)

hypermethylated_AvsB_A <- as.data.frame(AvsB$gene_id[AvsB$hypermethylated=="A"])
colnames(hypermethylated_AvsB_A)<-"gene_id"
hypermethylated_AvsB_B <- as.data.frame(AvsB$gene_id[AvsB$hypermethylated=="B"])
colnames(hypermethylated_AvsB_B)<-"gene_id"

hypermethylated_BvsC_B <- as.data.frame(BvsC$gene_id[BvsC$hypermethylated=="B"])
colnames(hypermethylated_BvsC_B)<-"gene_id"
hypermethylated_BvsC_C <- as.data.frame(BvsC$gene_id[BvsC$hypermethylated=="C"])
colnames(hypermethylated_BvsC_C)<-"gene_id"

hypermethylated_CvsD_C <- as.data.frame(CvsD$gene_id[CvsD$hypermethylated=="C"])
colnames(hypermethylated_CvsD_C)<-"gene_id"
hypermethylated_CvsD_D <- as.data.frame(CvsD$gene_id[CvsD$hypermethylated=="D"])
colnames(hypermethylated_CvsD_D)<-"gene_id"

hypermethylated_DvsE_D <- as.data.frame(DvsE$gene_id[DvsE$hypermethylated=="D"])
colnames(hypermethylated_DvsE_D)<-"gene_id"
hypermethylated_DvsE_E <- as.data.frame(DvsE$gene_id[DvsE$hypermethylated=="E"])
colnames(hypermethylated_DvsE_E)<-"gene_id"

hypermethylated_EvsF_E <- as.data.frame(EvsF$gene_id[EvsF$hypermethylated=="E"])
colnames(hypermethylated_EvsF_E)<-"gene_id"
hypermethylated_EvsF_F <- as.data.frame(EvsF$gene_id[EvsF$hypermethylated=="F"])
colnames(hypermethylated_EvsF_F)<-"gene_id"

hypermethylated_AvsE_A <- as.data.frame(AvsE$gene_id[AvsE$hypermethylated=="A"])
colnames(hypermethylated_AvsE_A)<-"gene_id"
hypermethylated_AvsE_E <- as.data.frame(AvsE$gene_id[AvsE$hypermethylated=="E"])
colnames(hypermethylated_AvsE_E)<-"gene_id"

hypermethylated_BvsF_B <- as.data.frame(BvsF$gene_id[BvsF$hypermethylated=="B"])
colnames(hypermethylated_BvsF_B)<-"gene_id"
hypermethylated_BvsF_F <- as.data.frame(BvsF$gene_id[BvsF$hypermethylated=="F"])
colnames(hypermethylated_BvsF_F)<-"gene_id"

all <- rbind(hypermethylated_AvsB_A,hypermethylated_AvsB_B,hypermethylated_BvsC_B,
             hypermethylated_BvsC_C,hypermethylated_CvsD_C,hypermethylated_CvsD_D,
             hypermethylated_DvsE_D,hypermethylated_DvsE_E,hypermethylated_EvsF_E,
             hypermethylated_EvsF_F,hypermethylated_AvsE_A,hypermethylated_AvsE_E,
             hypermethylated_BvsF_B,hypermethylated_BvsF_F)
all <- as.data.frame(all[!duplicated(all),]) #1516 total genes which change across development
colnames(all) <- "gene_id"

all$`Hypermethylated in worker brain compared to ovaries` <- 0
all$`Hypermethylated in worker brain compared to ovaries`[all$gene_id %in% hypermethylated_AvsB_A$gene_id] <- 1
all$`Hypermethylated in ovaries compared to worker brain` <- 0
all$`Hypermethylated in ovaries compared to worker brain`[all$gene_id %in% hypermethylated_AvsB_B$gene_id] <- 1

all$`Hypermethylated in ovaries compared to larvae` <- 0
all$`Hypermethylated in ovaries compared to larvae`[all$gene_id %in% hypermethylated_BvsC_B$gene_id] <- 1
all$`Hypermethylated in larvae compared to ovaries` <- 0
all$`Hypermethylated in larvae compared to ovaries`[all$gene_id %in% hypermethylated_BvsC_C$gene_id] <- 1

all$`Hypermethylated in larvae compared to pupae` <- 0
all$`Hypermethylated in larvae compared to pupae`[all$gene_id %in% hypermethylated_CvsD_C$gene_id] <- 1
all$`Hypermethylated in pupae compared to larvae` <- 0
all$`Hypermethylated in pupae compared to larvae`[all$gene_id %in% hypermethylated_CvsD_D$gene_id] <- 1

all$`Hypermethylated in pupae compared to male brain` <- 0
all$`Hypermethylated in pupae compared to male brain`[all$gene_id %in% hypermethylated_DvsE_D$gene_id] <- 1
all$`Hypermethylated in male brain compared to pupae` <- 0
all$`Hypermethylated in male brain compared to pupae`[all$gene_id %in% hypermethylated_DvsE_E$gene_id] <- 1

all$`Hypermethylated in male brain compared to sperm` <- 0
all$`Hypermethylated in male brain compared to sperm`[all$gene_id %in% hypermethylated_EvsF_E$gene_id] <- 1
all$`Hypermethylated in sperm compared to male brain` <- 0
all$`Hypermethylated in sperm compared to male brain`[all$gene_id %in% hypermethylated_EvsF_F$gene_id] <- 1

all$`Hypermethylated in worker brain compared to male brain` <- 0
all$`Hypermethylated in worker brain compared to male brain`[all$gene_id %in% hypermethylated_AvsE_A$gene_id] <- 1
all$`Hypermethylated in male brain compared to worker brain` <- 0
all$`Hypermethylated in male brain compared to worker brain`[all$gene_id %in% hypermethylated_AvsE_E$gene_id] <- 1

all$`Hypermethylated in ovaries compared to sperm` <- 0
all$`Hypermethylated in ovaries compared to sperm`[all$gene_id %in% hypermethylated_BvsF_B$gene_id] <- 1
all$`Hypermethylated in sperm compared to ovaries` <- 0
all$`Hypermethylated in sperm compared to ovaries`[all$gene_id %in% hypermethylated_BvsF_F$gene_id] <- 1

upset(all, nsets =14, order.by = "freq",
      text.scale = 1,
      point.size = 3,
      scale.sets = "identity",
      mainbar.y.label =NULL)
grid.text("Intersection Size",x = 0.4, y=0.60, gp=gpar(fontsize=14), rot = 90)


head(all)
write.table(all, file="all_diff_meth_geneIDs_with_category.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

all %>%
  map( function(x) table(x) )

#`Hypermethylated in worker brain compared to ovaries` 1220 
#`Hypermethylated in ovaries compared to worker brain` 2 

#`Hypermethylated in ovaries compared to larvae` 1 
#`Hypermethylated in larvae compared to ovaries` 285 

#`Hypermethylated in larvae compared to pupae` 4 
#`Hypermethylated in pupae compared to larvae` 1 

#`Hypermethylated in pupae compared to male brain` 4 
#`Hypermethylated in male brain compared to pupae` 153 

#`Hypermethylated in male brain compared to sperm` 2 
#`Hypermethylated in sperm compared to male brain` 47 

#`Hypermethylated in worker brain compared to male brain` 12 
#`Hypermethylated in male brain compared to worker brain` 25 

#`Hypermethylated in ovaries compared to sperm` 0
#`Hypermethylated in sperm compared to ovaries` 134 

# Goodness of fit
observed = c(0, 134)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions
chisq.test(x = observed,
           p = expected)
# AvsB X-squared = 1214, df = 1, p-value < 2.2e-16
# BvsC X-squared = 282.01, df = 1, p-value < 2.2e-16
# CvsD X-squared = 1.8, df = 1, p-value = 0.1797
# DvsE X-squared = 141.41, df = 1, p-value < 2.2e-16
# EvsF X-squared = 41.327, df = 1, p-value = 1.288e-10
# AvsE X-squared = 4.5676, df = 1, p-value = 0.03258
# BvsF X-squared = 134, df = 1, p-value < 2.2e-16

