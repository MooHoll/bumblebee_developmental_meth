## -------------------------------------------------------------------------
# Check meth levels of diff methylated genes
## ------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB")
library(readr)
library(data.table)
library(ggplot2)
library(scales)
library(ggpubr)
library(ggthemes)

# Weighted methylation and levels of all genes for each stage
weighted_meth <- read_delim("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_new_annotation/weighted_meth_by_group_bew_annot.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)
weighted_meth<- weighted_meth[weighted_meth$Feature=="gene",]
colnames(weighted_meth)[2] <- "gene_id"

# Read in differentially methylated genes
AvsB_diff_meth_genes_final <- read_csv("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/AvsB_diff_meth_genes_final.txt")
BvsC_diff_meth_genes_final <- read_csv("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/BvsC_diff_meth_genes_final.txt")
CvsD_diff_meth_genes_final <- read_csv("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/CvsD_diff_meth_genes_final.txt")
DvsE_diff_meth_genes_final <- read_csv("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/DvsE_diff_meth_genes_final.txt")
EvsF_diff_meth_genes_final <- read_csv("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/EvsF_diff_meth_genes_final.txt")
AvsE_diff_meth_genes_final <- read_csv("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/AvsE_diff_meth_genes_final.txt")
BvsF_diff_meth_genes_final <- read_csv("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/BvsF_diff_meth_genes_final.txt")

# Pull out the methylation levels
AvsB_all_data <- weighted_meth[weighted_meth$gene_id %in% AvsB_diff_meth_genes_final$gene_id,]
BvsC_all_data <- weighted_meth[weighted_meth$gene_id %in% BvsC_diff_meth_genes_final$gene_id,]
CvsD_all_data <- weighted_meth[weighted_meth$gene_id %in% CvsD_diff_meth_genes_final$gene_id,]
DvsE_all_data <- weighted_meth[weighted_meth$gene_id %in% DvsE_diff_meth_genes_final$gene_id,]
EvsF_all_data <- weighted_meth[weighted_meth$gene_id %in% EvsF_diff_meth_genes_final$gene_id,]
AvsE_all_data <- weighted_meth[weighted_meth$gene_id %in% AvsE_diff_meth_genes_final$gene_id,]
BvsF_all_data <- weighted_meth[weighted_meth$gene_id %in% BvsF_diff_meth_genes_final$gene_id,]

# Make a graph for each
AvsB_all_data <- AvsB_all_data[AvsB_all_data$Stage=="repro_brain" |
                                 AvsB_all_data$Stage=="repro_ovaries",]

a <- ggplot(AvsB_all_data, aes(x=bins, fill=Stage))+
  geom_bar(stat="count", position = "dodge")+
  ggtitle("Worker Brain vs Ovaries")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("high", "medium","low","none"),
                   labels = c("High", "Medium","Low","None"),
                   limits =c("high", "medium","low","none"))

BvsC_all_data <- BvsC_all_data[BvsC_all_data$Stage=="larvae" |
                                 BvsC_all_data$Stage=="repro_ovaries",]

b <- ggplot(BvsC_all_data, aes(x=bins, fill=Stage))+
  geom_bar(stat="count", position = "dodge")+
  ggtitle("Ovaries vs Larvae")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("high", "medium","low","none"),
                   labels = c("High", "Medium","Low","None"),
                   limits =c("high", "medium","low","none"))


CvsD_all_data <- CvsD_all_data[CvsD_all_data$Stage=="larvae" |
                                 CvsD_all_data$Stage=="pupae",]

c <- ggplot(CvsD_all_data, aes(x=bins, fill=Stage))+
  geom_bar(stat="count", position = "dodge")+
  ggtitle("Larvae vs Pupae")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("high", "medium","low","none"),
                   labels = c("High", "Medium","Low","None"),
                   limits =c("high", "medium","low","none"))

DvsE_all_data <- DvsE_all_data[DvsE_all_data$Stage=="pupae" |
                                 DvsE_all_data$Stage=="male_brain",]

d <- ggplot(DvsE_all_data, aes(x=bins, fill=Stage))+
  geom_bar(stat="count", position = "dodge")+
  ggtitle("Pupae vs Male Brain")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("high", "medium","low","none"),
                   labels = c("High", "Medium","Low","None"),
                   limits =c("high", "medium","low","none"))

EvsF_all_data <- EvsF_all_data[EvsF_all_data$Stage=="male_brain" |
                                 EvsF_all_data$Stage=="sperm",]

e <- ggplot(EvsF_all_data, aes(x=bins, fill=Stage))+
  geom_bar(stat="count", position = "dodge")+
  ggtitle("Male Brain vs Sperm")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("high", "medium","low","none"),
                   labels = c("High", "Medium","Low","None"),
                   limits =c("high", "medium","low","none"))

AvsE_all_data <- AvsE_all_data[AvsE_all_data$Stage=="repro_brain" |
                                 AvsE_all_data$Stage=="male_brain",]

f <- ggplot(AvsE_all_data, aes(x=bins, fill=Stage))+
  geom_bar(stat="count", position = "dodge")+
  ggtitle("Worker Brain vs Male Brain")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("high", "medium","low","none"),
                   labels = c("High", "Medium","Low","None"),
                   limits =c("high", "medium","low","none"))

BvsF_all_data <- BvsF_all_data[BvsF_all_data$Stage=="repro_ovaries" |
                                 BvsF_all_data$Stage=="sperm",]

g <- ggplot(BvsF_all_data, aes(x=bins, fill=Stage))+
  geom_bar(stat="count", position = "dodge")+
  ggtitle("Ovaries vs Sperm")+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("high", "medium","low","none"),
                   labels = c("High", "Medium","Low","None"),
                   limits =c("high", "medium","low","none"))


all <- ggarrange(a,b,c,d,e,f,g, ncol=3, nrow=3, common.legend = TRUE, legend="right")

annotate_figure(all, 
                left = text_grob("Number of Genes", 
                                 color = "black", rot = 90, size=20),
                bottom = text_grob("Methylation Level", 
                                   color = "black", size =20,
                                   hjust = 0.75 ))




