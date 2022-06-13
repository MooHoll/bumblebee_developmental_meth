## -------------------------------------------------------------------------
# Re-make meth over feature graph for diff levels of methylation
## -------------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_new_annotation")
library(readr)
library(doBy)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Hmisc)
library(scales)
library(ggpubr)
library(ggthemes)
## -------------------------------------------------------------------------
annotation <- read_delim("weighted_meth_annotation_by_stage_new_annotation.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
# removed feature column as running again with all genes, would need to add back in 
#annotation <- annotation[,-c(1,3,4)]
annotation <- annotation[,-c(3,4,5)]

melted_annot <- melt(annotation, id.vars = c("feature","gene_id"))
#colnames(melted_annot) <- c("ID","Stage","Weighted_Methylation")
colnames(melted_annot) <- c("Feature","ID","Stage","Weighted_Methylation")

# Remove rows where NA in one sex
melted_annot <- melted_annot[!is.na(melted_annot$Weighted_Methylation),]

# Remove mRNA and CDS not really informative
melted_annot <- melted_annot[!melted_annot$Feature == "CDS",]
melted_annot <- melted_annot[!melted_annot$Feature == "RNA",]
melted_annot <- melted_annot[!melted_annot$Feature == "mRNA",]
melted_annot <- melted_annot[!melted_annot$Feature == "transcript",]

## -------------------------------------------------------------------------
#### Define summary function (ref:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 
## -------------------------------------------------------------------------
# Normal graph with all information

summary_all<-summarySE(melted_annot, measurevar = "Weighted_Methylation", 
                       #groupvars = "Stage")
                        groupvars = c("Feature","Stage"))
summary_all$Stage <- factor(summary_all$Stage, 
                            levels = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))

ggplot(summary_all, aes(x=Feature, y=Weighted_Methylation, fill=Stage))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Methylation-ci, ymax=Weighted_Methylation+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Weighted Methylation Level")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

# Stats
library(multcomp)
head(melted_annot)
model1<-lm(Weighted_Methylation ~ Stage * Feature , data=melted_annot)
model2<-lm(Weighted_Methylation ~ Stage + Feature , data=melted_annot)
anova(model1,model2) # Significant interaction
summary.lm(model1) # Almost everything sig

## -------------------------------------------------------------------------
# Frequency of features being high/low methylated
head(melted_annot)
meth_low <- melted_annot[melted_annot$Weighted_Methylation < 0.3,]
meth_medium <- melted_annot[melted_annot$Weighted_Methylation > 0.3 &
                              melted_annot$Weighted_Methylation < 0.7,]
meth_high <- melted_annot[melted_annot$Weighted_Methylation > 0.7,]
meth_none <- melted_annot[melted_annot$Weighted_Methylation ==0,]

melted_annot$bins<-"low"
melted_annot$bins[melted_annot$Weighted_Methylation > 0.3 &
                    melted_annot$Weighted_Methylation < 0.7] <-"medium"
melted_annot$bins[melted_annot$Weighted_Methylation > 0.7] <-"high"
melted_annot$bins[melted_annot$Weighted_Methylation ==0] <-"none"
head(melted_annot)
write.table(melted_annot, file="weighted_meth_by_group_bew_annot.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")

melted_annot$combined <- paste0(melted_annot$Feature, "_", melted_annot$Stage)
melted_meth_stuff_2 <- melted_annot
melted_meth_stuff_2$counts <- with(melted_meth_stuff_2, 
                                  ave(bins, Feature, Stage, bins, FUN=length))
plot_data <- melted_meth_stuff_2[!duplicated(melted_meth_stuff_2),]
plot_data$Stage <- as.factor(plot_data$Stage)
plot_data$Stage <- factor(plot_data$Stage, 
                            levels = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))


plot_data_prom <- subset(plot_data, bins =="low")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b1<- ggplot(plot_data_prom, aes(x=Feature, fill=Stage, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("Genomic Feature")+
#  ylab("Count")+
  ggtitle("Low Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

plot_data_prom <- subset(plot_data, bins =="medium")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b2<- ggplot(plot_data_prom, aes(x=Feature, fill=Stage, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  #xlab("Genomic Feature")+
  #ylab("Count")+
  ggtitle("Medium Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

plot_data_prom <- subset(plot_data, bins =="high")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b3<- ggplot(plot_data_prom, aes(x=Feature, fill=Stage, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("")+
  ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title.y=element_text(size=12),
        axis.title.x = element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))


## -------------------------------------------------------------------------
plot_data_prom <- subset(plot_data, bins =="none")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b4 <- ggplot(plot_data_prom, aes(x=Feature, fill=Stage, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("Genomic Feature")+
  #  ylab("Count")+
  ggtitle("No Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))

plot_data_prom <- subset(plot_data, bins =="high")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b3_1<- ggplot(plot_data_prom, aes(x=Feature, fill=Stage, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("")+
 # ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"),
                    labels=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"),
                    limits =c("repro_brain","repro_ovaries","larvae","pupae", "male_brain","sperm"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","intergenic","lnc_RNA"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","Intergenic","lnc RNA"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","lnc_RNA","intergenic"))


## -------------------------------------------------------------------------
# Slightly diff figures
levels_across_feature <- ggarrange(b3_1,b2,b1,b4,
                                   ncol=2, nrow=2, common.legend = TRUE, legend="right")

annotate_figure(levels_across_feature, 
                 left = text_grob("Count", 
                                 color = "black", rot = 90, size=14),
                bottom = text_grob("Genomic Feature", 
                                   color = "black", size =12,
                                   hjust = 0.75 ))
