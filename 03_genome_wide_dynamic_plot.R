#---------------------------------------------------
# Genome wide methylation plots for genes only
#---------------------------------------------------

library(readr)
library(ggplot2)
setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_genes")

# This is a file containing the weighted meth level for each gene for each sample
weighted_meth_all_samples <- read_delim("weighted_meth_by_group_genes_only.txt", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)
weighted_meth_all_samples<-weighted_meth_all_samples[!is.na(weighted_meth_all_samples$bins),]

# Summary function
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
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Mean levels of each category of methylation
summary_overall<-summarySE(weighted_meth_all_samples, measurevar = "Weighted_Methylation", 
                        groupvars = c("bins","Stage"))

ggplot(summary_overall, aes(x=Stage, y=Weighted_Methylation, colour=bins))+
  geom_point(size=3)+
  geom_line(aes(group=bins), size=2)+
  geom_errorbar(aes(ymin=Weighted_Methylation-ci, ymax=Weighted_Methylation+ci),
                width=.2,size=2)+
  scale_x_discrete(breaks=c("repro_brain","repro_ovaries","larvae",
                            "pupae","male_brain","sperm"),
                   limits=c("repro_brain","repro_ovaries","larvae",
                            "pupae","male_brain","sperm"),
                   labels=c("Worker Brain","Worker Ovaries","Larvae Head",
                            "Pupae Head","Male Brain","Sperm"))+
  scale_color_discrete(breaks=c("high","medium","low","none"),
                       limits=c("high","medium","low","none"),
                       labels=c("High","Medium","Low","None"))+
  xlab("Developmental Stage")+
  ylab("Weighted Methylation")+
  theme_bw()+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14,angle = 45,hjust=1),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title = element_blank())

just_meth <- weighted_meth_all_samples[!weighted_meth_all_samples$bins=="none",]
summary_none_category<-summarySE(just_meth, measurevar = "Weighted_Methylation", 
                           groupvars = c("Stage"))

# Mean level of all categories excluding none methylated
ggplot(summary_none_category, aes(x=Stage, y=Weighted_Methylation))+
  geom_point(size=3)+
  geom_line(aes(group=1),size=2)+
  geom_errorbar(aes(ymin=Weighted_Methylation-ci, ymax=Weighted_Methylation+ci),
                width=.2,size=2)+
  scale_x_discrete(breaks=c("repro_brain","repro_ovaries","larvae",
                            "pupae","male_brain","sperm"),
                   limits=c("repro_brain","repro_ovaries","larvae",
                            "pupae","male_brain","sperm"),
                   labels=c("Worker Brain","Worker Ovaries","Larvae Head",
                            "Pupae Head","Male Brain","Sperm"))+
  xlab("Developmental Stage")+
  ylab("Weighted Methylation")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18,angle = 45,hjust=1),
        axis.title=element_text(size=22),
        legend.text=element_text(size=22),
        legend.title = element_blank())

# Now need to do the same but for absolubte counts
# Should therefore take an average across replicates so not losing counts for missing data etc.
head(weighted_meth_all_samples)
library(doBy)
counts_per_group <- summaryBy(Weighted_Methylation ~ ID+Stage, data = weighted_meth_all_samples, FUN=mean) 
length(unique(counts_per_group$ID)) # 10568 genes
counts_per_group$meth_category<- "None"
counts_per_group$meth_category[(counts_per_group$Weighted_Methylation.mean > 0.005 & counts_per_group$Weighted_Methylation.mean <= 0.3)] <- "Low"
counts_per_group$meth_category[(counts_per_group$Weighted_Methylation.mean > 0.3 & counts_per_group$Weighted_Methylation.mean <= 0.7)] <- "Medium"
counts_per_group$meth_category[(counts_per_group$Weighted_Methylation.mean > 0.7 & counts_per_group$Weighted_Methylation.mean <= 1)] <- "High"

counts_per_group_summary <- as.data.frame(table(counts_per_group[,c(4,2)]))
counts_per_group_summary$meth_category = factor(counts_per_group_summary$meth_category, levels=c("High","Medium","Low","None"))

ggplot(counts_per_group_summary, aes(x=Stage, y=Freq, colour=meth_category))+
  geom_point(size=3)+
  geom_line(aes(group=meth_category), size=2)+
  scale_x_discrete(breaks=c("repro_brain","repro_ovaries","larvae",
                            "pupae","male_brain","sperm"),
                   limits=c("repro_brain","repro_ovaries","larvae",
                            "pupae","male_brain","sperm"),
                   labels=c("Worker Brain","Worker Ovaries","Larvae Head",
                            "Pupae Head","Male Brain","Sperm"))+
  xlab("Developmental Stage")+
  ylab("Gene Count")+
  theme_bw()+
  theme(axis.text.y=element_text(size=13),
        axis.text.x=element_text(size=18,angle = 45,hjust=1),
        axis.title=element_text(size=22),
        legend.position = "none",
        strip.text.y = element_text(size = 18))+
  facet_grid(rows = vars(meth_category),scales = "free")+
  scale_colour_manual(breaks=c("High","Medium","Low","None"),
                      values=c("purple4","purple2","mediumpurple2","grey55"))

