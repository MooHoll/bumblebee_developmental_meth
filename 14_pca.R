## -------------------------------------------------------------------------
## Making PCA
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/coverage_reports")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Get a methylkit object for all samples
sample.list <- list("A1.CpG_report.merged_CpG_evidence.cov",
                    "A2.CpG_report.merged_CpG_evidence.cov",
                    "A3.CpG_report.merged_CpG_evidence.cov",
                    "A4.CpG_report.merged_CpG_evidence.cov",
                    "B1.CpG_report.merged_CpG_evidence.cov",
                    "B2.CpG_report.merged_CpG_evidence.cov",
                    "B3.CpG_report.merged_CpG_evidence.cov",
                    "B4.CpG_report.merged_CpG_evidence.cov",
                    "C2.CpG_report.merged_CpG_evidence.cov",
                    "C3.CpG_report.merged_CpG_evidence.cov",
                    "C4.CpG_report.merged_CpG_evidence.cov",
                    "D1.CpG_report.merged_CpG_evidence.cov",
                    "D2.CpG_report.merged_CpG_evidence.cov",
                    "D3.CpG_report.merged_CpG_evidence.cov",
                    "D4.CpG_report.merged_CpG_evidence.cov",
                    "E1.CpG_report.merged_CpG_evidence.cov",
                    "E2.CpG_report.merged_CpG_evidence.cov",
                    "F1.CpG_report.merged_CpG_evidence.cov",
                    "F2.CpG_report.merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list("A1","A2","A3","A4","B1","B2","B3","B4",
                                    "C2","C3","C4","D1","D2","D3","D4",
                                    "E1","E2","F1","F2"),
                   assembly="bter_1.0",
                   treatment=c(0,0,0,0,1,1,1,1,2,2,2,3,3,3,3,4,4,5,5),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)
rm(CPGRaw)
normalized <- normalizeCoverage(filtered_data)
rm(filtered_data)
meth_all_data <- unite(normalized, destrand=F)
#meth_all_data <- unite(normalized, destrand=F, min.per.group = 2L) 
rm(normalized)
nrow(meth_all_data)
# For unite all = 60,368
# For unite 3L = Not enough memory on laptop

df_meth_all <- getData(meth_all_data)

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]
g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]
k <- df_meth_all[,29:30]
l <- df_meth_all[,32:33]
a1 <- df_meth_all[,35:36]
b1 <- df_meth_all[,38:39]
c1 <- df_meth_all[,41:42]
d1 <- df_meth_all[,44:45]
e1 <- df_meth_all[,47:48]
f1 <- df_meth_all[,50:51]
g1 <- df_meth_all[,53:54]
h1 <- df_meth_all[,56:57]
k1 <- df_meth_all[,59:60]


bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h,k,l,a1,b1,c1,d1,e1,f1,g1,h1,k1)) {
  df <- df[complete.cases(df),]
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) 
# For unite all = 1,189
# For unite 3L = 

subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(subset_methBase, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)

PCA_data1$Stage <- c(rep("Worker Brain",4), rep("Worker Ovaries",4),
                     rep("Larvae",3), rep("Pupae",4), 
                     rep("Male Brain",2), rep("Sperm",2))


percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=Stage))+
  geom_point(size=14)+
  #geom_text_repel(aes(label=sample), size=12,show.legend=FALSE, 
  #                point.padding = 2, box.padding = 1)+
  theme_bw()+
  xlab(paste0("PC1:",percentage[1],"variance")) +
  ylab(paste0("PC2:",percentage[2],"variance")) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22),
        legend.text=element_text(size=22),
        legend.title=element_blank())+
  scale_colour_manual(breaks=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                    values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"))


# Add in scree plot to see what's going on with clusters
PCASamples(subset_methBase, screeplot=T, obj.return = T)

# Take a look at other PCs
ggplot(PCA_data1, aes(PC3, PC4, colour=Stage))+
  geom_point(size=14)+
  theme_bw()+
  xlab(paste0("PC7:",percentage[7],"variance")) +
  ylab(paste0("PC8:",percentage[8],"variance")) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank())+
  scale_colour_manual(breaks=c("Worker Brain","Worker Ovaries","Larvae","Pupae", "Male Brain","Sperm"),
                      values=c("#D55E00","#E69F00","#F0E442","#009E73","#0072B2","#CC79A7"))


