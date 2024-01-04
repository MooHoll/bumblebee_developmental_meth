## -------------------------------------------------------------------------
## Binomial test for number of CpGs methylated
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/coverage_reports")
library(grid)
library(readr)
library(ggplot2)
library(ggdendro)
library(ggfortify)
library(ggrepel)
library(methylKit)

## -------------------------------------------------------------------------
# Process the below in chunks for the binomial as the computer doesn't have enough memory

# Get a methylkit object for all samples
sample.list <- list(#"A1.CpG_report.merged_CpG_evidence.cov",
 #"A2.CpG_report.merged_CpG_evidence.cov",
  #"A3.CpG_report.merged_CpG_evidence.cov",
  #"A4.CpG_report.merged_CpG_evidence.cov")
  #"B1.CpG_report.merged_CpG_evidence.cov",
  #"B2.CpG_report.merged_CpG_evidence.cov",
  #"B3.CpG_report.merged_CpG_evidence.cov",
  #"B4.CpG_report.merged_CpG_evidence.cov")
  #"C2.CpG_report.merged_CpG_evidence.cov",
  #"C3.CpG_report.merged_CpG_evidence.cov",
  #"C4.CpG_report.merged_CpG_evidence.cov")
  #"D1.CpG_report.merged_CpG_evidence.cov",
  #"D2.CpG_report.merged_CpG_evidence.cov",
  #"D3.CpG_report.merged_CpG_evidence.cov",
  #"D4.CpG_report.merged_CpG_evidence.cov")
  #"E1.CpG_report.merged_CpG_evidence.cov",
  #"E2.CpG_report.merged_CpG_evidence.cov")
  "F1.CpG_report.merged_CpG_evidence.cov",
  "F2.CpG_report.merged_CpG_evidence.cov")

CPGRaw <- methRead(sample.list, 
                   sample.id = list(#"A1","A2","A3","A4"),
                     #"B1","B2","B3","B4"),
                     #"C2","C3","C4"),
                     #"D1","D2","D3","D4"),
                     #"E1","E2"),
                     "F1","F2"),
                   assembly="bter_1.0",
                   treatment=c(0,0),
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
meth_all_data <- unite(normalized, destrand=F, min.per.group = 1L)
rm(normalized)
nrow(meth_all_data)
# A = 9627086
# B = 11458877
# C = 7938691
# D = 11945788
# E = 11197117
# F = 4136208

df_meth_all <- getData(meth_all_data)

a <- df_meth_all[,c(1:2,5:6)]
b <- df_meth_all[,c(1:2,8:9)]
#c <- df_meth_all[,c(1:2,11:12)]
#d <- df_meth_all[,c(1:2,14:15)]

datalist <- list(a,b)

sample.id = c(#"A1","A2","A3","A4")
  #"B1","B2","B3","B4")
  #"C2","C3","C4")
  #"D1","D2","D3","D4")
  #"E1","E2")
  "F1","F2")
names(datalist) <- sample.id

# Probability of success = average non-conversion rate
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater") $p.value}

for (i in seq_along(datalist)) {  
  colnames(datalist[[i]]) <- c("chr","position","coverage", "Ccount")
  datalist[[i]] <- datalist[[i]][!is.na(datalist[[i]]$coverage),]
  datalist[[i]]$pVal <- mapply(bt, datalist[[i]]$Ccount, datalist[[i]]$coverage)
  datalist[[i]]$FDR <- p.adjust(datalist[[i]]$pVal, method = "BH", n = length(datalist[[i]]$pVal))
  datalist[[i]]$methylated <- 0
  datalist[[i]]$methylated[datalist[[i]]$FDR < 0.05] <- 1
  datalist[[i]] <- datalist[[i]][,c(1,2,7)]
  colnames(datalist[[i]]) <- c("chr","position",paste0(names(datalist[i])))
}


# Number methylated CpGs per sample
table(datalist[[1]][,3])
table(datalist[[2]][,3])
#table(datalist[[3]][,3])
#table(datalist[[4]][,3])

# NOTE: number of CpGs within the genome = 25055682
# A1 0.108642
(27221 / 25055682)*100
# A2 0.07012781
(17571 / 25055682)*100
# A3 0.2032673
(50930 / 25055682)*100
# A4 0.1035693
(25950 / 25055682)*100
# B1 0.1104181
(27666 / 25055682)*100
# B2 0.1510276
(37841 / 25055682)*100
# B3 0.1031463
(25844 / 25055682)*100
# B4 0.1386991
(34752 / 25055682)*100
# C2 0.059
(14931 / 25055682)*100
# C3 0.074
(18542 / 25055682)*100
# C4 0.038
(9694 / 25055682)*100
# D1 0.2921094
(73190 / 25055682)*100
# D2 0.1242393
(31129 / 25055682)*100
# D3 0.2972779
(74485 / 25055682)*100
# D4 0.1862811
(46674 / 25055682)*100
# E1 0.33044
(82794 / 25055682)*100
# E2 0.1243989
(31169 / 25055682)*100
# F1 0.02179146
(5460 / 25055682)*100
# F2 0.017509
(4387 / 25055682)*100
