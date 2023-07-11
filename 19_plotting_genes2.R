## -------------------------------------------------------------------------
## Making diff meth region views
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
# will put them all together in one dataframe at the end

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

for (i in seq_along(datalist)) {  
  colnames(datalist[[i]]) <- c("chr","position","coverage", "Ccount")
  datalist[[i]] <- datalist[[i]][!is.na(datalist[[i]]$coverage),]
  datalist[[i]]$prop <- datalist[[i]]$Ccount / datalist[[i]]$coverage
  datalist[[i]] <- datalist[[i]][,c(1,2,5)]
  colnames(datalist[[i]]) <- c("chr","position",paste0(names(datalist[i])))
}

#all_A = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)
#all_B = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)
#all_C = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)
#all_D = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)
#all_E = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)
all_F = Reduce(function(...) merge(..., all=T, by = c("chr","position")), datalist)

# ------------------------------------

# Remove NA rows
#all_A1 <- all_A[complete.cases(all_A),] # 495395
#all_B1 <- all_B[complete.cases(all_B),] # 3915926
#all_C1 <- all_C[complete.cases(all_C),] # 538564
#all_D1 <- all_D[complete.cases(all_D),] # 4018926
#all_E1 <- all_E[complete.cases(all_E),] # 5898813
#all_F1 <- all_F[complete.cases(all_F),] # 509739

# Incase of any mess ups
all_A1 <- all_A
all_B1 <- all_B
all_C1 <- all_C
all_D1 <- all_D
all_E1 <- all_E
all_F1 <- all_F

# Make an average per replicate set for plotting
all_A1$repro_brain <- rowMeans(all_A1[,c(3:6)], na.rm=T)
all_A1 <- all_A1[,c(1,2,7)]
all_B1$ovaries <- rowMeans(all_B1[,c(3:6)], na.rm=T)
all_B1 <- all_B1[,c(1,2,7)]
all_C1$larvae <- rowMeans(all_C1[,c(3:5)], na.rm=T)
all_C1 <- all_C1[,c(1,2,6)]
all_D1$pupae <- rowMeans(all_D1[,c(3:6)], na.rm=T)
all_D1 <- all_D1[,c(1,2,7)]
all_E1$male_brain <- rowMeans(all_E1[,c(3:4)], na.rm=T)
all_E1 <- all_E1[,c(1,2,5)]
all_F1$sperm <- rowMeans(all_F1[,c(3:4)], na.rm=T)
all_F1 <- all_F1[,c(1,2,5)]

# ------------------------------------

# Put it all together
all_final <- list(all_A1, all_B1, all_C1, all_D1, all_E1, all_F1)
all = Reduce(function(...) merge(..., all=T, by = c("chr","position")), all_final)
nrow(all) #12220545

# Remove rows where everyone is 0
all1 <- all[rowSums(all[,c(3:8)], na.rm = T)>0,] #5239031

# Write out file
write.table(all, file="meth_prop_all_samples.txt", quote=F, col.names=T, row.names=F, sep="\t")

# ------------------------------------
# Now need to change the file to fit the requirements of: http://maplab.imppc.org/methylation_plotter/

# List of regions to look at for each comparison:
# A vs B
# LOC105667080 664 - 771 NW_003570706.1 - hyper B
# LOC100649773 2835574 - 2835966 NC_015769.1 - hyper A 

subset1 <- all1[,c(1:4)]
subset <- subset1[subset1$chr == "NC_015772.1" & (subset1$position > 6860174 & subset1$position < 6860466 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/AvsB_B_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)
# This one
subset <- subset1[subset1$chr == "NC_015762.1" & (subset1$position > 2025446 & subset1$position < 2025783 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/AvsB_A_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)


# B vs C
# LOC100644833 1273152 - 1273257 NW_003566384.1 - hyper C
# LOC100647694 11509318 - 11509371 NC_015767.1 - hyper B

subset1 <- all1[,c(1,2,4,5)]
subset <- subset1[subset1$chr == "NC_015762.1" & (subset1$position > 494091 & subset1$position < 494691 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/BvsC_C_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# This one
subset <- subset1[subset1$chr == "NC_015767.1" & (subset1$position > 2917797 & subset1$position < 2918106 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/BvsC_C1_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# C vs D
# LOC100649539 14632699 - 14632809 	NC_015770.1 - hyper D
# LOC110120112 9609 - 11422 NW_003565741.1 - hyper C

subset1 <- all1[,c(1,2,5,6)]
subset <- subset1[subset1$chr == "NC_015770.1" & (subset1$position > 14632699 & subset1$position < 14632809 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/CvsD_D_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# Promising one
subset <- subset1[subset1$chr == "NW_003565741.1" & (subset1$position > 9609 & subset1$position < 11422 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/CvsD_C_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# D vs E
# LOC100651568 16123978 - 16124682 NC_015768.1 - hyper E
# LOC100650012 7875972 - 7876127 NC_015775.1 - hyper D

# Promising one
subset1 <- all1[,c(1,2,6,7)]
subset <- subset1[subset1$chr == "NC_015768.1" & (subset1$position > 16123978 & subset1$position < 16124682 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/DvsE_E_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

subset <- subset1[subset1$chr == "NC_015775.1" & (subset1$position > 7875972 & subset1$position < 7876127 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/DvsE_D_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# E vs F
# LOC105665864 573339 - 573697 NC_015768.1 - hyper F
# LOC100650138 1604628 - 1604814 NW_003566384.1 - hyper E

# this one
subset1 <- all1[,c(1,2,7,8)]
subset <- subset1[subset1$chr == "NW_003569910.1" & (subset1$position > 1450 & subset1$position < 1644 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/EvsF_F_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

subset <- subset1[subset1$chr == "NC_015779.1" & (subset1$position > 2477815 & subset1$position < 2478037 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/EvsF_F1_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# A vs E
# LOC100646972 5390322 - 5390401 NC_015762.1 - hyper E
# LOC105665864 573339 - 573697 NC_015768.1 - hyper A

subset1 <- all1[,c(1,2,3,7)]
subset <- subset1[subset1$chr == "NC_015762.1" & (subset1$position > 5390322 & subset1$position < 5390401 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/AvsE_E_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

#This one
subset <- subset1[subset1$chr == "NC_015769.1" & (subset1$position > 7958907 & subset1$position < 7959066 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/AvsE_E2_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# B vs F 
# LOC105665864 573339 - 573697 NC_015768.1 - hyper F
# LOC100646594 2008374 - 2008683 NC_015772.1 - hyper B

# This one
subset1 <- all1[,c(1,2,4,8)]
subset <- subset1[subset1$chr == "NC_015771.1" & (subset1$position > 2753188 & subset1$position < 2753621 ),] 
subset <- subset[,-1]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/BvsF_F_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

subset <- subset1[subset1$chr == "NC_015779.1" & (subset1$position > 363458 & subset1$position < 363825 ),] 
subset <- subset[,-1]
subset <- subset[-c(1:4),]
df_transpose = as.data.frame(t(subset))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="../files_for_online_plotter/BvsF_F1_hyper.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)


