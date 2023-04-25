# Overlapping gene lists 

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/JSD_index")
library(readr)
library(UpSetR)
library(grid)
library(ggplot2)


jsd <- read_delim("genes_with_JSDindex_new_annotation.txt", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)

diff_meth_genes <- read_delim("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/diff_meth_CpG_lists/all_diff_meth_geneIDs_with_category.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE)

weighted_meth <- read_delim("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_new_annotation/weighted_meth_by_group_bew_annot.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)
# Overlap JSD and diff meth
head(jsd)
head(diff_meth_genes)

jsd$diff_meth <- "no"
jsd$diff_meth[jsd$gene_id %in% diff_meth_genes$gene_id] <- "yes"

ggplot(jsd, aes(x=outlier, fill=diff_meth))+
  geom_bar(stat="count")+
  xlab("Highly Variable Gene")+
  ylab("Gene Count")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))+
  scale_fill_discrete(name = "Differentially\nMethylated")

# Sig diff?
table(jsd[,c(3:4)])
#          diff_meth
#outlier   no  yes
#no       6117 1238
#yes      627  200

6117+1238
1238/7355 #0.17

627+200
200/827 #0.24

# Where are these 200 genes in terms of developmental diffs
diff_variable_gene_list <- jsd$gene_id[jsd$outlier=="yes" & jsd$diff_meth=="yes"]
diff_variable_genes <- data.frame(diff_meth_genes[diff_meth_genes$gene_id %in% diff_variable_gene_list,])
row.names(diff_variable_genes) <- diff_variable_genes$gene_id
diff_variable_genes <- diff_variable_genes[,-1]

upset(diff_variable_genes, nsets =14, order.by = "freq",
      text.scale = 1,
      point.size = 3,
      scale.sets = "identity",
      mainbar.y.label =NULL)
grid.text("Intersection Size",x = 0.5, y=0.60, gp=gpar(fontsize=14), rot = 90)

# Pull out gene lists for GO
highly_variable <- data.frame(jsd$gene_id[jsd$outlier=="yes"])
colnames(highly_variable) <- "gene_id"
write.table(highly_variable, file="highly_variable_as_background.txt", sep="\t",
            quote = F, col.names = T, row.names = F)

head(diff_variable_gene_list)
diff_variable_gene_list <- data.frame(diff_variable_gene_list)
colnames(diff_variable_gene_list) <- "gene_id"
write.table(diff_variable_gene_list, file="highly_variable_and_diff_meth.txt", sep="\t",
            quote = F, col.names = T, row.names = F)
