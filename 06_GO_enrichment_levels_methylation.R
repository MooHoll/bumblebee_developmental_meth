#-----------------------------------------------
# Script created by Alun Jones, see paper Bebane et al. (2019) Neonics and bumblebees...
#-----------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/GO_analysis")
library(GOstats)
library(GSEABase)
library(treemap)

#-----------------------------------------------
# Read in background GO set and make compatible with GOstats

#GO_annotations <- read.table("./background_levels_GO_lists/all_genes/A_allgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/all_genes/B_allgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/all_genes/C_allgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/all_genes/D_allgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/all_genes/E_allgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/all_genes/F_allgenes_background_gene_set.txt")

#GO_annotations <- read.table("./background_levels_GO_lists/meth_genes/A_methgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/meth_genes/B_methgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/meth_genes/C_methgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/meth_genes/D_methgenes_background_gene_set.txt")
#GO_annotations <- read.table("./background_levels_GO_lists/meth_genes/E_methgenes_background_gene_set.txt")
GO_annotations <- read.table("./background_levels_GO_lists/meth_genes/F_methgenes_background_gene_set.txt")


GO_annotations[,3] <- paste("IEA")
names(GO_annotations) <- c("genes","GOIds","evi")
GO_annotations[,3] <- paste("IEA")
GO_annotations <- GO_annotations[c(2,3,1)]
GO_annotations <- GO_annotations[-1,] # remove the original colnames

#-----------------------------------------------
# Create necessary objects

GO_frame <- GOFrame(GO_annotations,organism = "Bombus terrestris")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))

#-----------------------------------------------
# Read in gene's of interest 

#my_genes <- read.table("./levels_meth_gene_lists/A_no_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/B_no_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/C_no_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/D_no_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/E_no_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/F_no_methylation_genes.txt", header = T)

#my_genes <- read.table("./levels_meth_gene_lists/A_high_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/B_high_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/C_high_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/D_high_methylation_genes.txt", header = T)
#my_genes <- read.table("./levels_meth_gene_lists/E_high_methylation_genes.txt", header = T)
my_genes <- read.table("./levels_meth_gene_lists/F_high_methylation_genes.txt", header = T)

my_genes <- as.data.frame(na.omit(my_genes$geneID))
colnames(my_genes) <- "genes"
my_genes <- as.vector(my_genes[,1])

#-----------------------------------------------
# Keep only genes with annotated GOs

my_genes <- my_genes[my_genes %in% universe]
length(my_genes)
# None:
# A 2781/3939
# B 2892/4045
# C 3010/4227
# D 2766/3927
# E 2909/4118
# F 3176/4384

# Highly meth:
# A 20/30
# B 0/1
# C 36/46
# D 5/9
# E 32/42
# F 401/448

#-----------------------------------------------
# Set up paramters for hypergeometric test

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = my_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      param_list <- c(param_list,parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.05)


#-----------------------------------------------
# Hypergeometric test

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}


GO_enrichment <- Hyper_G_test(param_list = param_list)


#-----------------------------------------------
# Choose what you want to look for, here the choice is biological process over-represented 

Result <- summary(GO_enrichment[["BP_over"]])

#-----------------------------------------------
# FDR correction

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.01,]

#-----------------------------------------------
# Make an output for REVIGO and write out

REVIGO <- Result_FDR[,1:2]

#write.table(REVIGO,"./Levels_FDR0.01/A_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/B_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/C_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/D_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/E_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/F_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./Levels_FDR0.01/A_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/B_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/C_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/D_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.01/E_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"./Levels_FDR0.01/F_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
