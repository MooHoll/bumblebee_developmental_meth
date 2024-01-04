#-----------------------------------------------
# Make a background GO set
#-----------------------------------------------

setwd("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/Overlaps")
library(GOstats)
library(GSEABase)
library(treemap)

highly_variable_as_background <- read_csv("~/Dropbox/Research/Leicester_postdoc/Projects/IDLE/Ben_Developmental_BB/GO_analysis/Revisions/Overlaps/highly_variable_as_background.txt")

Bumble_bee_ensemble_GO_terms <- read_delim("../Bombus_terrestris_HGD_go_annotation.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)
Bumble_bee_ensemble_GO_terms <- Bumble_bee_ensemble_GO_terms[,-2]
colnames(Bumble_bee_ensemble_GO_terms) <- c("gene_id", "GOIds")

highly_variable <- merge(highly_variable_as_background, Bumble_bee_ensemble_GO_terms, by = "gene_id")
length(unique(highly_variable$gene_id)) #777/827
write.table(highly_variable, file="highly_variable_background_gene_set.txt", sep="\t",quote = F, col.names = T, row.names = F)

#-----------------------------------------------
# Script created by Alun Jones, see paper Bebane et al. (2019) Neonics and bumblebees...
#-----------------------------------------------

#-----------------------------------------------
# Read in background GO set and make compatible with GOstats

GO_annotations <- read.table("./highly_variable_background_gene_set.txt")

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

my_genes <- read.table("./highly_variable_and_diff_meth.txt", header = T)

my_genes <- as.data.frame(na.omit(my_genes$gene_id))
colnames(my_genes) <- "genes"
my_genes <- as.vector(my_genes[,1])

#-----------------------------------------------
# Keep only genes with annotated GOs

my_genes <- my_genes[my_genes %in% universe]
length(my_genes)
# 190/200

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

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#-----------------------------------------------
# Make an output for REVIGO and write out

REVIGO <- Result_FDR[,1:2]

write.table(REVIGO,"./highly_variable_and_diff_meth_GOs.txt",row.names = F,sep = "\t",quote = F)

# Get table for supplementary
p_vals <- read_delim("./highly_variable_and_diff_meth_GOs.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

freq_info <- read_delim("./Revigo_BP_Table.tsv", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
freq_info <- freq_info[,c(1,2,5)]
colnames(freq_info)[1] <- "GOBPID"

both <- merge(p_vals, freq_info, by = "GOBPID")

write.table(both, file="info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)


