################################
#### Code by: Ruth De Paula ####
#### 2022 ######################
#### ruthbpaula@gmail.com ######
################################


##### NETWORK ANNOTATION #####

# Creating variables
FMR1_list = read.table("FMR1_genes", sep="\t", header = TRUE)
MBD5_list = read.table("MBD5_genes", sep="\t", header = TRUE)
MECP2_list = read.table("MECP2_genes", sep="\t", header = TRUE)
RAI1_list = read.table("RAI1_genes", sep="\t", header = TRUE)
UBE3A_list = read.table("UBE3A_genes", sep="\t", header = TRUE)
concat_all = read.table("concat_all.txt", sep="\t", header = TRUE)

# Making empty dataframe to store all padj from all tables
pValueMatrix = data.frame(matrix(nrow=length(unique(concat_all$baseMean)),ncol=5))

row.names(pValueMatrix) = unique(concat_all$baseMean)
pValueMatrix = pValueMatrix[-1,]

colnames(pValueMatrix) = c("FMR1","MBD5","MECP2","RAI1","UBE3A")

# Adding padj data into the empty dataframe
for(i in c(1:nrow(pValueMatrix))) {

  pValueMatrix[i,1] = ifelse(row.names(pValueMatrix)[i] %in% FMR1_list$baseMean, 1, 0)
  
  pValueMatrix[i,2] = ifelse(row.names(pValueMatrix)[i] %in% MBD5_list$baseMean, 1, 0)
  
  pValueMatrix[i,3] = ifelse(row.names(pValueMatrix)[i] %in% MECP2_list$baseMean, 1, 0)
  
  pValueMatrix[i,4] = ifelse(row.names(pValueMatrix)[i] %in% RAI1_list$baseMean, 1, 0)
  
  pValueMatrix[i,5] = ifelse(row.names(pValueMatrix)[i] %in% UBE3A_list$baseMean, 1, 0)
}

pValueMatrix$numDisorders = 0
for(i in c(1:nrow(pValueMatrix))){
  pValueMatrix[i, "numDisorders"] = sum(pValueMatrix[i, c(1:5)]==1, na.rm=TRUE)
}

# Use the following file to import into Cytoscape and merge gene name columns
write.table((pValueMatrix), file="pValueMatrix_numDisorders.txt", row.names = TRUE, col.names = TRUE, sep="\t")



### Version with log2Fold instead of zeros and ones

# Creating variables
FMR1_DE_2 = read.table("diff_FMR1_FXS_vs_ctrl_padj_log2Fold_2.txt", sep="\t", header = TRUE)
FMR1_DE_3 = read.table("diff_FMR1_FXS_vs_ctrl_padj_log2Fold_3.txt", sep="\t", header = TRUE)
MECP2_DE_1 = read.table("diff_MECP2_neurons_vs_ctrl_padj_log2Fold.txt", sep="\t", header = TRUE)
MECP2_DE_2 = read.table("diff_MECP2_NPC_vs_ctrl_padj_log2Fold.txt", sep="\t", header = TRUE)
MBD5_DE = read.table("diff_MBD5_MAND_vs_ctrl_padj_log2Fold.txt", sep="\t", header = TRUE)
RAI1_DE_1 = read.table("diff_RAI1_NPC_vs_ctrl_padj_log2Fold_1.txt", sep="\t", header = TRUE)
RAI1_DE_2 = read.table("diff_RAI1_NPC_vs_ctrl_padj_log2Fold_2.txt", sep="\t", header = TRUE)
RAI1_DE_3 = read.table("diff_RAI1_NPC_vs_ctrl_padj_log2Fold_3.txt", sep="\t", header = TRUE)
UBE3A_DE = read.table("diff_UBE3A_neurons_vs_ctrl_padj_log2Fold_1.txt", sep="\t", header = TRUE)
concat_all = read.table("genes/concat_all.txt", sep="\t", header = TRUE)

FMR1_DE_2 = read.table("diff_FMR1_FXS_vs_ctrl_2.txt", sep="\t", header = TRUE)
FMR1_DE_3 = read.table("diff_FMR1_FXS_vs_ctrl_3.txt", sep="\t", header = TRUE)
MECP2_DE_1 = read.table("diff_MECP2_neurons_vs_ctrl.txt", sep="\t", header = TRUE)
MECP2_DE_2 = read.table("diff_MECP2_NPC_vs_ctrl.txt", sep="\t", header = TRUE)
MBD5_DE = read.table("diff_MBD5_MAND_vs_ctrl.txt", sep="\t", header = TRUE)
RAI1_DE_1 = read.table("diff_RAI1_NPC_vs_ctrl_1.txt", sep="\t", header = TRUE)
RAI1_DE_2 = read.table("diff_RAI1_NPC_vs_ctrl_2.txt", sep="\t", header = TRUE)
RAI1_DE_3 = read.table("diff_RAI1_NPC_vs_ctrl_3.txt", sep="\t", header = TRUE)
UBE3A_DE = read.table("diff_UBE3A_neurons_vs_ctrl_1.txt", sep="\t", header = TRUE)
concat_all = read.table("genes/concat_all.txt", sep="\t", header = TRUE)

# Making empty dataframe to store all padj from all tables
pValueMatrix = data.frame(matrix(nrow=length(unique(concat_all$baseMean)),ncol=9))

row.names(pValueMatrix) = unique(concat_all$baseMean)
pValueMatrix = pValueMatrix[-1,]

colnames(pValueMatrix) = c("FMR1_DE_2", "FMR1_DE_3", "MECP2_DE_1", "MECP2_DE_2", "MBD5_DE", "RAI1_DE_1", "RAI1_DE_2", "RAI1_DE_3", "UBE3A_DE")

# Adding padj data into the empty dataframe
for(i in c(1:nrow(pValueMatrix))) {
  
  pValueMatrix[i,1] = ifelse(row.names(pValueMatrix)[i] %in% FMR1_DE_2$genes, FMR1_DE_2[which(FMR1_DE_2$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,2] = ifelse(row.names(pValueMatrix)[i] %in% FMR1_DE_3$genes, FMR1_DE_3[which(FMR1_DE_3$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,3] = ifelse(row.names(pValueMatrix)[i] %in% MECP2_DE_1$genes, MECP2_DE_1[which(MECP2_DE_1$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,4] = ifelse(row.names(pValueMatrix)[i] %in% MECP2_DE_2$genes, MECP2_DE_2[which(MECP2_DE_2$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,5] = ifelse(row.names(pValueMatrix)[i] %in% MBD5_DE$genes, MBD5_DE[which(MBD5_DE$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,6] = ifelse(row.names(pValueMatrix)[i] %in% RAI1_DE_1$genes, RAI1_DE_1[which(RAI1_DE_1$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,7] = ifelse(row.names(pValueMatrix)[i] %in% RAI1_DE_2$genes, RAI1_DE_2[which(RAI1_DE_2$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,8] = ifelse(row.names(pValueMatrix)[i] %in% RAI1_DE_3$genes, RAI1_DE_3[which(RAI1_DE_3$genes==row.names(pValueMatrix)[i]), 3], NA)
  
  pValueMatrix[i,9] = ifelse(row.names(pValueMatrix)[i] %in% UBE3A_DE$genes, UBE3A_DE[which(UBE3A_DE$genes==row.names(pValueMatrix)[i]), 3], NA)
  
}

# Use the following file to import into Cytoscape and merge gene name columns
write.table((pValueMatrix), file="pValueMatrix_log2Fold.txt", row.names = TRUE, col.names = TRUE, sep="\t")


#________________________________________________________________

### USE THE SAME CODE FOR REACTOME PATHWAY INTERSECTION

# Creating variables
FMR1_list = read.table("FMR1_genes_reactome_ids", sep="\t", header = TRUE)
MBD5_list = read.table("MBD5_genes_reactome_ids", sep="\t", header = TRUE)
MECP2_list = read.table("MECP2_genes_reactome_ids", sep="\t", header = TRUE)
RAI1_list = read.table("RAI1_genes_reactome_ids", sep="\t", header = TRUE)
TCF4_list = read.table("TCF4_genes_reactome_ids", sep="\t", header = TRUE)
UBE3A_list = read.table("UBE3A_genes_reactome_ids", sep="\t", header = TRUE)
concat_all = read.table("reactome_concat_all.txt", sep="\t", header = TRUE)

# Making empty dataframe to store all padj from all tables
pValueMatrix = data.frame(matrix(nrow=length(unique(concat_all$baseMean)),ncol=6))

row.names(pValueMatrix) = unique(concat_all$baseMean)
pValueMatrix = pValueMatrix[-1,]

colnames(pValueMatrix) = c("FMR1","MBD5","MECP2","RAI1","TCF4","UBE3A")

# Adding padj data into the empty dataframe
for(i in c(1:nrow(pValueMatrix))) {

  pValueMatrix[i,1] = ifelse(row.names(pValueMatrix)[i] %in% FMR1_list$baseMean, 1, 0)
  
  pValueMatrix[i,2] = ifelse(row.names(pValueMatrix)[i] %in% MBD5_list$baseMean, 1, 0)
  
  pValueMatrix[i,3] = ifelse(row.names(pValueMatrix)[i] %in% MECP2_list$baseMean, 1, 0)
  
  pValueMatrix[i,4] = ifelse(row.names(pValueMatrix)[i] %in% RAI1_list$baseMean, 1, 0)
  
  pValueMatrix[i,5] = ifelse(row.names(pValueMatrix)[i] %in% TCF4_list$baseMean, 1, 0)
  
  pValueMatrix[i,6] = ifelse(row.names(pValueMatrix)[i] %in% UBE3A_list$baseMean, 1, 0)
}

pValueMatrix$numDisorders = 0
for(i in c(1:nrow(pValueMatrix))){
  pValueMatrix[i, "numDisorders"] = sum(pValueMatrix[i, c(1:6)]==1, na.rm=TRUE)
}

# Use the following file to import into Cytoscape and merge gene name columns
write.table((pValueMatrix), file="pValueMatrix_reactome.txt", row.names = TRUE, col.names = TRUE, sep="\t")
