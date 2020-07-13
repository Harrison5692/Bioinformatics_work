# DE analysis for Harry

library("DESeq2")
library("ggplot2")
library("sva")
library('writexl')
library('pheatmap')
library('VennDiagram')
library('RColorBrewer')
library("biomaRt")
library("dplyr")


#=============================================
#========Preparing tables for Analysis========
#=============================================

#===========read in the counts ===============

normalized_counts = read.csv(file = 'cntsRUVNormalized-04-23-19.csv')


#========== read in analysis info ============

analysis_info = readxl::read_xlsx('Meta_info.xlsx')


#==========make a loop adding x to every sample in the analysis info =====

analysis_info$PID <- paste("X", analysis_info$PID, sep="")

#========prepare meta_data file ============

meta_data = analysis_info

colnames(meta_data) = c('PID', 'Phenotype', 'Gender')

meta_data[,3] = chartr('12', 'FM', meta_data$Gender)


#==========get gene ID from ensemble ID======

ensb_id = normalized_counts$X

write.table(ensb_id, file = "ensble_id.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


ens_geneID = read.table(file ='Human_ensemble_to_Gsymbol_genes_15k.txt', header = TRUE, fill = TRUE)

GO_Mes = read.table(file ='GO_Mesenchymal_cell_prolif.txt')
Nakamura_up = read.table(file = 'Nakamura_adipogenesis_early_up.txt')
Nakamura_dn = read.table(file = 'Nakamura_adipogenesis_early_dn.txt')
KEGG_apopt = read.table(file ='KEGG_apoptosis_genes.txt')
KEGG_cellcycle = read.table(file ='KEGG_cell_cycle.txt')
KEGG_PPAR = read.table(file ='KEGG_PPAR_signaling.txt')
KEGG_WNT = read.table(file ='KEGG_WNT_Signaling.txt')

ens_geneID2 = ens_geneID

ens_geneID2 = ens_geneID2[!(is.na(ens_geneID2$V2) | ens_geneID2$V2==""), ]

new_counts <- normalized_counts[is.element(normalized_counts$X, ens_geneID2$V1),]

new_counts$Gene_symbol = ens_geneID2$V2

new_counts = new_counts[,-1]

new_counts2 <- new_counts[,is.element(colnames(new_counts), meta_data$PID)]
new_counts2$Gene_symbol = new_counts$Gene_symbol

counts_GO_Mes = new_counts2[is.element(new_counts2$Gene_symbol, GO_Mes$V1),] # lost 8 genes out of 44
counts_Nakamura_up = new_counts2[is.element(new_counts2$Gene_symbol, Nakamura_up$V1),] # lost 15 genes out of 62
counts_Nakamura_dn = new_counts2[is.element(new_counts2$Gene_symbol, Nakamura_dn$V1),] # lost 6 genes out of 37
counts_KEGG_apopt = new_counts2[is.element(new_counts2$Gene_symbol, KEGG_apopt$V1),] # lost 13 genes out of 87
counts_KEGG_cellcycle = new_counts2[is.element(new_counts2$Gene_symbol, KEGG_cellcycle$V1),] # lost 4 genes out of 124
counts_KEGG_PPAR = new_counts2[is.element(new_counts2$Gene_symbol, KEGG_PPAR$V1),] # lost 21 genes out of 69
counts_KEGG_WNT = new_counts2[is.element(new_counts2$Gene_symbol, KEGG_WNT$V1),] # lost 28 genes out of 150

#----------make gene symbols into rownames and remove genesymbol column
row.names(counts_GO_Mes) = counts_GO_Mes$Gene_symbol
counts_GO_Mes = counts_GO_Mes[,-39]

row.names(counts_KEGG_apopt) = counts_KEGG_apopt$Gene_symbol
counts_KEGG_apopt = counts_KEGG_apopt[,-39]

row.names(counts_KEGG_cellcycle) = counts_KEGG_cellcycle$Gene_symbol
counts_KEGG_cellcycle = counts_KEGG_cellcycle[,-39]

row.names(counts_KEGG_PPAR) = counts_KEGG_PPAR$Gene_symbol
counts_KEGG_PPAR = counts_KEGG_PPAR[,-39]

row.names(counts_KEGG_WNT) = counts_KEGG_WNT$Gene_symbol
counts_KEGG_WNT = counts_KEGG_WNT[,-39]

row.names(counts_Nakamura_dn) = counts_Nakamura_dn$Gene_symbol
counts_Nakamura_dn = counts_Nakamura_dn[,-39]

row.names(counts_Nakamura_up) = counts_Nakamura_up$Gene_symbol
counts_Nakamura_up = counts_Nakamura_up[,-39]


#===============================================
#=====begin DE analysis for every geneset=======
#===============================================

intersect(meta_data$PID, colnames(counts_GO_Mes))

#  sample X21630 is missing from counts matrix
#meta_data = meta_data[-15,]


#-------------create deseq object for GO Mesenchymal geneset-------

deseq_obj_GO_Mes <- 
  DESeqDataSetFromMatrix(
    countData = counts_GO_Mes, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Gender # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 


#-------------create deseq object for Nakamura UP-------

deseq_obj_Naka_up <- 
  DESeqDataSetFromMatrix(
    countData = counts_Nakamura_up, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Gender # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 

#-------------create deseq object for Nakamura UP-------

deseq_obj_Naka_dn <- 
  DESeqDataSetFromMatrix(
    countData = counts_Nakamura_dn, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Gender # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 

#-------------create deseq object for KEGG apoptosis-------

deseq_obj_KEGG_apopt <- 
  DESeqDataSetFromMatrix(
    countData = counts_KEGG_apopt, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Gender # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 

#-------------create deseq object for KEGG cell cycle-------

deseq_obj_KEGG_cellcycle <- 
  DESeqDataSetFromMatrix(
    countData = counts_KEGG_cellcycle, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Gender # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 

#-------------create deseq object for KEGG PPAR-------

deseq_obj_KEGG_PPAR <- 
  DESeqDataSetFromMatrix(
    countData = counts_KEGG_PPAR, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Gender # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 

#-------------create deseq object for KEGG WNT-------

deseq_obj_KEGG_WNT <- 
  DESeqDataSetFromMatrix(
    countData = counts_KEGG_WNT, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Gender # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 


#------------now running deseq on all objects-----------

deseq_obj_GO_Mes = DESeq(deseq_obj_GO_Mes)

deseq_obj_Naka_up = DESeq(deseq_obj_Naka_up)

deseq_obj_Naka_dn = DESeq(deseq_obj_Naka_dn)

deseq_obj_KEGG_apopt = DESeq(deseq_obj_KEGG_apopt)

deseq_obj_KEGG_cellcycle = DESeq(deseq_obj_KEGG_cellcycle)

deseq_obj_KEGG_PPAR = DESeq(deseq_obj_KEGG_PPAR)

deseq_obj_KEGG_WNT = DESeq(deseq_obj_KEGG_WNT)

#===========for PCA=============

# normalize_rlog

rlog_counts_GO_Mes <- 
  deseq_obj_GO_Mes %>%
  rlogTransformation()

rlog_counts_KEGG_apopt <- 
  deseq_obj_KEGG_apopt %>%
  rlogTransformation()

rlog_counts_KEGG_cellcycle <- 
  deseq_obj_KEGG_cellcycle %>%
  rlogTransformation()

rlog_counts_KEGG_PPAR <- 
  deseq_obj_KEGG_PPAR %>%
  rlogTransformation()

rlog_counts_KEGG_WNT <- 
  deseq_obj_KEGG_WNT %>%
  rlogTransformation()

rlog_counts_Naka_up <- 
  deseq_obj_Naka_up %>%
  rlogTransformation()

rlog_counts_Naka_dn <- 
  deseq_obj_Naka_dn %>%
  rlogTransformation()

rlogt_counts_GO_Mes = assay(rlog_counts_GO_Mes)
rlogt_counts_Naka_up = assay(rlog_counts_Naka_up)
rlogt_counts_Naka_dn = assay(rlog_counts_Naka_dn)
rlogt_counts_KEGG_apopt = assay(rlog_counts_KEGG_apopt)
rlogt_counts_KEGG_cellcycle = assay(rlog_counts_KEGG_cellcycle)
rlogt_counts_KEGG_PPAR = assay(rlog_counts_KEGG_PPAR)
rlogt_counts_KEGG_WNT = assay(rlog_counts_KEGG_WNT)



#======================================================
#=================== for PCA plots=====================
#======================================================

x = 1

CV_table = as.data.frame(x)

for (x in 1:nrow(rlogt_counts_KEGG_WNT)){
  
  
  CV_table[x] = sd(rlogt_counts_KEGG_WNT[x,1:38])/mean(as.numeric(rlogt_counts_KEGG_WNT[x,1:38]))
  
}
counts_plus_CV = as.data.frame(rlogt_counts_KEGG_WNT)

counts_plus_CV$CV = as.numeric(CV_table)


# get top 1000 rlog CVs (not 1000 here become of small gene subset)
top_1000 = top_n(counts_plus_CV, 1000)

top_1000_genes = top_1000[,-39]

pca_soln <- 
  top_1000_genes %>% # gets counts
  t %>% # input into PCA should be gene x subject
  prcomp()


# add this to your metadata file for ggplot
meta_data$PC1 <- pca_soln$x[ ,1 ]
meta_data$PC2 <- pca_soln$x[ ,2 ]

#----Nake up
meta_data$PC3 <- pca_soln$x[ ,1 ]
meta_data$PC4 <- pca_soln$x[ ,2 ]

#----Nake dn
meta_data$PC5 <- pca_soln$x[ ,1 ]
meta_data$PC6 <- pca_soln$x[ ,2 ]

#----KEGG apopt
meta_data$PC7 <- pca_soln$x[ ,1 ]
meta_data$PC8 <- pca_soln$x[ ,2 ]

#----KEGG cellcycle
meta_data$PC9 <- pca_soln$x[ ,1 ]
meta_data$PC10 <- pca_soln$x[ ,2 ]

#----KEGG PPAR
meta_data$PC11 <- pca_soln$x[ ,1 ]
meta_data$PC12 <- pca_soln$x[ ,2 ]

#----KEGG WNT
meta_data$PC13 <- pca_soln$x[ ,1 ]
meta_data$PC14 <- pca_soln$x[ ,2 ]


#====================================
#=======Plotting PCAs================
#====================================

ggplot(data = meta_data, aes(x = PC1, y = PC2, color = Phenotype, shape = Gender)) +
  geom_point(size = 2) +
  ggtitle('GO_Mesenchymal') +
  scale_alpha_continuous('Phenotype')


ggplot(data = meta_data, aes(x = PC3, y = PC4, color = Phenotype, shape = Gender)) +
  geom_point(size = 2) +
  ggtitle('Nakamura UP') +
  scale_alpha_continuous('Phenotype')

ggplot(data = meta_data, aes(x = PC5, y = PC6, color = Phenotype, shape = Gender)) +
  geom_point(size = 2) +
  ggtitle('Nakamura DN') +
  scale_alpha_continuous('Phenotype')

ggplot(data = meta_data, aes(x = PC7, y = PC8, color = Phenotype, shape = Gender)) +
  geom_point(size = 2) +
  ggtitle('KEGG apoptosis') +
  scale_alpha_continuous('Phenotype')

ggplot(data = meta_data, aes(x = PC9, y = PC10, color = Phenotype, shape = Gender)) +
  geom_point(size = 2) +
  ggtitle('KEGG cell cycle') +
  scale_alpha_continuous('Phenotype')

ggplot(data = meta_data, aes(x = PC11, y = PC12, color = Phenotype, shape = Gender)) +
  geom_point(size = 2) +
  ggtitle('KEGG PPAR') +
  scale_alpha_continuous('Phenotype')

ggplot(data = meta_data, aes(x = PC13, y = PC14, color = Phenotype, shape = Gender)) +
  geom_point(size = 2) +
  ggtitle('KEGG WNT') +
  scale_alpha_continuous('Phenotype')

#=============extracting results ===================================

res_GO_Mes = results(deseq_obj_GO_Mes, contrast = c('Phenotype', 'Ob', 'NW'))
res_Naka_up = results(deseq_obj_Naka_up, contrast = c('Phenotype', 'Ob', 'NW'))
res_Naka_dn = results(deseq_obj_Naka_dn, contrast = c('Phenotype', 'Ob', 'NW'))
res_KEGG_apopt = results(deseq_obj_KEGG_apopt, contrast = c('Phenotype', 'Ob', 'NW'))
res_KEGG_cellcycle = results(deseq_obj_KEGG_cellcycle, contrast = c('Phenotype', 'Ob', 'NW'))
res_KEGG_PPAR = results(deseq_obj_KEGG_PPAR, contrast = c('Phenotype', 'Ob', 'NW'))
res_KEGG_WNT = results(deseq_obj_KEGG_WNT, contrast = c('Phenotype', 'Ob', 'NW'))

rest_GO_Mes = as.data.frame(res_GO_Mes)
rest_Naka_up = as.data.frame(res_Naka_up)
rest_Naka_dn = as.data.frame(res_Naka_dn)
rest_KEGG_apopt = as.data.frame(res_KEGG_apopt)
rest_KEGG_cellcycle = as.data.frame(res_KEGG_cellcycle)
rest_KEGG_PPAR = as.data.frame(res_KEGG_PPAR)
rest_KEGG_WNT = as.data.frame(res_KEGG_WNT)

sub_rest_GO_Mes = subset(rest_GO_Mes, pvalue < .05)
sub_rest_Naka_up = subset(rest_Naka_up, pvalue < .05)
sub_rest_Naka_dn = subset(rest_Naka_dn, pvalue < .05)
sub_rest_KEGG_apopt = subset(rest_KEGG_apopt, pvalue < .05)
sub_rest_KEGG_cellcycle = subset(rest_KEGG_cellcycle, pvalue < .05)
sub_rest_KEGG_PPAR = subset(rest_KEGG_PPAR, pvalue < .05)
sub_rest_KEGG_WNT = subset(rest_KEGG_WNT, pvalue < .05)

#===========write results tables=============

writexl::write_xlsx(
  list(GO_Mes = rest_GO_Mes, Naka_up = rest_Naka_up, Naka_dn = rest_Naka_dn, KEGG_apoptosis = rest_KEGG_apopt,
       KEGG_cellcycle = rest_KEGG_cellcycle, KEGG_PPAR = rest_KEGG_PPAR, KEGG_WNT = rest_KEGG_WNT),
  path = 'Harrys_DE_results_allpaths_20-5-15.xlsx'
  
)

write.csv(sub_rest_GO_Mes, file = 'sub_rest_GO_Mes.csv', quote = FALSE )
write.csv(sub_rest_Naka_up, file = 'sub_rest_Naka_up.csv', quote = FALSE )
write.csv(sub_rest_Naka_dn, file = 'sub_rest_Naka_dn.csv', quote = FALSE )
write.csv(sub_rest_KEGG_apopt, file = 'sub_rest_KEGG_apopt.csv', quote = FALSE )
write.csv(sub_rest_KEGG_cellcycle, file = 'sub_rest_KEGG_cellcycle.csv', quote = FALSE )
write.csv(sub_rest_KEGG_PPAR, file = 'sub_rest_KEGG_PPAR.csv', quote = FALSE )
write.csv(sub_rest_KEGG_WNT, file = 'sub_rest_KEGG_WNT.csv', quote = FALSE )


#=========volcano plots

ggplot(data = rest_GO_Mes, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point( size = 1, alpha = .6, color = 'black') +
  geom_point(data = sub_rest_GO_Mes, size = 1, alpha = 1, color = 'green') +
  labs(title = "Ob vs NW (GO_mesenchymal)", x = "log 2 Fold Change", y = "-log10(P value)", colour="Datasets") +
  geom_vline(xintercept = c(1.5,-1.5)) +
  geom_hline(yintercept = 1.3)


ggplot(data = rest_Naka_up, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point( size = 1, alpha = .6, color = 'black') +
  geom_point(data = sub_rest_Naka_up, size = 1, alpha = 1, color = 'green') +
  labs(title = "Ob vs NW (Nakamura_up)", x = "log 2 Fold Change", y = "-log10(P value)", colour="Datasets") +
  geom_vline(xintercept = c(1.5,-1.5)) +
  geom_hline(yintercept = 1.3)


ggplot(data = rest_Naka_dn, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point( size = 1, alpha = .6, color = 'black') +
  geom_point(data = sub_rest_Naka_dn, size = 1, alpha = 1, color = 'green') +
  labs(title = "Ob vs NW (Nakamura_dn)", x = "log 2 Fold Change", y = "-log10(P value)", colour="Datasets") +
  geom_vline(xintercept = c(1.5,-1.5)) +
  geom_hline(yintercept = 1.3)


ggplot(data = rest_KEGG_apopt, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point( size = 1, alpha = .6, color = 'black') +
  geom_point(data = sub_rest_KEGG_apopt, size = 1, alpha = 1, color = 'green') +
  labs(title = "Ob vs NW (KEGG_apoptosis)", x = "log 2 Fold Change", y = "-log10(P value)", colour="Datasets") +
  geom_vline(xintercept = c(1.5,-1.5)) +
  geom_hline(yintercept = 1.3)


ggplot(data = rest_KEGG_cellcycle, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point( size = 1, alpha = .6, color = 'black') +
  geom_point(data = sub_rest_KEGG_cellcycle, size = 1, alpha = 1, color = 'green') +
  labs(title = "Ob vs NW (KEGG_cell_cycle)", x = "log 2 Fold Change", y = "-log10(P value)", colour="Datasets") +
  geom_vline(xintercept = c(1.5,-1.5)) +
  geom_hline(yintercept = 1.3)


ggplot(data = rest_KEGG_PPAR, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point( size = 1, alpha = .6, color = 'black') +
  geom_point(data = sub_rest_KEGG_PPAR, size = 1, alpha = 1, color = 'green') +
  labs(title = "Ob vs NW (KEGG_PPAR)", x = "log 2 Fold Change", y = "-log10(P value)", colour="Datasets") +
  geom_vline(xintercept = c(1.5,-1.5)) +
  geom_hline(yintercept = 1.3)

ggplot(data = rest_KEGG_WNT, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point( size = 1, alpha = .6, color = 'black') +
  geom_point(data = sub_rest_KEGG_WNT, size = 1, alpha = 1, color = 'green') +
  labs(title = "Ob vs NW (KEGG_WNT)", x = "log 2 Fold Change", y = "-log10(P value)", colour="Datasets") +
  geom_vline(xintercept = c(1.5,-1.5)) +
  geom_hline(yintercept = 1.3)

