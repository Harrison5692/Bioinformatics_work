# Heatmap generation



#---------table prep




raw <- Mus.musculus_mm10_RAW_COUNT_all_expressions

tpm <- Mus.musculus_mm10_TPM_all_expressions

DESeq2_count <- DESeq2_normalized


#-------------------40 genes working with deseq2 counts--------

counts_40_deseq2 <- DESeq2_count[is.element(row.names(DESeq2_count), Drugs_working$genes),]

#----dropping groups not needed

counts_40_deseq2 = counts_40_deseq2[,-c(7:8)]




#--------------------40 genes working with TPM------------------

counts_40_tpm <- tpm[is.element(tpm$GENE_SYMBOL, Drugs_working$genes),]

counts_40_tpm <- counts_40_tpm[,-1]

counts_40_tpm = counts_40_tpm[,order(names(counts_40_tpm))]
counts_40_tpm = counts_40_tpm[,-c(8,11)]
counts_40_dropa485_tpm = counts_40_tpm[,-c(7:9)]




#---------------------just the 40 genes working raw-------------

counts_40 <- raw[is.element(raw$GENE_SYMBOL, Drugs_working$genes),]

counts_40 <- counts_40[,-1]

counts_40 = counts_40[,order(names(counts_40))]
counts_40 = counts_40[,-11]
counts_40_dropa485 = counts_40[,-c(7:9)]


#--------------------- 11 side effect genes as a result of the inhibitor 

side_effect 



tpm_xxx = de_TGFB_PFCBP1[!de_TGFB_PFCBP1$gene_id %in% DEG_DESq2_genialis_final$X,] 

tpm_side_fx <- tpm[is.element(tpm$GENE_SYMBOL, tpm_xxx$gene_id),]

tpm_side_fx <- tpm_side_fx[,-1]
tpm_side_fx <- tpm_side_fx[,order(names(tpm_side_fx))]
tpm_side_fx <- tpm_side_fx[,-11]
tpm_side_fx_drop_reps <- tpm_side_fx[,-c(1:3,7:9)]


#----------------------------all genes changing in DMSO vs TGFb---------------------------

#-----tpm------------

counts_548 <- tpm[is.element(tpm$GENE_SYMBOL, DEG_DESq2_genialis_final$X),]

counts_548 <- counts_548[,-1]


counts_548 = counts_548[,order(names(counts_548))]
counts_548 = counts_548[,-c(8, 11)]
counts_548_dropa485 = counts_548[,-c(4,7:11)]

df_548 = scale(counts_548_dropa485)


#---------DESeq2------------


counts_548_deseq2 <- DESeq2_count[is.element(row.names(DESeq2_count), DEG_DESq2_genialis_final$X),]

counts_548_deseq2_drop <- counts_548_deseq2[,-c(7:10)]








# here we are going to scale the counts
# this will calculate the mean and the standard deviation of the entire vector 
# then, "scale" each element by those values by subtracting the mean and dividing by the standard deviation
# note that scale(x, scale = FALSE) will only subtract mean and not devide by the SD
# The scale function is the same as (x- mean(x) / SD(x)) (this is possibly wrong...

# This is not changing the data, rather changing the scale ( the axis values when plotting)
# Think of grabbing the axis at the two ends and stretching or compressing it. This is scale
# in contrast, log transforming the data changes the values. The impact of log is "stronger" for larger values
# and more minimal for smaller values

# This provides standerdization of the data. The values it creates are known under several different names, one them being Z-scores

#------------raw----------#


df_40 = scale(counts_40)
df_40_drop = scale(counts_40_dropa485)

#-----------TPM-----------#

#working drugs

df_tpm_40 = scale(counts_40_tpm)
df_tpm_40_drop = scale(counts_40_dropa485_tpm)


#side fx

df_tpm_s_fx_11_drop = scale(tpm_side_fx_drop_reps)
df_tpm_s_fx_11_drop_tgfb1 = scale(tpm_side_fx_drop_reps2)

#-------------------------#
#----------DESeq2---------#
#-------------------------#

#-----40 working genes-----------

df_DESeq2_40 = scale(counts_40_deseq2)


#--------548 tgfb effects

counts_548_deseq2_drop_scaled <- scale(counts_548_deseq2_drop)


counts_548_deseq2_scaled <- scale(counts_548_deseq2)



heatmap(df_40, scale = "row", col = cm.colors(256))


######################################################

library("gplots")


heatmap.2(as.matrix(df_40), Colv=FALSE, scale = 'row', trace = 'none', margins = c(9,10), col = bluered(100)) 
          

heatmap.2(as.matrix(df_40_drop), Colv = FALSE, scale = 'row', trace = 'none', margins = c(10,10), cexRow = .3 , col = bluered(100) ) 




#########################################################

library("pheatmap")

#----------------------------------
#----------raw counts---------------
#----------------------------------


pheatmap(counts_40_dropa485, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete" )

pheatmap(counts_40_dropa485, scale = "none", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete" )

pheatmap(counts_40_dropa485, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "correlation", clustering_method = "complete" )

pheatmap(counts_40_dropa485, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "manhattan", clustering_method = "complete" )


#------------Z scores on raw-------------
pheatmap(df_40_drop, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete" , main = "Z-score of raw counts and row clustering")



# here we cluster by column as well

pheatmap(counts_40_dropa485, scale = "row", cluster_cols = TRUE, fontsize = 5, clustering_distance_rows = "correlation", clustering_method = "complete", main = "Z-score normalizing, column and row clustering" )


pheatmap(df_40_drop, scale = "none", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete" )

#---------------------------------#
#-----------TPM-------------------#
#---------------------------------#

#----------Z score on TPM --------

#----------40 gense from drug working

pheatmap(df_548, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete", labels_row = FALSE, main = 'TPM normalized')

pheatmap(counts_40_dropa485_tpm, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete")

pheatmap(df_tpm_40_drop, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete", main = "Z-scores of TPM normalized reads")

#-----------11 genes from drug side fx

pheatmap(df_tpm_s_fx_11_drop, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete", main = "Z-scores of TPM normalized reads")

#------dropping tgfb1

#------------------------------#
#----------DESeq2--------------#
#------------------------------#

pheatmap(df_DESeq2_40, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete", main = "Z-scores of DESeq2 normalized reads", ) #, cutree_rows = 6)



pheatmap(counts_548_deseq2_drop_scaled, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete", main = "Z-scores of DESeq2 normalized reads") #, cutree_rows = 6)








#scale set to 'row' for values being centered and scaled in the row direction
# "row wise scaling"
# removing the mean (centering) and dividing by the standard deviation (scaling)


#Columns not clusterd
#distance measure used in clustering rows -> euclidean distance
#clustering method is complete





#-------------------genes changing from DMSO vs TGFb------------------------------------------

#pheatmap(counts_548, scale = "row", cluster_cols = FALSE, show_rownames = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete")

#pheatmap(df_548, scale = "row", cluster_cols = FALSE, show_rownames = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete", cutree_rows = 6)





pheatmap(counts_548_deseq2_drop_scaled, scale = "row", cluster_cols = FALSE, show_rownames = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete",  main = "Z-score of DESeq2 normalization")



pheatmap(counts_548_deseq2_scaled, scale = "row", cluster_cols = FALSE, show_rownames = FALSE, clustering_distance_rows = "euclidean", clustering_method = "complete", main = "Z-score of DESeq2 normalization")









pheatmap(df_tpm_s_fx_11_drop_tgfb1, scale = "row", cluster_cols = FALSE, fontsize = 5, clustering_distance_rows = "euclidean", clustering_method = "complete", main = "Z-scores of TPM normalized reads")





#scale set to 'row' for values being centered and scaled in the row direction
#Columns not clusterd
#





library("ggplot2")



expr <- readRDS(paste0(system.file(package = "ComplexHeatmap"),
                       "/extdata/gene_expression.rds"))
mat <- as.matrix(expr[, grep("cell", colnames(expr))])
type <- gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(df = data.frame(type = type))

Heatmap(mat, name = "expression", km = 5, top_annotation = ha, 
        top_annotation_height = unit(4, "mm"), 
        show_row_names = FALSE, show_column_names = FALSE) +
  Heatmap(expr$length, name = "length", width = unit(5, "mm"),
          col = circlize::colorRamp2(c(0, 100000), c("white", "orange"))) +
  Heatmap(expr$type, name = "type", width = unit(5, "mm")) +
  Heatmap(expr$chr, name = "chr", width = unit(5, "mm"),
          col = circlize::rand_color(length(unique(expr$chr))))

