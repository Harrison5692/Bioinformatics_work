#------------------------------------
#--------------Bargraphs for DEG------
#------------------------------------
library(gplots)



# table of chosen genes within deseq results table


#------------most significant in tgfb_vs_pfCBP1


counts_TGFB_vs_PFCBP1 <- de_TGFB_PFCBP1[is.element(de_TGFB_PFCBP1$gene_id, Drugs_working$genes),]

counts_TGFB_vs_PFCBP1_and <- de_TGFB_PFCBP1[is.element(Drugs_working$genes, de_TGFB_PFCBP1$gene_id),]


#----------- 548 from DMSO_vs_TGFb


counts_40_in_dmso_tgfb <- DEG_DESq2_genialis_final[is.element(DEG_DESq2_genialis_final$X, Drugs_working$genes),]

counts40_from_dmso_tgfb <- order(counts_40_in_dmso_tgfb[dgfdgfd$gene_id,])



#--------------------sorting by log2fc-------------------------


up_regulated_deg_working <- subset(counts_TGFB_vs_PFCBP1, logfc > 0)


down_regulated_def_working <- subset(counts_TGFB_vs_PFCBP1, logfc < 0)



up_regulated_deg_working <- subset(counts_TGFB_vs_PFCBP1_and, logfc > 0)


down_regulated_def_working <- subset(counts_TGFB_vs_PFCBP1_and, logfc < 0)




counts_40_in_dmso_tgfb_up <- subset(counts_40_in_dmso_tgfb, dgfdgfd$logfc) 








#------------------------------------------------------------------
  
  
  
bjhbjhj <- counts_40_deseq22[order(counts_40_deseq22$logfc),]  # Test
dgfdgfd <- counts_TGFB_vs_PFCBP1[order(counts_TGFB_vs_PFCBP1$logfc),]






up_regulated_deg_working <- up_regulated_deg_working[order(up_regulated_deg_working$logfc),]

down_regulated_def_working <- down_regulated_def_working[order(down_regulated_def_working$logfc),]


#------------------------------------------------------------------

barplot(bjhbjhj[,2], names.arg = bjhbjhj$gene_id, horiz = TRUE)
?barplot

barplot2(bjhbjhj[,2], names.arg = bjhbjhj$gene_id, horiz = TRUE, cex.names = .2, las = 1)



barplot2(up_regulated_deg_working[,2], names.arg = up_regulated_deg_working$gene_id, horiz = TRUE, cex.names = .6, las = 1, col = "dodgerblue3", space = 0, main = "Up-regulated log2FC in PFCBP1")#, pos = 45)

dev.off()


barplot2(down_regulated_def_working[,2], names.arg = down_regulated_def_working$gene_id, horiz = TRUE, 
         cex.names = .9, las = 1, col = "firebrick2",  space = 0, main = "Down-regulated Log2FC in PFCBP1")


barplot2(counts40_from_dmso_tgfb[,3], names.arg = counts_40_in_dmso_tgfb$X, horiz = TRUE, cex.names = .4, las = 1)



#-------------------fpkm expression table

# cuffnorm
cuffnorm_new <- cuffnorm_all_fpkm_exprs_raw[,order(names(cuffnorm_all_fpkm_exprs_raw))]
#----------------------
# FPKM _feature counts
cuffnorm_new2 <-  Mus.musculus_mm10_FPKM_all_expressions[,-1]

cuffnorm_new2 <- cuffnorm_new2[,order(names(cuffnorm_new2))]

cuffnorm_new2 <- cuffnorm_new2[,-c(2,4,6,8,10,12,14,16,18,20,22,24)]
#-------------------------------------


#bargraph_fpkm <- cuffnorm_new
bargraph_fpkm <- cuffnorm_new2

#bargraph_fpkm = bargraph_fpkm[,order(names(bargraph_fpkm))]


Postn <- c(bargraph_fpkm['Postn', 1], bargraph_fpkm['Postn', 2], bargraph_fpkm['Postn', 3], bargraph_fpkm['Postn', 4], bargraph_fpkm['Postn', 5], bargraph_fpkm['Postn', 6],bargraph_fpkm['Postn', 7],bargraph_fpkm['Postn', 8],bargraph_fpkm['Postn', 9],bargraph_fpkm['Postn', 10], bargraph_fpkm['Postn', 11], bargraph_fpkm['Postn', 12] )

Sertad4 <- c(bargraph_fpkm['Sertad4', 1], bargraph_fpkm['Sertad4', 2], bargraph_fpkm['Sertad4', 3], bargraph_fpkm['Sertad4', 4], bargraph_fpkm['Sertad4', 5], bargraph_fpkm['Sertad4', 6],bargraph_fpkm['Sertad4', 7],bargraph_fpkm['Sertad4', 8],bargraph_fpkm['Sertad4', 9],bargraph_fpkm['Sertad4', 10], bargraph_fpkm['Sertad4', 11],bargraph_fpkm['Sertad4', 12])

Gapdh <- c(bargraph_fpkm['Gapdh', 1], bargraph_fpkm['Gapdh', 2], bargraph_fpkm['Gapdh', 3], bargraph_fpkm['Gapdh', 4], bargraph_fpkm['Gapdh', 5], bargraph_fpkm['Gapdh', 6],bargraph_fpkm['Gapdh', 7],bargraph_fpkm['Gapdh', 8],bargraph_fpkm['Gapdh', 9],bargraph_fpkm['Gapdh', 10])

Smad4 <- c(bargraph_fpkm['Smad4', 1], bargraph_fpkm['Smad4', 2], bargraph_fpkm['Smad4', 3], bargraph_fpkm['Smad4', 4], bargraph_fpkm['Smad4', 5], bargraph_fpkm['Smad4', 6],bargraph_fpkm['Smad4', 7],bargraph_fpkm['Smad4', 8],bargraph_fpkm['Smad4', 9],bargraph_fpkm['Smad4', 10])

Pdgfb <- c(bargraph_fpkm['Pdgfb', 1], bargraph_fpkm['Pdgfb', 2], bargraph_fpkm['Pdgfb', 3], bargraph_fpkm['Pdgfb', 4], bargraph_fpkm['Pdgfb', 5], bargraph_fpkm['Pdgfb', 6],bargraph_fpkm['Pdgfb', 7],bargraph_fpkm['Pdgfb', 8],bargraph_fpkm['Pdgfb', 9],bargraph_fpkm['Pdgfb', 10])

Sox4 <- c(bargraph_fpkm['Sox4', 1], bargraph_fpkm['Sox4', 2], bargraph_fpkm['Sox4', 3], bargraph_fpkm['Sox4', 4], bargraph_fpkm['Sox4', 5], bargraph_fpkm['Sox4', 6],bargraph_fpkm['Sox4', 7],bargraph_fpkm['Sox4', 8],bargraph_fpkm['Sox4', 9],bargraph_fpkm['Sox4', 10])

Brd4 <- c(bargraph_fpkm['Brd4', 1], bargraph_fpkm['Brd4', 2], bargraph_fpkm['Brd4', 3], bargraph_fpkm['Brd4', 4], bargraph_fpkm['Brd4', 5], bargraph_fpkm['Brd4', 6],bargraph_fpkm['Brd4', 7],bargraph_fpkm['Brd4', 8],bargraph_fpkm['Brd4', 9],bargraph_fpkm['Brd4', 10])

Crebbp <- c(bargraph_fpkm['Crebbp', 1], bargraph_fpkm['Crebbp', 2], bargraph_fpkm['Crebbp', 3], bargraph_fpkm['Crebbp', 4], bargraph_fpkm['Crebbp', 5], bargraph_fpkm['Crebbp', 6],bargraph_fpkm['Crebbp', 7],bargraph_fpkm['Crebbp', 8],bargraph_fpkm['Crebbp', 9],bargraph_fpkm['Crebbp', 10])

Acta2 <- c(bargraph_fpkm['Acta2', 1], bargraph_fpkm['Acta2', 2], bargraph_fpkm['Acta2', 3], bargraph_fpkm['Acta2', 4], bargraph_fpkm['Acta2', 5], bargraph_fpkm['Acta2', 6],bargraph_fpkm['Acta2', 7],bargraph_fpkm['Acta2', 8],bargraph_fpkm['Acta2', 9],bargraph_fpkm['Acta2', 10], bargraph_fpkm['Acta2', 11],bargraph_fpkm['Acta2', 12])



myColors2 <- ifelse(colnames(bargraph_fpkm)==c("DMSO_1_Group1", "DMSO_2_Group1", "DMSO_3_Group1"), rgb(0.1,0.1,0.7,0.5) , 
                    ifelse(colnames(bargraph_fpkm)==c("TGFB_1_Group4", "TGFB_2_Group4", "TGFB_3_Group4"), rgb(0.1,0.6,0.7,0.2) ,
                           ifelse(colnames(bargraph_fpkm)==c("TGFB_A485_1_Group2", "TGFB_A485_2_Group2", "TGFb_A485_3_Group2"), rgb(0.10,0.6,0.2,0.4),               
                                  "purple" ) ))

myColors2 <- ifelse(colnames(bargraph_fpkm)==c("DMSO_1", "DMSO_2", "DMSO_3"), rgb(0.1,0.1,0.7,0.5) , 
                    ifelse(colnames(bargraph_fpkm)==c("TGFB_1", "TGFB_2", "TGFB_3"), rgb(0.1,0.6,0.7,0.2) ,
                           ifelse(colnames(bargraph_fpkm)==c("TGFB_A485_1", "TGFB_A485_2", "TGFb_A485_3"), rgb(0.10,0.6,0.2,0.4),               
                                  "purple" ) ))

barplot2(Postn, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Postn", legend.text = TRUE, ylab = "FPKM expression", ylim = c(0, 150))
barplot2(Sertad4, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Sertad4", ylab = "FPKM expression", ylim = c(0,50))
barplot2(Gapdh, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Gapdh", ylab = "FPKM expression", ylim = c(0,15))

barplot2(Smad4, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Smad4", legend.text = TRUE, ylab = "FPKM expression", ylim = c(0, 50))
barplot2(Pdgfb, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Pdgfb", ylab = "FPKM expression", ylim = c(0,10))
barplot2(Sox4, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Sox4", ylab = "FPKM expression", ylim = c(0,40))
barplot2(Brd4, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Brd4", ylab = "FPKM expression", ylim = c(0,15))
barplot2(Crebbp, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Crebbp", ylab = "FPKM expression", ylim = c(0,6))
barplot2(Acta2, names.arg = colnames(bargraph_fpkm), cex.names = .5, col = myColors2, main = "Acta2", ylab = "FPKM expression", ylim = c(0,4000))


#-------------Rlog transformed

DESeq2_count <- DESeq2_normalized_log

myColors <- ifelse(colnames(DESeq2_count)==c("DMSO_1", "DMSO_2", "DMSO_3"), rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(DESeq2_count)==c("TGFB_1", "TGFB_2", "TGFB_3"), rgb(0.1,0.6,0.7,0.2) ,
                          ifelse(colnames(DESeq2_count)==c("TGFb_A485_1", "TGFb_A485_3"), rgb(0.10,0.6,0.2,0.4),               
                                 "purple" ) ))

Sertad4 <- c(DESeq2_count['Sertad4', 1], DESeq2_count['Sertad4', 2], DESeq2_count['Sertad4', 3], DESeq2_count['Sertad4', 4], DESeq2_count['Sertad4', 5], DESeq2_count['Sertad4', 6],DESeq2_count['Sertad4', 7],DESeq2_count['Sertad4', 8],DESeq2_count['Sertad4', 9],DESeq2_count['Sertad4', 10])

Postn <- c(DESeq2_count['Postn', 1], DESeq2_count['Postn', 2], DESeq2_count['Postn', 3], DESeq2_count['Postn', 4], DESeq2_count['Postn', 5], DESeq2_count['Postn', 6],DESeq2_count['Postn', 7],DESeq2_count['Postn', 8],DESeq2_count['Postn', 9],DESeq2_count['Postn', 10])

Gapdh <- c(DESeq2_count['Gapdh', 1], DESeq2_count['Gapdh', 2], DESeq2_count['Gapdh', 3], DESeq2_count['Gapdh', 4], DESeq2_count['Gapdh', 5], DESeq2_count['Gapdh', 6],DESeq2_count['Gapdh', 7],DESeq2_count['Gapdh', 8],DESeq2_count['Gapdh', 9],DESeq2_count['Gapdh', 10])

Smad4 <- c(DESeq2_count['Smad4', 1], DESeq2_count['Smad4', 2], DESeq2_count['Smad4', 3], DESeq2_count['Smad4', 4], DESeq2_count['Smad4', 5], DESeq2_count['Smad4', 6],DESeq2_count['Smad4', 7],DESeq2_count['Smad4', 8],DESeq2_count['Smad4', 9],DESeq2_count['Smad4', 10])

Pdgfb <- c(DESeq2_count['Pdgfb', 1], DESeq2_count['Pdgfb', 2], DESeq2_count['Pdgfb', 3], DESeq2_count['Pdgfb', 4], DESeq2_count['Pdgfb', 5], DESeq2_count['Pdgfb', 6],DESeq2_count['Pdgfb', 7],DESeq2_count['Pdgfb', 8],DESeq2_count['Pdgfb', 9],DESeq2_count['Pdgfb', 10])

Sox4 <- c(DESeq2_count['Sox4', 1], DESeq2_count['Sox4', 2], DESeq2_count['Sox4', 3], DESeq2_count['Sox4', 4], DESeq2_count['Sox4', 5], DESeq2_count['Sox4', 6],DESeq2_count['Sox4', 7],DESeq2_count['Sox4', 8],DESeq2_count['Sox4', 9],DESeq2_count['Sox4', 10])

Brd4 <- c(DESeq2_count['Brd4', 1], DESeq2_count['Brd4', 2], DESeq2_count['Brd4', 3], DESeq2_count['Brd4', 4], DESeq2_count['Brd4', 5], DESeq2_count['Brd4', 6],DESeq2_count['Brd4', 7],DESeq2_count['Brd4', 8],DESeq2_count['Brd4', 9],DESeq2_count['Brd4', 10])

Crebbp <- c(DESeq2_count['Crebbp', 1], DESeq2_count['Crebbp', 2], DESeq2_count['Crebbp', 3], DESeq2_count['Crebbp', 4], DESeq2_count['Crebbp', 5], DESeq2_count['Crebbp', 6],DESeq2_count['Crebbp', 7],DESeq2_count['Crebbp', 8],DESeq2_count['Crebbp', 9],DESeq2_count['Crebbp', 10])



barplot2(Postn, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Postn", ylab = "Log Transformed expression", ylim = c(0,18))
barplot2(Sertad4, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Sertad4", ylab = "Log Transformed expression", ylim = c(0,18))
barplot2(Gapdh, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Gapdh", ylab = "Log Transformed expression", ylim = c(0,18))

barplot2(Smad4, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Smad4", legend.text = TRUE, ylab = "Log Transformed expression", ylim = c(0, 12))
barplot2(Pdgfb, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Pdgfb", ylab = "Log Transformed expression", ylim = c(0,8))
barplot2(Sox4, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Sox4", ylab = "Log Transformed expression", ylim = c(0,10))
barplot2(Brd4, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Brd4", ylab = "Log Transformed expression", ylim = c(0,10))
barplot2(Crebbp, names.arg = colnames(DESeq2_count), cex.names = .5, col = myColors2, main = "Crebbp", ylab = "Log Transformed expression", ylim = c(0,10))


#----------------TMM

tmm <- tmm

myColors <- ifelse(colnames(tmm)==c("DMSO_1", "DMSO_2", "DMSO_3"), rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(tmm)==c("TGFB_1", "TGFB_2", "TGFB_3"), rgb(0.1,0.6,0.7,0.2) ,
                          ifelse(colnames(tmm)==c("TGFb_A485_1", "TGFb_A485_3"), rgb(0.10,0.6,0.2,0.4),               
                                 "purple" ) ))

Sertad4 <- c(tmm['Sertad4', 1], tmm['Sertad4', 2], tmm['Sertad4', 3], tmm['Sertad4', 4], tmm['Sertad4', 5], tmm['Sertad4', 6],tmm['Sertad4', 7],tmm['Sertad4', 8],tmm['Sertad4', 9],tmm['Sertad4', 10])

Postn <- c(tmm['Postn', 1], tmm['Postn', 2], tmm['Postn', 3], tmm['Postn', 4], tmm['Postn', 5], tmm['Postn', 6],tmm['Postn', 7],tmm['Postn', 8],tmm['Postn', 9],tmm['Postn', 10])

Gapdh <- c(tmm['Gapdh', 1], tmm['Gapdh', 2], tmm['Gapdh', 3], tmm['Gapdh', 4], tmm['Gapdh', 5], tmm['Gapdh', 6],tmm['Gapdh', 7],tmm['Gapdh', 8],tmm['Gapdh', 9],tmm['Gapdh', 10])

Smad4 <- c(tmm['Smad4', 1], tmm['Smad4', 2], tmm['Smad4', 3], tmm['Smad4', 4], tmm['Smad4', 5], tmm['Smad4', 6],tmm['Smad4', 7],tmm['Smad4', 8],tmm['Smad4', 9],tmm['Smad4', 10])

Pdgfb <- c(tmm['Pdgfb', 1], tmm['Pdgfb', 2], tmm['Pdgfb', 3], tmm['Pdgfb', 4], tmm['Pdgfb', 5], tmm['Pdgfb', 6],tmm['Pdgfb', 7],tmm['Pdgfb', 8],tmm['Pdgfb', 9],tmm['Pdgfb', 10])

Sox4 <- c(tmm['Sox4', 1], tmm['Sox4', 2], tmm['Sox4', 3], tmm['Sox4', 4], tmm['Sox4', 5], tmm['Sox4', 6],tmm['Sox4', 7],tmm['Sox4', 8],tmm['Sox4', 9],tmm['Sox4', 10])

Brd4 <- c(tmm['Brd4', 1], tmm['Brd4', 2], tmm['Brd4', 3], tmm['Brd4', 4], tmm['Brd4', 5], tmm['Brd4', 6],tmm['Brd4', 7],tmm['Brd4', 8],tmm['Brd4', 9],tmm['Brd4', 10])

Crebbp <- c(tmm['Crebbp', 1], tmm['Crebbp', 2], tmm['Crebbp', 3], tmm['Crebbp', 4], tmm['Crebbp', 5], tmm['Crebbp', 6],tmm['Crebbp', 7],tmm['Crebbp', 8],tmm['Crebbp', 9],tmm['Crebbp', 10])



barplot2(Postn, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Postn", ylab = "TMM expression", ylim = c(0, 1800))
barplot2(Sertad4, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Sertad4", ylab = "TMM expression", ylim = c(0,150))
barplot2(Gapdh, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Gapdh", ylab = "TMM expression", ylim = c(0, 1000))

barplot2(Smad4, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Smad4", legend.text = TRUE, ylab = "TMM expression", ylim = c(0, 150))
barplot2(Pdgfb, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Pdgfb", ylab = "TMM expression", ylim = c(0,25))
barplot2(Sox4, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Sox4", ylab = "TMM Transformed expression", ylim = c(0,100))
barplot2(Brd4, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Brd4", ylab = "TMM Transformed expression", ylim = c(0,80))
barplot2(Crebbp, names.arg = colnames(tmm), cex.names = .5, col = myColors2, main = "Crebbp", ylab = "TMM Transformed expression", ylim = c(0,40))

