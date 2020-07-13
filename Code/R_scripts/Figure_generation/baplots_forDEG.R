#------------------------------------
#--------------Bargraphs for DEG------
#------------------------------------

counts_40_deseq22 <- de_TGFB_PFCBP1[is.element(Drugs_working$genes, de_TGFB_PFCBP1$gene_id),]


up_regulated_deg_working <- subset(counts_40_deseq22, logfc > 0)


down_regulated_def_working <- subset(counts_40_deseq22, logfc < 0)

#------------------------------------------------------------------
  
  
  
bjhbjhj <- counts_40_deseq22[order(counts_40_deseq22$logfc),]  # Test


up_regulated_deg_working <- up_regulated_deg_working[order(up_regulated_deg_working$logfc),]

down_regulated_def_working <- down_regulated_def_working[order(down_regulated_def_working$logfc),]


#------------------------------------------------------------------

barplot(bjhbjhj[,2], names.arg = bjhbjhj$gene_id, horiz = TRUE)
?barplot

barplot2(bjhbjhj[,2], names.arg = bjhbjhj$gene_id, horiz = TRUE, cex.names = .2, las = 1)



barplot2(up_regulated_deg_working[,2], names.arg = up_regulated_deg_working$gene_id, horiz = TRUE, cex.names = .4, las = 1)


barplot2(down_regulated_def_working[,2], names.arg = down_regulated_def_working$gene_id, horiz = TRUE, 
         cex.names = .9, las = 1, col = "red", border = "red4", cex.axis = 1.5,  space = 0)
axis(side=1)



Library(ggplot2)
install.packages("ggplot2")

