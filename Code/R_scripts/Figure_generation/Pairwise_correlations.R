###---------------replicate correlations---------------####

logcounts <- DESeq2_count

logcounts <- scale(logcounts)

#--------DMSO---------

pairs(logcounts[,1:3], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson', line.main = "Pairwise Replicate Correlations (Pearson)")


#-------TGFb----------

pairs(logcounts[,4:6], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson')

#-------TGFb+A485-----

pairs(logcounts[,1:6], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson')

#------TGFb+PFCBO1-----

pairs(logcounts[,10:12], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson')

pairs(logcounts[,c(10,12)], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson')

pairs(logcounts[,10:11], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson')


#-----------TGFB and DMSO---------------------------------


pairs(logcounts[,1:6], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson', line.main = "Pairwise Replicate Correlations (Pearson)")



pairs(logcounts[,4:9], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson', line.main = "Pairwise Replicate Correlations (Pearson)")


pairs(logcounts[,c(1:3,10:12)], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, method = 'pearson', line.main = "Pairwise Replicate Correlations (Pearson)")



#---------------------------------------------------------




overlap7 <- DEG_EdgeR_Genialis[is.element(DEG_EdgeR_Genialis$gene_id, de_TGFB_PFCBP1$gene_id),]




















