
DMSO_TGFB_DESeq2 <- DMSO_vs_TGFb

#-----------P-values

### CHANGING COLOR TO RED FOR LOG2FC > 1 IN MAGNITUDE ###

DMSO_TGFB_DESeq2$color[DMSO_TGFB_DESeq2$logfc > 1 & DMSO_TGFB_DESeq2$pvalue < .05] = 'blue'

DMSO_TGFB_DESeq2$color[DMSO_TGFB_DESeq2$logfc < -1 & DMSO_TGFB_DESeq2$pvalue < .05] = 'red'

DMSO_TGFB_DESeq2[is.na(DMSO_TGFB_DESeq2)] <- 'black'


plot(DMSO_TGFB_DESeq2$logfc, -log10(DMSO_TGFB_DESeq2$pvalue), 
     pch=20, main = 'DMSO vs TGFb Volcano Plot', xlim = c(-3, 3),
     ylab= '-log10 P_Value', xlab= 'Log2 fold change DMSO vs TGFb',
     col = DMSO_TGFB_DESeq2$color)
abline(h=(DMSO_TGFB_DESeq2$pvalue=1.3), v =c(-1,1), col='blue')


#-----------FDR values


DMSO_TGFB_DESeq2$color[DMSO_TGFB_DESeq2$logfc > 1 & DMSO_TGFB_DESeq2$fdr < .05] = 'blue'

DMSO_TGFB_DESeq2$color[DMSO_TGFB_DESeq2$logfc < -1 & DMSO_TGFB_DESeq2$fdr < .05] = 'red'

DMSO_TGFB_DESeq2[is.na(DMSO_TGFB_DESeq2)] <- 'black'


plot(DMSO_TGFB_DESeq2$logfc, -log10(DMSO_TGFB_DESeq2$fdr), 
     pch=20, main = 'DMSO vs TGFb Volcano Plot', xlim = c(-3, 3),
     ylab= '-log10 FDR', xlab= 'Log2 fold change DMSO vs TGFb',
     col = DMSO_TGFB_DESeq2$color)
abline(h=(DMSO_TGFB_DESeq2$fdr=1.3), v =c(-1,1), col='blue')


# Histograms
breaks75 = 75
breaks20 = 20

hist(DMSO_TGFB_DESeq2$pvalue, main = "DMSO vs TGFb P-values", xlab = 'Range', breaks = breaks75)

hist(DMSO_TGFB_DESeq2$fdr, xlim = c(0,1), main = "DMSO vs TGFb P-values", xlab = 'Range', breaks = breaks75)

hist(DMSO_TGFB_DESeq2$pvalue, main = "DMSO vs TGFb P-values", xlab = 'Range', breaks = breaks20)

hist(DMSO_TGFB_DESeq2$fdr, xlim = c(0,1), main = "DMSO vs TGFb FDR", xlab = 'Range', breaks = breaks20)



#=============================
#=====TGFB_vs_A485============
#=============================


TGFB_A485 <- `de_file(1).tab`
TGFB_A485_DESeq2 <- TGFB_A485

### CHANGING COLOR TO RED FOR LOG2FC > 1 IN MAGNITUDE ###

TGFB_A485_DESeq2$color[TGFB_A485_DESeq2$logfc > 1 & TGFB_A485_DESeq2$pvalue < .05] = 'blue'

TGFB_A485_DESeq2$color[TGFB_A485_DESeq2$logfc < -1 & TGFB_A485_DESeq2$pvalue < .05] = 'red'

TGFB_A485_DESeq2[is.na(TGFB_A485_DESeq2)] <- 'black'


plot(TGFB_A485_DESeq2$logfc, -log10(TGFB_A485_DESeq2$pvalue), 
     pch=20, main = 'TGFb vs A485 Volcano Plot', xlim = c(-3, 3),
     ylab= '-log10 P_Value', xlab= 'Log2 fold change DMSO vs TGFb',
     col = TGFB_A485_DESeq2$color)
abline(h=(TGFB_A485_DESeq2$pvalue=1.3), v =c(-1,1), col='blue')


#-----------FDR values


TGFB_A485_DESeq2$color[TGFB_A485_DESeq2$logfc > 1 & TGFB_A485_DESeq2$fdr < .05] = 'blue'

TGFB_A485_DESeq2$color[TGFB_A485_DESeq2$logfc < -1 & TGFB_A485_DESeq2$fdr < .05] = 'red'

TGFB_A485_DESeq2[is.na(TGFB_A485_DESeq2)] <- 'black'


plot(TGFB_A485_DESeq2$logfc, -log10(TGFB_A485_DESeq2$fdr), 
     pch=20, main = 'TGFB vs A485 Volcano Plot', xlim = c(-3, 3), ylim = c(0,2),
     ylab= '-log10 FDR', xlab= 'Log2 fold change TGFB vs A485',
     col = TGFB_A485_DESeq2$color)
abline(h=(TGFB_A485_DESeq2$fdr=1.3), v =c(-1,1), col='blue')


# Histograms
breaks75 = 75
breaks20 = 20

hist(TGFB_A485_DESeq2$pvalue, main = "TGFB vs A485 P-values", xlab = 'Range', breaks = breaks75)

hist(TGFB_A485_DESeq2$fdr, xlim = c(0,1), main = "TGFB vs A485 P-values", xlab = 'Range', breaks = breaks75)

hist(TGFB_A485_DESeq2$pvalue, main = "TGFB vs A485 P-values", xlab = 'Range', breaks = breaks20)

hist(TGFB_A485_DESeq2$fdr, xlim = c(0,1), main = "TGFB vs A485 FDR", xlab = 'Range', breaks = breaks20)

#=============================
#=====TGFB_vs_PFCBP1============
#=============================


TGFB_PFCBP1 <- TGFB_PFCBP1_DESeq2

### CHANGING COLOR TO RED FOR LOG2FC > 1 IN MAGNITUDE ###

TGFB_PFCBP1$color[TGFB_PFCBP1$logfc > 1 & TGFB_PFCBP1$pvalue < .05] = 'blue'

TGFB_PFCBP1$color[TGFB_PFCBP1$logfc < -1 & TGFB_PFCBP1$pvalue < .05] = 'red'

TGFB_PFCBP1[is.na(TGFB_PFCBP1)] <- 'black'


plot(TGFB_PFCBP1$logfc, -log10(TGFB_PFCBP1$pvalue), 
     pch=20, main = 'TGFB vs PFCBP1 Volcano Plot', xlim = c(-3, 3),
     ylab= '-log10 P_Value', xlab= 'Log2 fold change TGFB vs PFCBP1',
     col = TGFB_PFCBP1$color)
abline(h=(TGFB_PFCBP1$pvalue=1.3), v =c(-1,1), col='blue')


#-----------FDR values


TGFB_PFCBP1$color[TGFB_PFCBP1$logfc > 1 & TGFB_PFCBP1$fdr < .05] = 'blue'

TGFB_PFCBP1$color[TGFB_PFCBP1$logfc < -1 & TGFB_PFCBP1$fdr < .05] = 'red'

TGFB_PFCBP1[is.na(TGFB_PFCBP1)] <- 'black'


plot(TGFB_PFCBP1$logfc, -log10(TGFB_PFCBP1$fdr), 
     pch=20, main = 'TGFB vs PFCBP1 Volcano Plot', xlim = c(-3, 3), ylim = c(0,2),
     ylab= '-log10 FDR', xlab= 'Log2 fold change TGFB vs PFCBP1',
     col = TGFB_PFCBP1$color)
abline(h=(TGFB_PFCBP1$fdr=1.3), v =c(-1,1), col='blue')


# Histograms
breaks75 = 75
breaks20 = 20

hist(TGFB_PFCBP1$pvalue, main = "TGFB vs PFCBP1 P-values", xlab = 'Range', breaks = breaks75)

hist(TGFB_PFCBP1$fdr, xlim = c(0,1), main = "TGFB vs PFCBP1 P-values", xlab = 'Range', breaks = breaks75)

hist(TGFB_PFCBP1$pvalue, main = "TGFB vs PFCBP1 P-values", xlab = 'Range', breaks = breaks20)

hist(TGFB_PFCBP1$fdr, xlim = c(0,1), main = "TGFB vs PFCBP1 FDR", xlab = 'Range', breaks = breaks20)

