#### box plots #######

#aims:

#global expression or overall abundance

# booxplots of fold change values?

#box plots for gene sets?





########  means ############

myColors_fpkm_means <- ifelse(colnames(cuff_fpkm_means)==c("DMSO"), rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(cuff_fpkm_means)==c("TGFb"), rgb(0.1,0.6,0.7,0.2) ,
                          ifelse(colnames(cuff_fpkm_means)==c("A485"), rgb(0.10,0.6,0.2,0.4),               
                                 "purple" ) ))

#-------------------means FPKM
names(cuff_fpkm_means) <- c("DMSO", "TGFb", "A485", "PFCBP1")

cuff_fpkm_means <- cuffnorm_exprs_fpkm_means

cuff_fpkm_means <- cuff_fpkm_means[,c("DMSO", "TGFb", "A485", "PFCBP1")]


#----everthing---#
boxplot(cuff_fpkm_means, ylim=c(y1=-2,y2=5000),
        col=myColors_fpkm_means, main="fpkm (cuffnorm) Boxplots", xlab= "Sample groups", ylab="Average Total Expression ")
abline(h=median(data.matrix(cuff_fpkm_means)), col='blue')

#----zoom------#
boxplot(cuff_fpkm_means, ylim=c(y1=-4,y2=80),
        col=myColors_fpkm_means, main="fpkm (cuffnorm) Boxplots", xlab= "Sample groups", ylab="Average Total Expression ")
abline(h=median(data.matrix(cuff_fpkm_means)), col='red')

boxplot(cuff_fpkm_means, ylim=c(y1=5,y2=20),
        col=myColors_fpkm_means, main="fpkm (cuffnorm) Boxplots", xlab= "Sample groups", ylab="Average Total Expression ")
abline(h=median(data.matrix(cuff_fpkm_means)), col='red')

#-------log everything-----#
boxplot(log(cuff_fpkm_means), ylim=c(y1=-3,y2=10),
        col=("gray80"), main="fpkm (cuffnorm) Boxplots", xlab= "Sample groups", ylab="Average Total Expression ")
abline(h=median(data.matrix(log(cuff_fpkm_means))), col='blue')


#-------log zoom---------#


boxplot(log(cuff_fpkm_means), ylim=c(y1=2,y2=3),
        col=("gray80"), main="log(fpkm) Boxplots", xlab= "Sample groups", ylab="Average Total Expression ")
abline(h=median(data.matrix(log(cuff_fpkm_means))), col='blue')



#####------------------ All exprs FPKM

myColors2 <- ifelse(colnames(cuffnorm_new)==c("DMSO_1_Group1", "DMSO_2_Group1", "DMSO_3_Group1"), rgb(0.1,0.1,0.7,0.5) , 
                   ifelse(colnames(cuffnorm_new)==c("TGFB_1_Group4", "TGFB_2_Group4", "TGFB_3_Group4"), rgb(0.1,0.6,0.7,0.2) ,
                          ifelse(colnames(cuffnorm_new)==c("TGFB_A485_1_Group2", "TGFB_A485_2_Group2", "TGFb_A485_3_Group2"), rgb(0.10,0.6,0.2,0.4),               
                                 "purple" ) ))

#-----all
boxplot(cuffnorm_new, ylim=c(y1=-1,y2=30),
col=myColors2, main="fpkm Boxplots", xlab= "Sample means", ylab="FPKM expression abundance")
abline(h=median(data.matrix(cuffnorm_new)), col='red')
dev.off()

#----zoom
boxplot(cuffnorm_new, ylim=c(y1=-.5,y2=1.5),
        col=myColors2, main="fpkm Boxplots", xlab= "Sample means", ylab="FPKM expression abundance")
abline(h=median(data.matrix(cuffnorm_new)), col='red')



FPKM_featureCounts <- Mus.musculus_mm10_FPKM_all_expressions[,-1]
FPKM_featureCounts <- FPKM_featureCounts[,order(names(FPKM_featureCounts))]

boxplot(FPKM_featureCounts, ylim=c(y1=-1,y2=60),
        col=("gray80"), main="fpkm (post-featureCounts)  Boxplots", xlab= "Sample means", ylab="overal Gene Expression ")
abline(h=median(data.matrix(FPKM_featureCounts)), col='blue')
               
#------------------ CPM 



#####DROP 8
boxplot(log(cpm_counts), ylim=c(y1=-4,y2=12),
        col=("gray80"), main="cpm Boxplots", las = 1, xlab= "Sample means", ylab="Average overal Expression ")
abline(h=median(cpm_counts), col='blue')


###### DROP 7, 8
boxplot(logcounts, ylim=c(y1=-1,y2=2),
        col=("gray80"), main="log cpm Boxplots", las = 1, xlab= "Sample means", ylab="Average Overal Expression ")
abline(h=median(logcounts), col='blue')


#------------------ DESeq2

#####DROP 7, 8
myColors = levels(DESeq2_count[,1:3]==DESeq2_count[,1:3], rgb(0.1,0.1,0.7,0.5))

myColors <- ifelse(colnames(DESeq2_count)==c("DMSO_1", "DMSO_2", "DMSO_3"), rgb(0.1,0.1,0.7,0.5) , 
            ifelse(colnames(DESeq2_count)==c("TGFB_1", "TGFB_2", "TGFB_3"), rgb(0.1,0.6,0.7,0.2) ,
            ifelse(colnames(DESeq2_count)==c("TGFb_A485_1", "TGFb_A485_3"), rgb(0.10,0.6,0.2,0.4),               
                     "purple" ) ))

#zoom
boxplot(DESeq2_count, ylim=c(y1=2.75,y2=3.25),
        col=myColors, main=" log transformed DESeq2", xlab= "Sample means", ylab="DESeq2 Expression ")
abline(h=median(DESeq2_count), col='red')

#all
boxplot(DESeq2_count, ylim=c(y1=-3,y2=18),
        col=myColors, main=" log transformed DESeq2", xlab= "Sample means", ylab="overal Gene Expression ")
abline(h=median(DESeq2_count), col='red')

###########  log2FC #################





