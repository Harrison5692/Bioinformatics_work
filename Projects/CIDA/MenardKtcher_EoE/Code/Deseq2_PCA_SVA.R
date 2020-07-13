#===============================================================================
#------------P20003 Menard Katcher EoE Project DE code--------------------------
# ==============================================================================
library("dplyr")
library("DESeq2")
library("ggplot2")

#This is for new PCA plots without the control

EoE2_raw_counts <- read.csv(file = "~/Projects/Calies_Menard-Katcher/Data/Counts_matrices/EoE2_raw_counts.txt", sep = '', row.names = 1,header = TRUE)
EoE2_raw_counts <- EoE2_raw_counts[,-8]
EoE2_raw_counts2 <- EoE2_raw_counts[,-c(1,2,3,4,5,6,7,8)]

# Pre-filtering out 2/3 of samples without a single expression count

EoE2_raw_counts <- EoE2_raw_counts[rowSums(EoE2_raw_counts==0)<18,]

meta_data <- read.csv("~/Projects/Calies_Menard-Katcher/data/Meta_data/Meta_data_condensed.csv", header = TRUE, row.names = 1)
meta_data <- meta_data[-8,]
row.names(meta_data) <- toupper(row.names(meta_data))
colnames(EoE2_raw_counts) <- toupper(colnames(EoE2_raw_counts))
all(rownames(meta_data)==colnames(EoE2_raw_counts))

# load counts from Matt, then put into DESeq2 format
deseq_obj <- 
  DESeqDataSetFromMatrix(
    countData = EoE2_raw_counts, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Age..years., # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 


# (need to study this)
plotDispEsts(deseq_obj)

# normalize_vst
vst_counts <-
  deseq_obj %>%
  varianceStabilizingTransformation()

# normalize_rlog

rlog_counts <- 
  deseq_obj %>%
  rlogTransformation()


# get Coefficient of variation of all genes
rlog_counts_t = assay(rlog_counts)

CV_table = as.data.frame(x)

for (x in 1:nrow(rlog_counts_t)){
  
  
  CV_table[x] = sd(rlog_counts_t[x,1:27])/mean(as.numeric(rlog_counts_t[x,1:27]))
  
  }

counts_plus_CV2 = as.data.frame(rlog_counts_t)

counts_plus_CV2$CV = as.numeric(CV_table)



# get top 1000 rlog CVs
top_1000 = top_n(counts_plus_CV2, 1000)

top_1000_genes = top_1000[,-28]

##################
# run PCA#########
##################

pca_soln <- 
  top_1000_genes %>% # gets counts
  t %>% # input into PCA should be gene x subject
  prcomp()


# add this to your metadata file for ggplot
meta_data$PC1 <- pca_soln$x[ ,1 ]
meta_data$PC2 <- pca_soln$x[ ,2 ]

# finally, plot
ggplot(data = meta_data, aes(x = PC1, y = PC2, color = meta_data[,3], shape = Phenotype)) +
  geom_point(size = 2) +
  ggtitle('Age')
  

ggplot(data = meta_data, aes(x = PC1, y = PC2, color = meta_data[,5], shape = Phenotype)) +
  geom_point(size = 2) +
  ggtitle('PPI Drug')

ggplot(data = meta_data, aes(x = PC1, y = PC2, color = meta_data[,6], shape = Phenotype)) +
  geom_point(size = 2) +
  ggtitle('STC Drug')

ggplot(data = meta_data, aes(x = PC1, y = PC2, color = meta_data[,7], shape = Phenotype)) +
  geom_point(size = 2) +
  ggtitle('Peak eos/hpf')

ggplot(data = meta_data, aes(x = PC1, y = PC2,shape = Phenotype, color = Phenotype)) +
  geom_point(size = 2) +
  ggtitle('Phenotype') 
 # geom_text(show.legend = FALSE)
  
ggplot(data = meta_data, aes(x = PC1, y = PC2, color = meta_data[,9], shape = Phenotype)) +
  geom_point(size = 2) +
  ggtitle('RIN')

##########################
## SVA analysis###########
##########################

library("sva")

#remove rows with zero counts
rlog_t2 = rlog_counts_t[rowSums(rlog_counts_t[, -1] > 0) != 0, ]

svobj = sva(rlog_counts_t, # should be transformed count matrix of genes x subject
    mod = model.matrix( ~ Phenotype + Age..years., data = meta_data),
    mod0 = model.matrix( ~ Age..years., data = meta_data))


str(svobj)
svobj$sv
meta_data$SV1 = svobj$sv[,1]
meta_data$SV2 = svobj$sv[,2]
meta_data$SV3 = svobj$sv[,3]
meta_data$SV4 = svobj$sv[,4]


#####################################################
# now to redo DE analyis with the surrogate variables
#####################################################

#check tables for accuracy
all(rownames(meta_data)==colnames(EoE2_raw_counts))

EoE2_raw_counts

# load counts again but this time include Surrogate variables from meta data in the design
deseq_obj2 <- 
  DESeqDataSetFromMatrix(
    countData = EoE2_raw_counts, # DESeq2 wants this in gene (rows) x subject (columns) form
    colData = meta_data,
    design =  ~ Phenotype + Age..years. + SV1 + SV2 + SV3 + SV4, # just keep this ~ Group for now
  ) %>%
  estimateSizeFactors() %>% # these steps important for normalization
  estimateDispersions() 

plotDispEsts(deseq_obj2)

deseq_obj2 = DESeq(deseq_obj2)

resultsNames(deseq_obj2)

rlog_new = rlog(deseq_obj2)

rlog_table = assay(rlog_new)

#Here we view differences between phenotype

#----------getting the results
res1 <- results(deseq_obj2, contrast = c("Phenotype", "Control", "EoE"))
res2 <- results(deseq_obj2, contrast = c("Phenotype", "Control", "FS-EoE"))
res3 <- results(deseq_obj2, contrast = c("Phenotype", "EoE", "FS-EoE"))
res3 <- results(deseq_obj2, contrast = c("Phenotype", "EoE", "FS-EoE"))
#age_something <- results(deseq_obj2, contrast = c( "Phenotype_EoE_vs_control", "Age..years."))

Age_DE <- results(deseq_obj2, name = c("Age..years."))
SV1_DE <- results(deseq_obj2, name = c("SV1"))
SV2_DE <- results(deseq_obj2, name = c("SV2"))
SV3_DE <- results(deseq_obj2, name = c("SV3"))
SV4_DE <- results(deseq_obj2, name = c("SV4"))


res1_control_EoE = as.data.frame((res1))
res2_control_FSEoE = as.data.frame((res2))
res3_EoE_FSEoE = as.data.frame((res3))
age_de = as.data.frame((Age_DE))
sv1_de = as.data.frame((SV1_DE))
sv2_de = as.data.frame((SV2_DE))
sv3_de = as.data.frame((SV3_DE))
sv4_de = as.data.frame((SV4_DE))


## Order by adjusted p-value
res1_control_EoE <- res1_control_EoE[order(res1_control_EoE$padj), ]
res2_control_FSEoE <- res2_control_FSEoE[order(res2_control_FSEoE$padj), ]

age_de <- age_de[order(age_de$padj), ]
sv1_de <- sv1_de[order(sv1_de$padj), ]
sv2_de <- sv2_de[order(sv2_de$padj), ]
sv3_de <- sv3_de[order(sv3_de$padj), ]
sv4_de <- sv4_de[order(sv4_de$padj), ]



top_control_EoE = subset(res1_control_EoE, padj < .1)
top_control_FSEoE = subset(res2_control_FSEoE, padj < .1)
top_EoE_FSEoE = subset(res3_EoE_FSEoE, padj < .1)
top_age = subset(age_de, padj < .1)
top_sv1 = subset(sv1_de, padj < .1)
top_sv2 = subset(sv2_de, padj < .1)
top_sv3 = subset(sv3_de, padj < .1)
top_sv4 = subset(sv4_de, padj < .1)


save(deseq_obj2, meta_data, EoE2_raw_counts, file = "~/Projects/Calies_Menard-Katcher/Code/DESeq_objects_calies_postfiltering.RData")
load("DESeq_objects_calies_postfiltering.RData")

session_info = sessionInfo()

write.csv(top_control_EoE, file = 'top_control_EoE.csv', quote = FALSE )
write.csv(top_control_FSEoE, file = 'top_control_FSEoE.csv', quote = FALSE )
write.csv(top_EoE_FSEoE, file = 'top_EoE_FSEoE.csv', quote = FALSE )
write.csv(top_age, file = 'top_age.csv', quote = FALSE )

gene_list_Age = rownames(top_age)
gene_list_cont_EoE = rownames(top_control_EoE)
gene_list_cont_FSEoE = rownames(top_control_FSEoE)
gene_list_EoE_FSEoE = rownames(top_EoE_FSEoE)
gene_list_All = rownames(rlog_t2)
  
write.table(gene_list_Age, file ='gene_list_Age_de.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep= " ")
write.table(gene_list_cont_EoE, file ='gene_list_cont_EoE_de.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep= " ")
write.table(gene_list_cont_FSEoE, file ='gene_list_cont_FSEoE_de.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep= " ")
write.table(gene_list_EoE_FSEoE, file ='gene_list_EoE_FSEoE_de.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep= " ")
write.table(gene_list_All, file ='gene_list_all.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep= " ")

#full ranked list
age_de_all = rownames(age_de)
write.table(age_de_all, file ='gene_list_Age_de_all.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep= " ")

control_EoE_all = rownames(res1_control_EoE)
write.table(control_EoE_all, file ='control_EoE_all.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep= " ")


par(mfrow=c(2,2))

gene_1 = plot(y = rlog_table['BTBD17',], x = meta_data$Phenotype, main = 'BTBD17') #, ylim = c(-1.928,-1.920))
gene_2 = plot(y = rlog_table['LINC00298',], x = meta_data$Phenotype, main ='LINC00298')
gene_3 = plot(y = rlog_table['FAM186A',], x = meta_data$Phenotype, main ='FAM186A')
gene_4 = plot(y = rlog_table['RPS17',], x = meta_data$Phenotype, main ='RPS17')
gene_5 = plot(y = rlog_table['TNFSF12-TNFSF13',], x = meta_data$Phenotype, main ='TNFSF12-TNFSF13')
gene_6 = plot(y = rlog_table['FGFBP1',], x = meta_data$Phenotype, main ='FGFBP1')
gene_7 = plot(y = rlog_table['IGLL1',], x = meta_data$Phenotype, main ='IGLL1')
gene_8 = plot(y = rlog_table['MIR643',], x = meta_data$Phenotype, main ='MIR643')
gene_9 = plot(y = rlog_table['CDC42P3',], x = meta_data$Phenotype, main ='CDC42P3')
gene_10 = plot(y = rlog_table['SNORD26',], x = meta_data$Phenotype, main ='SNORD26')
gene_11 = plot(y = rlog_table['KIR3DX1',], x = meta_data$Phenotype, main ='KIR3DX1')
gene_12 = plot(y = rlog_table['KCNS1',], x = meta_data$Phenotype, main ='KCNS1')
gene_13 = plot(y = rlog_table['ABHD16B',], x = meta_data$Phenotype, main ='ABHD16B')
gene_14 = plot(y = rlog_table['PLK2',], x = meta_data$Phenotype, main ='PLK2')
gene_15 = plot(y = rlog_table['TUNAR',], x = meta_data$Phenotype, main ='TUNAR')
?plot



rlog_counts_new2 <- rlog_table[is.element(rownames(rlog_table), rownames(top_control_EoE)),]
library('pheatmap')
pheatmap(rlog_counts_new2, colv = FALSE, scale = "row", labels_row = ' ', cluster_cols = FALSE,
         clustering_distance_rows = 'euclidean', clustering_method = 'complete')

library('writexl')

writexl::write_xlsx(
  list(EoE_vs_Control = res1, FS_vs_Control = res2, EoE_vs_FS = res3, Age = age_de),
  path = "my_output.xlsx")


#===========overlap genes=========#

# Load library
library(VennDiagram)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- c("blue", "orange")

venn.diagram(
  x = list(gene_list_cont_EoE, gene_list_cont_FSEoE),
  category.names = c("cont_vs_EoE" , "cont_vs_FSEoE"),
  filename = 'Con_EoE_vs_con_FSEoE_fullmodel.png',
  output=TRUE,
  fill = myCol,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 200,
  compression = "lzw",
  lwd = 2,
  lty = 0,
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055)

)

