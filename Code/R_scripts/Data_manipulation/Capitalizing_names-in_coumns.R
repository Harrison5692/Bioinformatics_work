###capitalizing names in each row of a column#####################


cuffnorm_Group1_vs_Group2[,1] <- toupper(cuffnorm_Group1_vs_Group2[,1])

cuffnorm_Group1_vs_Group3[,1] <- toupper(cuffnorm_Group1_vs_Group3[,1])

cuffnorm_Group1_vs_Group4[,1] <- toupper(cuffnorm_Group1_vs_Group4[,1])

cuffnorm_Group2_vs_Group3[,1] <- toupper(cuffnorm_Group2_vs_Group3[,1])

cuffnorm_Group2_vs_Group4[,1] <- toupper(cuffnorm_Group2_vs_Group4[,1])

cuffnorm_Group3_vs_Group4[,1] <- toupper(cuffnorm_Group3_vs_Group4[,1])

write.table(cuffnorm_Group1_vs_Group2, file = "C:/Users/Harrison/Desktop/GSEA_inputs_drop7_8/cuffnorm_Group1_vs_Group2.gct", row.names = FALSE, quote = FALSE, sep = "\t" )

write.table(cuffnorm_Group1_vs_Group3, file = "C:/Users/Harrison/Desktop/GSEA_inputs_drop7_8/cuffnorm_Group1_vs_Group3.gct", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(cuffnorm_Group1_vs_Group4, file = "C:/Users/Harrison/Desktop/GSEA_inputs_drop7_8/cuffnorm_Group1_vs_Group4.gct", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(cuffnorm_Group2_vs_Group3, file = "C:/Users/Harrison/Desktop/GSEA_inputs_drop7_8/cuffnorm_Group2_vs_Group3.gct", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(cuffnorm_Group2_vs_Group4, file = "C:/Users/Harrison/Desktop/GSEA_inputs_drop7_8/cuffnorm_Group2_vs_Group4.gct", row.names = FALSE, quote = FALSE, sep = "\t")

write.table(cuffnorm_Group3_vs_Group4, file = "C:/Users/Harrison/Desktop/GSEA_inputs_drop7_8/cuffnorm_Group3_vs_Group4.gct", row.names = FALSE, quote = FALSE, sep = "\t")
