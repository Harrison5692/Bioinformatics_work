merge_genesets_new <- merge(merged_genesets$, merged_genesets[2], all=TRUE, no.dups = TRUE)

geneset_list = c(as.character(merged_genesets[,1]),as.character(merged_genesets[,2]),as.character(merged_genesets[,3]),
                 as.character(merged_genesets[,4]),as.character(merged_genesets[,5]),as.character(merged_genesets[,6]))

write.table(geneset_list, file='~/geneset_list.txt')

table <- gsea_DOCA_8_D8vsS8
table2 <- gsea_report_for_DOCA_8_1536192150770
table3 <- gsea_ITF_I8vsS8
table4 <- gsea_report_for_ITF_8_1536192150770
table5 <- gsea_Sham_I8vsS8
table6 <- GSEA_sham_D8vsS8

table[is.na(table)] <- 0
table2[is.na(table2)] <- 0
table3[is.na(table3)] <- 0
table4[is.na(table4)] <- 0
table5[is.na(table5)] <- 0
table6[is.na(table6)] <- 0


geneset_list=unique(geneset_list)

#DOCA SETS
inflams=which(as.numeric(table[,1])==1)
inv = which(as.numeric(table[,1])==2)
mes=which(as.numeric(table[,1])==3)

inflams2=which(as.numeric(table2[,1])==1)
inv2 = which(as.numeric(table2[,1])==2)
mes2=which(as.numeric(table2[,1])==3)

doca_flam = unique(c(as.character(table[inflams,2]),as.character(table2[inflams2,2])))
doca_inv = unique(c(as.character(table[inv,2]),as.character(table2[inv2,2])))
doca_mes = unique(c(as.character(table[mes,2]),as.character(table2[mes2,2])))


#ITF SETS
inflams3=which(as.numeric(table3[,1])==1)
inv3 = which(as.numeric(table3[,1])==2)
mes3=which(as.numeric(table3[,1])==3)

inflams4=which(as.numeric(table4[,1])==1)
inv4 = which(as.numeric(table4[,1])==2)
mes4=which(as.numeric(table4[,1])==3)

itf_flam = unique(c(as.character(table3[inflams3,2]),as.character(table4[inflams4,2])))
itf_inv = unique(c(as.character(table3[inv3,2]),as.character(table4[inv4,2])))
itf_mes = unique(c(as.character(table3[mes3,2]),as.character(table4[mes4,2])))


#SHAM SETS
inflams5=which(as.numeric(table5[,1])==1)
inv5 = which(as.numeric(table5[,1])==2)
mes5=which(as.numeric(table5[,1])==3)

inflams6=which(as.numeric(table6[,1])==1)
inv6 = which(as.numeric(table6[,1])==2)
mes6=which(as.numeric(table6[,1])==3)

sham_flam = unique(c(as.character(table5[inflams5,2]),as.character(table6[inflams6,2])))
sham_inv = unique(c(as.character(table5[inv5,2]),as.character(table6[inv6,2])))
sham_mes = unique(c(as.character(table5[mes5,2]),as.character(table6[mes6,2])))

inflam_list = unique(c(as.character(table[inflams,2]),as.character(table2[inflams2,2]),as.character(table3[inflams3,2])),
                     c(as.character(table4[inflams,2]),as.character(table5[inflams2,2]),as.character(table6[inflams3,2])))
inv_list = unique(c(as.character(table[inv,2]), as.character(table2[inv2,2]), as.character(table3[inv3,2])),
                  c(as.character(table4[inv,2]), as.character(table5[inv2,2]), as.character(table6[inv3,2])))
mes_list = unique(c(as.character(table[mes,2]), as.character(table2[mes2,2]), as.character(table3[mes3,2])),
                  c(as.character(table4[mes,2]), as.character(table5[mes2,2]), as.character(table6[mes3,2])))

listy_list=c(as.character(table[inv,2]),as.character(table[inflams,2]),as.character(table[mes,2]),
             as.character(table2[inv2,2]),as.character(table2[inflams2,2]),as.character(table2[mes2,2]),
             as.character(table3[inv3,2]),as.character(table3[inflams3,2]),as.character(table3[mes3,2]),
             as.character(table[inv4,2]),as.character(table4[inflams,2]),as.character(table4[mes,2]),
             as.character(table2[inv5,2]),as.character(table5[inflams2,2]),as.character(table5[mes2,2]),
             as.character(table3[inv6,2]),as.character(table6[inflams3,2]),as.character(table6[mes3,2]))

listy_list=unique(listy_list)


write.table(listy_list,file='~/table_of_genesets.txt',sep='\t')

gen_mat <- matrix(nrow=403, ncol=9)
colnames(gen_mat) = c('DOCA_flam','DOCA_inv','DOCA_mes', 'ITF_flam','ITF_inv','ITF_mes','SHAM_flam','SHAM_inv','SHAM_mes')
rownames(gen_mat) = geneset_list

df_rows=c()
for(i in 1:length(doca_flam)){
  row=which(as.character(geneset_list)==as.character(doca_flam[i]))
  
  df_rows=c(df_rows,row)
}

gen_mat[df_rows,1]=1 

di_rows=c()
for(i in 1:length(doca_inv)){
  row=which(as.character(geneset_list)==as.character(doca_inv[i]))
  
  di_rows=c(di_rows,row)
}

gen_mat[di_rows,2]=2 

dm_rows=c()
for(i in 1:length(doca_mes)){
  row=which(as.character(geneset_list)==as.character(doca_mes[i]))
  
  dm_rows=c(dm_rows,row)
}

gen_mat[dm_rows,3]=3

if_rows=c()
for(i in 1:length(itf_flam)){
  row=which(as.character(geneset_list)==as.character(itf_flam[i]))
  
  if_rows=c(if_rows,row)
}

gen_mat[if_rows,4]=1 


ii_rows=c()
for(i in 1:length(itf_inv)){
  row=which(as.character(geneset_list)==as.character(itf_inv[i]))
  
  ii_rows=c(ii_rows,row)
}

gen_mat[ii_rows,5]=2 



im_rows=c()
for(i in 1:length(itf_mes)){
  row=which(as.character(geneset_list)==as.character(itf_mes[i]))
  
  im_rows=c(im_rows,row)
}

gen_mat[im_rows,6]=3 


sf_rows=c()
for(i in 1:length(sham_flam)){
  row=which(as.character(geneset_list)==as.character(sham_flam[i]))
  
  sf_rows=c(sf_rows,row)
}

gen_mat[sf_rows,7]=1 


si_rows=c()
for(i in 1:length(sham_inv)){
  row=which(as.character(geneset_list)==as.character(sham_inv[i]))
  
  si_rows=c(si_rows,row)
}

gen_mat[si_rows,8]=2 



sm_rows=c()
for(i in 1:length(sham_mes)){
  row=which(as.character(geneset_list)==as.character(sham_mes[i]))
  
  sm_rows=c(sm_rows,row)
}

gen_mat[sm_rows,9]=3 


gen_mat[is.na(gen_mat)] <- 0


gen_mat_inf = cbind(gen_mat[,1],gen_mat[,4],gen_mat[,7])
gen_mat_inv = cbind(gen_mat[,2],gen_mat[,5],gen_mat[,8])
gen_mat_mes = cbind(gen_mat[,3],gen_mat[,6],gen_mat[,9])

gen_mat_inf <- gen_mat_inf[rowSums(gen_mat_inf > 0) >= 1, ]
gen_mat_inv <- gen_mat_inv[rowSums(gen_mat_inv > 0) >= 1, ]
gen_mat_mes <- gen_mat_mes[rowSums(gen_mat_mes > 0) >= 1, ]

colnames(gen_mat_inf)=c('DOCA','ITF','SHAM')

colnames(gen_mat_inv)=c('DOCA','ITF','SHAM')

colnames(gen_mat_mes)=c('DOCA','ITF','SHAM')

ord = order(gen_mat_inf[,1])
ord2 = order(gen_mat_inv[,1])
ord3 = order(gen_mat_mes[,1])

heatmap(gen_mat_inf[ord,],scale = "none", Rowv = NA, Colv = NA, col = c('white','green'),main='Inflammatory genesets heatmap')
heatmap(gen_mat_inv[ord2,],scale = "none", Rowv = NA, Colv = NA, col = c('white','red'), main='Invasive genesets heatmap')
heatmap(gen_mat_mes[ord3,],scale = "none", Rowv = NA, Colv = NA, col = c('white','yellow'), main='Mesenchymal genesets heatmap')


heatmap(gen_mat[gs_order,],scale = "none", Rowv = NA, Colv = NA, col = c('white','black'))


gen_mat_vs <- matrix(nrow=403, ncol=3)
colnames(gen_mat_vs) = c('DOCAvSHAM','DOCAvITF','ITFvSHAM')
rownames(gen_mat_vs) = geneset_list


gen_mat_vs_inf <- matrix(nrow=403, ncol=3)
colnames(gen_mat_vs_inf) = c('DOCAvSHAM','DOCAvITF','ITFvSHAM')
rownames(gen_mat_vs_inf) = geneset_list

inflams = unique(inflams)

df_rows=c()
for(i in 1:length(inflams)){
  row=which(as.character(geneset_list)==as.character(table[inflams[i],2]))
  
  df_rows=c(df_rows,row)
}

gen_mat_vs_inf[df_rows,1]=1 
gen_mat_vs[df_rows,1]=1 

di_rows=c()
for(i in 1:length(inflams2)){
  row=which(as.character(geneset_list)==as.character(table2[inflams2[i],2]))
  
  di_rows=c(di_rows,row)
}

gen_mat_vs_inf[di_rows,2]=1 
gen_mat_vs[di_rows,2]=1

is_rows=c()
for(i in 1:length(inflams3)){
  row=which(as.character(geneset_list)==as.character(table3[inflams3[i],2]))
  
  is_rows=c(is_rows,row)
}

gen_mat_vs_inf[is_rows,3]=1
gen_mat_vs[is_rows,3]=1 

gen_mat_vs_inf[is.na(gen_mat_vs_inf)] <- 0

gen_mat_vs_inf <- gen_mat_vs_inf[rowSums(gen_mat_vs_inf > 0) >= 1, ]

ord = order(gen_mat_vs_inf[,1])


heatmap(gen_mat_vs_inf[ord,],scale = "none", Rowv = NA, Colv = NA, col = c('white','green'),main='Inflammatory genesets heatmap')








gen_mat_vs_inv <- matrix(nrow=403, ncol=3)
colnames(gen_mat_vs_inv) = c('DOCAvSHAM','DOCAvITF','ITFvSHAM')
rownames(gen_mat_vs_inv) = geneset_list

df_rows=c()
for(i in 1:length(inv)){
  row=which(as.character(geneset_list)==as.character(table[inv[i],2]))
  
  df_rows=c(df_rows,row)
}

gen_mat_vs_inv[df_rows,1]=1
gen_mat_vs[df_rows,1]=2

di_rows=c()
for(i in 1:length(inv2)){
  row=which(as.character(geneset_list)==as.character(table2[inv2[i],2]))
  
  di_rows=c(di_rows,row)
}

gen_mat_vs_inv[di_rows,2]=1 
gen_mat_vs[di_rows,2]=2

is_rows=c()
for(i in 1:length(inv3)){
  row=which(as.character(geneset_list)==as.character(table3[inv3[i],2]))
  
  is_rows=c(is_rows,row)
}

gen_mat_vs_inv[is_rows,3]=1 
gen_mat_vs[is_rows,3]=2

gen_mat_vs_inv[is.na(gen_mat_vs_inv)] <- 0

gen_mat_vs_inv <- gen_mat_vs_inv[rowSums(gen_mat_vs_inv > 0) >= 1, ]

ord = order(gen_mat_vs_inv[,1])


heatmap(gen_mat_vs_inv[ord,],scale = "none", Rowv = NA, Colv = NA, col = c('white','red'),main='Invasive genesets heatmap')





gen_mat_vs_mes <- matrix(nrow=403, ncol=3)
colnames(gen_mat_vs_mes) = c('DOCAvSHAM','DOCAvITF','ITFvSHAM')
rownames(gen_mat_vs_mes) = geneset_list

df_rows=c()
for(i in 1:length(mes)){
  row=which(as.character(geneset_list)==as.character(table[mes[i],2]))
  
  df_rows=c(df_rows,row)
}

gen_mat_vs_mes[df_rows,1]=1 
gen_mat_vs[df_rows,1]=3

di_rows=c()
for(i in 1:length(mes2)){
  row=which(as.character(geneset_list)==as.character(table2[mes2[i],2]))
  
  di_rows=c(di_rows,row)
}

gen_mat_vs_mes[di_rows,2]=1 
gen_mat_vs[di_rows,2]=3

is_rows=c()
for(i in 1:length(mes3)){
  row=which(as.character(geneset_list)==as.character(table3[mes3[i],2]))
  
  is_rows=c(is_rows,row)
}

gen_mat_vs_mes[is_rows,3]=1 
gen_mat_vs[is_rows,3]=3

gen_mat_vs_mes[is.na(gen_mat_vs_mes)] <- 0

gen_mat_vs_mes <- gen_mat_vs_mes[rowSums(gen_mat_vs_mes > 0) >= 1, ]

ord = order(gen_mat_vs_mes[,1])


heatmap(gen_mat_vs_mes[ord,],scale = "none", Rowv = NA, Colv = NA, col = c('white','yellow'),main='Mesenchymal genesets heatmap')


gen_mat_vs[is.na(gen_mat_vs)] <- 0

ord= order(gen_mat_vs[,1])

heatmap(gen_mat_vs[ord,],scale='none', Rowv = NA, Colv = NA, col = c('white','green','red','yellow'),main='All inf/inv/mes genesets heatmap')


gen_mat_vs_1 <- gen_mat_vs[rowSums(gen_mat_vs > 0) >= 1, ]

ord=order(gen_mat_vs_1[,1])
heatmap(gen_mat_vs_1[ord,],scale = "none", Rowv = NA, Colv = NA, col = c('white','green','red','yellow'),main='only inf/inv/mes genesets heatmap')












