library.pa

####load in mapped gff for each sample

#must figure out how these gffs are generated. Look at pipeline from project folder. look at what generates them and follow the rabbit hole
BRD4_tgfb_TSS = read.delim('/home/harrison/projects/McKinsey_tgf/mappedFolder/RN6_TSS_ALL_-3000_+0_brd4_tgfb.gff', header = TRUE)
BRD4_tgfb_TTR = read.delim('/home/harrison/projects/McKinsey_tgf/mappedFolder/RN6_TTR_ALL_-0_+3000_brd4_tgfb.gff', header = TRUE)
BRD4_tgfb_TXN = read.delim('/home/harrison/projects/McKinsey_tgf/mappedFolder/RN6_TXN_ALL_-0_+0_brd4_tgfb.gff', header = TRUE)

BRD4_veh_TSS = read.delim('/home/harrison/projects/McKinsey_tgf/mappedFolder/RN6_TSS_ALL_-3000_+0_brd4_veh.gff', header = TRUE)
BRD4_veh_TTR = read.delim('/home/harrison/projects/McKinsey_tgf/mappedFolder/RN6_TTR_ALL_-0_+3000_brd4_veh.gff', header = TRUE)
BRD4_veh_TXN = read.delim('/home/harrison/projects/McKinsey_tgf/mappedFolder/RN6_TXN_ALL_-0_+0_brd4_veh.gff', header = TRUE)

#TGFB_control_TSS = `RN6_TSS_ALL_.3000_+0_tgfb_control`
#TGFB_control_TTR = `RN6_TTR_ALL_.0_+3000_tgfb_control`
#TGFB_control_TXN = `RN6_TXN_ALL_.0_+0_tgfb_control`

#VEH_control_TSS = `RN6_TSS_ALL_.3000_+0_veh_control`
#VEH_control_TTR = `RN6_TTR_ALL_.0_+3000_veh_control`
#VEH_control_TXN = `RN6_TXN_ALL_.0_+0_veh_control`


#------------load in active genes list or whatever genes of interest you have

#active_genes = `RN6_NRVM_POL2_TSS_ACTIVE_.1000_+1000`

#------------ post intersecting with the file above (done on a separate script)
#active_genes = gained_NM_Genes

#active_genes = lost_NM_Genes

#----------New updated meta 6/26/19 with 400 genes from Pol2 dynamics

active_genes = Gained_genes
#active_genes = Lost_genes

#-----------combine, TSS, TXN, TTR, in order for each different sample/timepoint


BRD4_tgfb  = cbind(BRD4_tgfb_TSS[,3:62],BRD4_tgfb_TXN[,3:202],BRD4_tgfb_TTR[,3:62])
rownames(BRD4_tgfb) = BRD4_tgfb_TTR[,1]

BRD4_veh = cbind(BRD4_veh_TSS[,3:62], BRD4_veh_TXN[,3:202],BRD4_veh_TTR[,3:62])
rownames(BRD4_veh) = rownames(BRD4_tgfb)

#TGFB_control = cbind(TGFB_control_TSS[,3:62],TGFB_control_TXN[,3:202],TGFB_control_TTR[,3:62])
#rownames(TGFB_control) = rownames(BRD4_tgfb)

#VEH_control = cbind(VEH_control_TSS[,3:62], VEH_control_TXN[,3:202],VEH_control_TTR[,3:62])
#rownames(VEH_control) = rownames(BRD4_tgfb)

#---------------Active genes metas--------------------------------------------------


#--------filter for active genes/genes of interest

active0 = c()
for(i in 1:length(active_genes[,1])){
      row = which(as.character(rownames(BRD4_tgfb))==as.character(active_genes[i,1]))
      active0=c(active0,row)
}


# for each peak/region you will have as many values as you have genes in your list
# take the average for each region, to create vectors of the same length for each sample
# the length of these vectors should be the same as the number of bins you have in total


BRD4_tgfb_meta = apply(BRD4_tgfb[active0,],2,mean,na.rm=TRUE)
BRD4_veh_meta = apply(BRD4_veh[active0,],2,mean,na.rm=TRUE)
#TGFB_control_meta = apply(TGFB_control[active0,],2,mean,na.rm=TRUE)
#VEH_control_meta = apply(VEH_control[active0,],2,mean,na.rm=TRUE)



#---------------plot the meta plots-------------------------------------------------


pdf(file='/storage/cylin/grail/projects/McKinsey/meta/all_meta.pdf', width = 10, height = 8)
plot(1:320, BRD4_tgfb_meta, type='l', col='red', ylim = c(0,.8), xaxt='n', xlab='', ylab='rpm/bp', main='RN6_NRVM_POL2_TSS_ACTIVE_.1000_+1000')
lines(1:320, BRD4_veh_meta, type='l', col='black')
lines(1:320, TGFB_control_meta, type='l', col='black')
lines(1:320, VEH_control_meta, type='l', col='grey')
axis(1, c(0,60,260,320), c('-3kb', 'TSS', 'END', '+3kb'))
legend('topright', legend=c('Vehicle','TGF-B'),lty=1, col = c('black','red'), bty=0, cex = 1)

pdf(file='/storage/cylin/grail/projects/McKinsey/meta/all_meta.pdf', width = 10, height = 8)
plot(1:320, BRD4_tgfb_meta, type='l', col='red', ylim = c(0,.8), xaxt='n', xlab='', ylab='rpm/bp', main='Gained')
lines(1:320, BRD4_veh_meta, type='l', col='black')
lines(1:320, TGFB_control_meta, type='l', col='black')
lines(1:320, VEH_control_meta, type='l', col='grey')
axis(1, c(0,60,260,320), c('-3kb', 'TSS', 'END', '+3kb'))
legend('topright', legend=c('BRD4_tgfb','BRD4_veh'),lty=1, col = c('red','black'), bty=0, cex = 1)

pdf(file='/storage/cylin/grail/projects/McKinsey/meta/all_meta.pdf', width = 10, height = 8)
plot(1:320, BRD4_tgfb_meta, type='l', col='red', ylim = c(0,.8), xaxt='n', xlab='', ylab='rpm/bp', main='Pol II Dynamics *on top* (Genes with Decreased Enhancer-Bound BRD4+TGF-B)')
lines(1:320, BRD4_veh_meta, type='l', col='black')
lines(1:320, TGFB_control_meta, type='l', col='black')
lines(1:320, VEH_control_meta, type='l', col='grey')
axis(1, c(0,60,260,320), c('-3kb', 'TSS', 'END', '+3kb'))
legend('topright', legend=c('BRD4_tgfb','BRD4_veh'),lty=1, col = c('red','black'), bty=0, cex = 1)




#---------------new plots


pdf(file='/home/harrison/projects/McKinsey_tgf/gained_meta3.pdf', width = 10, height = 8)
plot(1:320, BRD4_tgfb_meta, type='l', col='red', ylim = c(0,4), xaxt='n', xlab='', ylab='rpm/bp', main='Gained')
lines(1:320, BRD4_veh_meta, type='l', col='black')
#lines(1:320, TGFB_control_meta, type='l', col='black')
#lines(1:320, VEH_control_meta, type='l', col='grey')
axis(1, c(0,60,260,320), c('-3kb', 'TSS', 'END', '+3kb'))
legend('topright', legend=c('BRD4_tgfb','BRD4_veh'),lty=1, col = c('red','black'), bty=0, cex = 1)

pdf(file='/home/harrison/projects/McKinsey_tgf/lost_meta3.pdf', width = 10, height = 8)
plot(1:320, BRD4_tgfb_meta, type='l', col='red', ylim = c(0,3), xaxt='n', xlab='', ylab='rpm/bp', main='Lost')
lines(1:320, BRD4_veh_meta, type='l', col='black')
#lines(1:320, TGFB_control_meta, type='l', col='black')
#lines(1:320, VEH_control_meta, type='l', col='grey')
axis(1, c(0,60,260,320), c('-3kb', 'TSS', 'END', '+3kb'))
legend('topright', legend=c('BRD4_tgfb','BRD4_veh'),lty=1, col = c('red','black'), bty=0, cex = 1)

