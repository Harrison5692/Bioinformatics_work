##################### Here we have a script that creates barplots of selected genes in different groups ################


#E2f4 bar plot

####### Loading raw table data into a new table object
#######

#---FPKM

#cuffnorm
mean_table <- cuffnorm_all_fpkm_exprs_raw


  
#--------TMM
mean_table <- x
  
#------DESeq2

mean_table <-  DESeq2_count
  
#CPM
mean_table <- x 
  

#raw
mean_table <- x
  
  

####### If we were using normalized data 
#######

#norm_table <- read.delim('/storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_all_fpkm_exprs_norm.txt', header=TRUE, sep = '\t')


####### Global Stuff (This was already here. Must just mean overall expression from the sample?)
#######

####### Here we input the genes of choice
#######

# my chosen
#gene_name = c('Postn', 'Sertad4', 'Map3k7', 'Brd4', 'Polr2a', 'Crebbp', 'Ep300','Trp53', 'Tgfb1', 'Tgfb', "Tgfb1i1", "Tgfbi","Tgfbr1","Tgfbr2","Tgfbr3","Tgfbe3l","Tgfbrap1", 'Tgfb2', "Tgfb3", "Smad4", "Smad3", "Smad2", "Smad1", "Smad6", "Smad7", "Pdgfrl", "Pdgfrb", "Pdgfra", "Pdgfd", "Pdgfc", "Pdgfb", "Pdgfa", "Foxp3", "Fas", "Gata3", "Kras", "Vegfa", "Vegfb", "Mmp9", "Mmp2","Lif","Sox4","Zeb1","Zeb2","Twist1","Twist2")


# DE chosen 40
gene_name = row.names(cfounts_40)


####### here we make a pdf file of the data
#######

pdf('C:/Users/Harrison/Desktop/bargraph_cuffnorm2.pdf')

#pdf('C:/Users/Harrison/Desktop/test_bargraph_deseq2.pdf')
####### Not sure whats happening here
####### Par is for parameters and mfrow how something to do with splitting data
#######

par(mfrow=c(2,2))


####### Ok this function does all of the work
#######

####### For each object in the list, search for it in the 'mean' table (might need to change the name)
#######

for( gene in gene_name){ 
  gene_row = which(rownames(mean_table) == gene)
  
  
  ####### just to check
  #######
  
  print(gene_row)
  print(gene)
  
  
  ####### Here we create a 'container' objects for the following code
  ####### (not sure why there needs to be an if statement)
  #######
  
  if(length(gene_row) != 0){ 
    plot_DMSO = c()
    plot_TGFb =  c()
    plot_A485 =  c()
    plot_PFCBP1 = c()
    
    
    
    ####### Here we gather the replicates into the 'containers'
    #######
    
    plot_DMSO = c(mean_table[gene_row,1], mean_table[gene_row,2], mean_table[gene_row,3])
    plot_TGFb =  c(mean_table[gene_row, 4], mean_table[gene_row,5], mean_table[gene_row, 6])
    plot_A485 =  c(mean_table[gene_row, 7], mean_table[gene_row,8])
    plot_PFCBP1 = c(mean_table[gene_row, 9], mean_table[gene_row,10])
    
    
    ####### Now we find the mean of these seperate replicate groups
    #######
    
    mean_DMSO = mean(plot_DMSO)
    mean_TGFb = mean(plot_TGFb)
    mean_A485 = mean(plot_A485)
    mean_PFCBP1 = mean(plot_PFCBP1)
    
    ####### Standard Deviation
    #######
    
    sd_DMSO = sd(plot_DMSO)
    sd_TGFb = sd(plot_TGFb)
    sd_A485 = sd(plot_A485)
    sd_PFCBP1 = sd(plot_PFCBP1)
    
    ####### Standard error for error bars?
    #######
    
    se_DMSO = sd_DMSO/(sqrt(length(plot_DMSO)))
    se_TGFb = sd_TGFb/(sqrt(length(plot_TGFb)))
    se_A485 = sd_A485/(sqrt(length(plot_A485)))
    se_PFCBP1 = sd_PFCBP1/(sqrt(length(plot_PFCBP1)))
    
    
    ####### No idea whats happening from here just yet but it works
    #######
    
    maxy=max(mean_DMSO,mean_TGFb,mean_A485, mean_PFCBP1) + 2.5*max(sd_DMSO,sd_TGFb,sd_A485, sd_PFCBP1) 
    print(maxy) 
    
    centers <- barplot(height = c(mean_DMSO,mean_TGFb,mean_A485, mean_PFCBP1),
                       beside = true, las = 2,
                       ylim = c(0, maxy),
                       cex.names = 0.75, 
                       main = gene,
                       ylab = "Expression Signal",
                       border = "black", axes = TRUE)
    
    
    # 45 degree string rotation
    text(x = centers, y = par("usr")[3] - 1, srt = 45,
         adj = 1, labels = c('DMSO', 'TGFb','A485','PFCBP1'), xpd = TRUE)
    
    segments(centers, c(mean_DMSO,mean_TGFb,mean_A485, mean_PFCBP1) - c(se_DMSO,se_TGFb, se_A485, se_PFCBP1) * 2, centers,
             c(mean_DMSO,mean_TGFb,mean_A485, mean_PFCBP1) + c(se_DMSO,se_TGFb, se_A485, se_PFCBP1) * 2, lwd = 1.5)
    
    arrows(centers, c(mean_DMSO,mean_TGFb,mean_A485, mean_PFCBP1) - c(se_DMSO,se_TGFb, se_A485, se_PFCBP1) * 2, centers,
           c(mean_DMSO,mean_TGFb,mean_A485, mean_PFCBP1) + c(se_DMSO,se_TGFb, se_A485, se_PFCBP1) * 2, lwd = 1.5, angle = 90,
           code = 3, length = 0.05)
  }
}

dev.off()


