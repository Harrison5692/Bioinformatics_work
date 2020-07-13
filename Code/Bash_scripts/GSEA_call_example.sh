#!/usr/bin/bash

#This script runs gene set enrichment analysis on the terminal 

source activate py27_anaconda
#COMMAND LINE GSEA CALLS 
#metric possibilities include Signal2Noise, tTest, Ratio_of_Classes, Diff_of_Classes, log2_Ratio_of_Classes
#-Xmx means memory


java -Xmx4000m -cp /storage/cylin/home/cl6/gsea2-3.0_beta_2.jar xtools.gsea.Gsea -res /storage/cylin/grail/projects/HS_Ras_myc/RNA/cufflinks2/rasmc_rnaseq_cuffnorm/output/rasmc_rnaseq_RASMC_RNA_0H_vs_RASMC_RNA_PDGF_2H.gct -cls /storage/cylin/grail/projects/HS_Ras_myc/RNA/cufflinks2/rasmc_rnaseq_cuffnorm/output/rasmc_rnaseq_RASMC_RNA_0H_vs_RASMC_RNA_PDGF_2H.cls#RASMC_RNA_PDGF_2H_versus_RASMC_RNA_0H -gmx /storage/cylin/grail/annotations/gsea/c2.all.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label RASMC_RNA_PDGF_2H_versus_RASMC_RNA_0H -metric log2_Ratio_of_Classes -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -o /storage/storage/cylin/grail/projects/HS_Ras_myc/RNA/GSEA/output/ -gui false



