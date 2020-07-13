#!/usr/bin/bash

#This script runs a HISAT2, a genomic aligner.


cd /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/
#SBATCH --output=/storage/cylin/grail/slurm_out/hisat2_RASMC_RNA_0H_A_20190201_11h35m43s_%j.out # Standard output and error log
#SBATCH -e /storage/cylin/grail/slurm_out/hisat2_RASMC_RNA_0H_A_20190201_11h35m43s_%j.err # Standard output and error log
pwd; hostname; date



#===================
#PROCESSING RASMC_RNA_0H_A
echo "processing RASMC_RNA_0H_A"
hisat2 -p 16 --no-unal -x /storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Hisat2Index_ERCC/rn6_ercc -U /storage/cylin/grail/data/rasmc/Rattus_norvegicus_ERCC/cardiomyocyte/cardiomyoctye/RNA/20140328_1212/ATCACG-s_7_1_sequence.fastq.gz -S /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.sam
/usr/bin/samtools view -bS /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.sam > /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.bam
/usr/bin/samtools sort /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.bam /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.sorted
/usr/bin/samtools index /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.sorted.bam
rm /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.sam
rm /storage/cylin/grail/projects/HS_Ras_myc/RNA/bams_rn6_ercc/HS20140328_1212.RN6_ERCC.bam



