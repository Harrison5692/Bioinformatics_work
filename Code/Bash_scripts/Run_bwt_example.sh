#!/usr/bin/bash
#SBATCH --output=/storage/cylin/grail/slurm_out/bwt2_NC248_20190806_11h35m31s_%j.out # Standard output and error log
#SBATCH -e /storage/cylin/grail/slurm_out/bwt2_NC248_20190806_11h35m31s_%j.err # Standard output and error log
pwd; hostname; date



#mkdir /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/
#cp /storage/genialis/bcm.genialis.com/data/24956/26-NC248_R1.fastq.gz /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.rawFastq.gz
#gunzip /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.rawFastq.gz
#bowtie2  -p 34 -x /storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -U /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.rawFastq -S /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.sam
#/bin/rm -f /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.rawFastq
#samtools view -bS '/storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.sam' > '/storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.bam'
#cd /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/
#samtools sort '/storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.bam' '/storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.hg19.bwt2.sorted'
#samtools index '/storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.hg19.bwt2.sorted.bam'
#/bin/rm -f '/storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.sam'
#mv /storage/cylin/grail/BOWTIE_TEMP/bwt2_NC248_8026/NC248.hg19.bwt2* /storage/cylin/grail/projects/moez/original_pipeline/bams/
ln /storage/cylin/grail/projects/moez/original_pipeline/bams/NC248.hg19.bwt2* /storage/cylin/grail/bam/hg19/
