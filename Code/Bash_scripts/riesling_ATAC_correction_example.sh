#!/usr/bin/bash
cd /storage/cylin/grail/projects/Lagor/ATAC/riesling-pipeline/
#cp /grail/TONY/Bradner/Homo.sapiens/Multiple.myeloma/MM1.S/ATAC/20151005_5081/20150621-MM1S-GM2033_S1_R1.fastq.gz /grail/projects/riesling/fastq/MM1S_ATAC/MM1S_ATAC.PE1.fastq.gz
#cp /grail/TONY/Bradner/Homo.sapiens/Multiple.myeloma/MM1.S/ATAC/20151005_5081/20150621-MM1S-GM2033_S1_R2.fastq.gz /grail/projects/riesling/fastq/MM1S_ATAC/MM1S_ATAC.PE2.fastq.gz


#bowtie2 --end-to-end --sensitive --no-unal --no-discordant --mm --met-stderr --time -x /storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 /storage/cylin/grail/projects/Lagor/mouse/ATAC/PE_R1_mm10/mouse_liver_r1.PE1.fastq.gz -2 /storage/cylin/grail/projects/Lagor/mouse/ATAC/PE_R1_mm10/mouse_liver_r1.PE2.fastq.gz 2>/storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R1_mm10/mouse_liver_r1.bt2.log | samtools view -bS - >/storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R1_mm10/mouse_liver_r1.bt2.bam &

#bowtie2 --end-to-end --sensitive --no-unal --no-discordant --mm --met-stderr --time -x /storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 /storage/cylin/grail/projects/Lagor/mouse/ATAC/PE_R2_mm10/mouse_liver_r2.PE1.fastq.gz -2 /storage/cylin/grail/projects/Lagor/mouse/ATAC/PE_R2_mm10/mouse_liver_r2.PE2.fastq.gz 2>/storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R2_mm10/mouse_liver_r2.bt2.log | samtools view -bS - >/storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R2_mm10/mouse_liver_r2.bt2.bam &

#bowtie2 --end-to-end --sensitive --no-unal --no-discordant --mm --met-stderr --time -x /storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -1 /storage/cylin/grail/projects/Lagor/human/ATAC/PE_R1_hg38/human_liver_r1.PE1.fastq.gz -2 /storage/cylin/grail/projects/Lagor/human/ATAC/PE_R1_hg38/human_liver_r1.PE2.fastq.gz 2>/storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R1_hg38/human_liver_r1.bt2.log | samtools view -bS - >/storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R1_hg38/human_liver_r1.bt2.bam



#this line below does not work hence the code above. The line below produces the code above except it messed up the path to the bowtie index which ive corrected above
                             
#python 1-map-to-genome.py -i /storage/cylin/grail/projects/Lagor/mouse/ATAC/PE_R1_mm10/ -o /storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R1_mm10/ -g mm10 -v


#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R1_hg38/ -o /storage/cylin/grail/projects/Lagor/ATAC/riesling/filtered/PE_R1_hg38/ -g hg38 -v &

#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R1_mm10/ -o /storage/cylin/grail/projects/Lagor/ATAC/riesling/filtered/PE_R1_mm10/ -g mm10 -v &

#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/Lagor/ATAC/riesling/bams/PE_R2_mm10/ -o /storage/cylin/grail/projects/Lagor/ATAC/riesling/filtered/PE_R2_mm10/ -g mm10 -v 

python 3-call-peaks.py -i /storage/cylin/grail/projects/Lagor/ATAC/riesling/filtered/PE_R1_hg38/ -o /storage/cylin/grail/projects/Lagor/ATAC/riesling/peaks/PE_R1_hg38/ -g hs -v --skip-macs14 & 

python 3-call-peaks.py -i /storage/cylin/grail/projects/Lagor/ATAC/riesling/filtered/PE_R1_mm10/ -o /storage/cylin/grail/projects/Lagor/ATAC/riesling/peaks/PE_R1_mm10/ -g mm -v --skip-macs14 &


python 3-call-peaks.py -i /storage/cylin/grail/projects/Lagor/ATAC/riesling/filtered/PE_R2_mm10/ -o /storage/cylin/grail/projects/Lagor/ATAC/riesling/peaks/PE_R2_mm10 -g mm -v --skip-macs14 
