#!/usr/bin/bash
#SBATCH --output=/storage/cylin/grail/slurm_out/ROSE2_mm10_liver_CHIP_data_table_ROSE2_20190714_18h03m10s_%j.out # Standard output and error log
#SBATCH -e /storage/cylin/grail/slurm_out/ROSE2_mm10_liver_CHIP_data_table_ROSE2_20190714_18h03m10s_%j.err # Standard output and error log
#SBATCH -n 8
#SBATCH --mem 32768
pwd; hostname; date



cd /storage/cylin/bin/pipeline
python2 ROSE2_META.py -g MM10 -i /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/macsEnriched/mm10_27ac_liver_R1.bed,/storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/macsEnriched/mm10_27ac_liver_R2.bed -r /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/bams/ENCFF001KMH.nodup.bam,/storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/bams/ENCFF001KMI.nodup.bam -o /storage/cylin/grail/projects/Lagor/meta/mm10_27ac_liver_ROSE_merged -t 2500 --mask /storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Annotation/Masks/mm10_blacklist.bed



#python2 ROSE2_META.py -g MM10 -i /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/macsEnriched/mm10_27ac_liver_R2.bed -r /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/bams/ENCFF001KMI.nodup.bam -o /storage/cylin/grail/projects/Lagor/rose/meta/mm10_27ac_liver_R2_ROSE -t 2500 --mask /storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Annotation/Masks/mm10_blacklist.bed


#python2 ROSE2_META.py -g MM10 -i /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/macsEnriched/mm10_pol2_liver_R1.bed -r /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/bams/ENCFF001LNH.nodup.bam -o /storage/cylin/grail/projects/Lagor/rose/meta/mm10_pol2_liver_R1_ROSE -t 2500 --mask /storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Annotation/Masks/mm10_blacklist.bed
#python2 ROSE2_META.py -g MM10 -i /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/macsEnriched/mm10_pol2_liver_R2.bed -r /storage/cylin/grail/projects/Lagor/aquas_pipe_mouse_CHIP/bams/ENCFF001LNR.nodup.bam -o /storage/cylin/grail/projects/Lagor/rose/meta/mm10_pol2_liver_R2_ROSE -t 2500 --mask /storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Annotation/Masks/mm10_blacklist.bed
