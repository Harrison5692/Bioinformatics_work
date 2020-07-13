#!/usr/bin/bash

# This is for ranking enhancer and super enhancer activity based on most to least H3K27ac ChIP signal

#SBATCH --output=/storage/cylin/grail/slurm_out/ROSE2_HS_RASMC_CHIP_DATA_TABLE_ROSE2_20190307_15h33m33s_%j.out # Standard output and error log
#SBATCH -e /storage/cylin/grail/slurm_out/ROSE2_HS_RASMC_CHIP_DATA_TABLE_ROSE2_20190307_15h33m33s_%j.err # Standard output and error log
#SBATCH -n 8
#SBATCH --mem 32768
pwd; hostname; date



cd /storage/cylin/home/harrisos/pipeline
#python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_POL2_PDGF_2H_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1184.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_POL2_PDGF_2H_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1161.rn6.bwt2.sorted.bam
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_POL2_UNSTIM_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1165.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_POL2_UNSTIM_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1197.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_POL2_PDGF_24H_JQ1_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1189.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_POL2_PDGF_24H_JQ1_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1182.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_POL2_PDGF_24H_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1186.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_POL2_PDGF_24H_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1172.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_POL2_PDGF_2H_JQ1_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1187.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_POL2_PDGF_2H_JQ1_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1166.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_BRD4_PDGF_24H_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1190.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_BRD4_PDGF_24H_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1183.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_BRD4_PDGF_24H_REP2_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1174.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_BRD4_PDGF_24H_REP2_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1199.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_BRD4_PDGF_2H_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1188.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_BRD4_PDGF_2H_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1169.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_BRD4_PDGF_2H_REP2_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1198.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_BRD4_PDGF_2H_REP2_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1168.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_BRD4_UNSTIM_NEW_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1194.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_BRD4_UNSTIM_NEW_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1180.rn6.bwt2.sorted.bam &
python ROSE2_main.py -g RN6 -i /storage/cylin/grail/projects/HS_rasmc/CHIP/macsEnriched/RASMC_BRD4_UNSTIM_REP1_peaks.bed -r /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1178.rn6.bwt2.sorted.bam -o /storage/cylin/grail/projects/HS_rasmc/CHIP/rose/RASMC_BRD4_UNSTIM_REP1_ROSE -t 2500 -c /storage/cylin/grail/projects/HS_rasmc/CHIP/bams/1170.rn6.bwt2.sorted.bam
