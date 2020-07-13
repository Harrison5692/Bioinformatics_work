#!/usr/bin/bash
cd ~/riesling-pipeline/

#These steps are according to the Gordon Lab's reisling pipeline  (https://github.com/GordonLab/riesling-pipeline) to aid identifying open chromatin and or enhancers.

#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151020_5222/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_24H_REP1/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151020_5223/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_24H_REP2/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151019_5200/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_2H_REP2/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151020_5221/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_24H_REP2/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151020_5220/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_24H_REP1/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151019_5202/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_UNSTIM_REP2/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151019_5199/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_2H_REP1/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151020_5224/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_2H_REP1/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151020_5225/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_2H_REP2/ -g rn6 -v &
#python 1-map-to-genome.py -i /storage/cylin/grail/data/rasmc/ATAC/20151019_5201/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_UNSTIM_REP1/ -g rn6 -v &

#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_24H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_24H_REP1/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_24H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_24H_REP2/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_2H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_2H_REP2/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_24H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_24H_REP2/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_24H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_24H_REP1/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_UNSTIM_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_UNSTIM_REP2/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_2H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_2H_REP1/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_2H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_2H_REP1/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_PDGF_JQ1_2H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_2H_REP2/ -g rn6 -v &
#python 2-sanitize-bam.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/RASMC_ATAC_UNSTIM_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_UNSTIM_REP1/ -g rn6 -v &

##python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_24H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_JQ1_24H_REP1/ -g mm -v --skip-macs2
######python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_24H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_JQ1_24H_REP2/ -g mm -v --skip-macs2 &
##python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_2H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_2H_REP2/ -g mm -v --skip-macs2
######python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_24H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_24H_REP2/ -g mm -v --skip-macs2 &
##python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_24H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_24H_REP1/ -g mm -v --skip-macs2 &
##python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_UNSTIM_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_UNSTIM_REP2/ -g mm -v --skip-macs2 &
##python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_2H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_2H_REP1/ -g mm -v --skip-macs2 &
#python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_2H_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_JQ1_2H_REP1/ -g mm -v --skip-macs2 &
python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_PDGF_JQ1_2H_REP2/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_PDGF_JQ1_2H_REP2/ -g mm -v --skip-macs2
##python 3-call-peaks.py -i /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/RASMC_ATAC_UNSTIM_REP1/ -o /storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/RASMC_ATAC_UNSTIM_REP1/ -g mm -v --skip-macs2 &
