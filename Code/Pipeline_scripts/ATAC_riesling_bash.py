#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2017 YOUR NAME HERE and Charles Lin lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for processing of YOUR PROJECT HERE




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/home/harrisos/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re

from collections import defaultdict
import subprocess
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'HS_rasmc/ATAC/'
genome ='rn6'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%s/riesling/peaks/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%s/riesling/wiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#atac-seq
#atac_data_file = '%sdata_tables/RASMC_ATAC_TABLE.txt' % (projectFolder)
atac_data_file = '/storage/cylin/grail/projects/HS_rasmc/data_tables/HS_RASMC_ATAC_TABLE.txt'



#This section sanity checks each data table and makes sure both bam and .bai files are accessible

dataDict=pipeline_dfci.loadDataTable(atac_data_file)

#for data file
pipeline_dfci.summary(atac_data_file)

#assumes macs has already been run and formatted
#    run_macs(chip_data_file)???

#    sys.exit()

#===============================================================
#=================ATAC bash attempt1 ===========================

namesList=dataDict.keys()


#temp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
#randTicker = random.randint(o,10000)



bashFileName= '%s1riesling_ATAC_temp.sh' % (projectFolder)
bashFile = open(bashFileName, 'w')

#ts = time.time()
#timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')

cmd = '#!/usr/bin/bash'
bashFile.write(cmd+'\n')

#cmd = '#SBATCH --output=/storage/cylin/grail/slurm_out/ATAC_%s_%s' % (projectName, timestamp) + '_%j.out # Standard output and error log'
#bashFile.write(cmd+'\n')

#cmd = '#SBATCH -e /storage/cylin/grail/slurm_out/ATAC_%s_%s' % (projectName,timestamp) + '_%j.err # Standard output and error log'
#bashFile.write(cmd+'\n')

#cmd = 'pwd; hostname; date'
#bashFile.write(cmd+'\n\n\n\n')

bashFile.write('cd %s' % (projectFolder))
bashFile.write('\n\n')


for name in namesList:

    fq_path=dataDict[name]['fastq']
    uid=dataDict[name]['uniqueID']
    path=fq_path.split(uid)[0]
    bam_dir='/storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/'+ name +'/'
    cmd='#python 1-map-to-genome.py -i %s -o %s -g rn6 -v' % (path,bam_dir)

    bashFile.write(cmd+' &\n')

    
for name in namesList:
    
    bam_dir='/storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/bams/'+ name +'/'
    filt_bam='/storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/' + name +'/'
    cmd= '#python 2-sanitize-bam.py -i %s -o %s -g rn6 -v' % (bam_dir, filt_bam)

    bashFile.write(cmd+' &\n')

for name in namesList:
    
    filt_bam='/storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/filtered/bams/' + name +'/'
    peaks='/storage/cylin/grail/projects/HS_rasmc/ATAC/riesling/peaks/' + name + '/'
    cmd='#python 4-call-peaks.py -i %s -o %s -g mm -v --skip-macs2' % (filt_bam, peaks)
    bashFile.write(cmd+' &\n')

    

bashFile.close()


#=====================Plot bams=======================================
'''
def main():

 # for BRD4                                                                                                                                                
                                              
    dataFile = chip_data_file

    figureGFFPath = '/storage/cylin/grail/projects/HS_rasmc/CHIPrx/gff/figure_1.gff'

    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)



    plotName = 'rasmc_all_figure_1_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=Tru\
e,rxGenome = '')                 

#=======================================END========================================


if __name__=="__main__":
    main()
'''
