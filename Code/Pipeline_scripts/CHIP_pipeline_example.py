#!/usr/bin/python
#pipeline_template.py

'''
The MIT License (MIT)

Copyright (c) 2015 Charles Lin

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

#generic pipeline template for human data


#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys
sys.path.append('/storage/cylin/home/harrisos/pipeline/')

import pipeline_dfci
import utils
import string
import time
import datetime

from collections import defaultdict


#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'moez'
dataFile = '/storage/cylin/grail/projects/moez/original_pipeline/HEPG2_VPR_TABLE.txt' 

#project folders
projectFolder = '/storage/cylin/grail/projects/%s/original_pipeline' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names
gffFolder ='/%s/gff/' % (projectFolder)
macsFolder = '/%s/macsFolder/' % (projectFolder)
macsEnrichedFolder = '/%s/macsEnriched/' % (projectFolder)
mappedEnrichedFolder = '/%s/mappedEnriched/' % (projectFolder)
mappedFolder = '/%s/mappedFolder/' % (projectFolder)
wiggleFolder = '/%s/wiggles/' % (projectFolder)
metaFolder = '/%s/meta/' % (projectFolder)
genome = 'hg19'
annotFile = '/storage/cylin/bin/pipeline/annotation/hg19_refseq.ucsc'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#=======================LOADING DATA ANNOTATION============================
#==========================================================================

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK


#LOADING THE DATA TABLE
dataDict = pipeline_dfci.loadDataTable(dataFile)

print(dataDict.keys())

pipeline_dfci.summary(dataFile)

print(dataDict)

#==========================================================================
#==========================CALLING BOWTIE==================================
#==========================================================================

#THIS SECTION CALLS BOWTIE ON RAW READ FILES TO GENERATE SORTED AND INDEXED BAMS IN THE BAM FOLDER

#namesList = [] <- fill this in if you want to only map a subset of the data. otherwise leave blank
namesList = []


#SET LAUNCH TO False to debug
pipeline_dfci.makeBowtieBashJobsSlurm(dataFile,namesList,launch=True,overwrite=True,pCount=34)

#==========================================================================
#=============================CALL MACS====================================
#==========================================================================
#namesList = []

#print(namesList)
#pipeline_dfci.callMacsSlurm(dataFile,macsFolder,namesList,overwrite=True,pvalue='1e-9')

#==========================================================================
#=======================FORMAT MACS OUTPUT=================================
#==========================================================================

#pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder, wigLink='')\

#==========================================================================
#===========================CALLING ROSE===================================
#==========================================================================
'''
namesList=['PRIMARY_ML2811_H3K27AC_A', 'PRIMARY_ML2811_H3K27AC_B', 'PRIMARY_ML2738_H3K27AC', 'PRIMARY_ML3040_H3K27AC']

parentFolder = utils.formatFolder('%s/rose' % (projectFolder),True)
bashFileName = '%sk27ac_rose_2500.sh' % (parentFolder)
pipeline_dfci.callRose2Slurm(dataFile,macsEnrichedFolder,parentFolder,namesList,extraMap = [],inputFile='',tss=2500,stitch='',bashFileName =bashFileName,mask='')
'''
#==========================================================================

#===================MAKE TRACK HUB=========================================

#==========================================================================

#analysis_name = ''

#dataFileList = [dataFile]

#chrom_sizes = '/storage/cylin/grail/genomes/chrom_sizes/hg19.chrom.sizes'

#project_folder = projectFolder

#pipeline_dfci.makeTrackHub(analysis_name,project_folder,chrom_sizes,dataFileList, wiggle_dir='',web_dir='/storage/cylin/web/Lin_Lab_Track_Hubs/',hub_name='',hub_short_lab='',hub_long_lab='',EMAIL='',fileType='bigWig',col='0,0,0')

