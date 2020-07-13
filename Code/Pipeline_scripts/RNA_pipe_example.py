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



projectName = 'RNA'
dataFile = '/storage/cylin/grail/projects/HS_Ras_myc/data_tables/HS_RASMC_RNA_DATA_TABLE.txt' 

#project folders                                                                                                                                                                            
projectFolder = '/storage/cylin/grail/projects/HS_Ras_myc/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER                                                                                            

#standard folder names                                                                                                                                                                      
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
genome = 'rn6'
annotFile = '/storage/cylin/bin/pipeline/annotation/%s_refseq.ucsc' % (genome)

#making folders                                                                                                                                                                             
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)

#==========================================================================                                             
#=======================LOADING DATA ANNOTATION============================                                             
#==========================================================================                                             

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK                                        


#LOADING THE DATA TABLE                                                                                                 
#dataDict = pipeline_dfci.loadDataTable(dataFile)

#print(dataDict.keys())

#pipeline_dfci.summary(dataFile)

#print(dataDict)

#==========================================================================                                                                                                                       
#=======================HISAT2 ALIGNMENT FASTQ TO BAM======================
#========================================================================== 



#pipeline_dfci.mapHisat(dataFile,namesList=[],useSRA=False,pCount=16,Launch=True)

#==========================================================================                                                                                                                                 
#============================Cufflinks=====================================                                                                                                                                 
#========================================================================== 

gtfFile='/storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf'
analysisName='rasmc_rnaseq'
cufflinksFolder='/storage/cylin/grail/projects/HS_Ras_myc/RNA/cufflinks'

pipeline_dfci.makeCuffTableSlurm(dataFile,analysisName,gtfFile,cufflinksFolder,groupList=[['RASMC_RNA_0H_A','RASMC_RNA_0H_B'],['RASMC_RNA_PDGF_2H_B','RASMC_RNA_PDGF_2H_C','RASMC_RNA_PDGF_2H_D'],['RASMC_RNA_PDGF_JQ1_2H_E','RASMC_RNA_PDGF_JQ1_2H_G','RASMC_RNA_PDGF_JQ1_2H_H'],['RASMC_RNA_PDGF_24H_A','RASMC_RNA_PDGF_24H_B','RASMC_RNA_PDGF_24H_D'],['RASMC_RNA_PDGF_JQ1_24H_E','RASMC_RNA_PDGF_JQ1_24H_F','RASMC_RNA_PDGF_JQ1_24H_H']],bashFileName = '')
