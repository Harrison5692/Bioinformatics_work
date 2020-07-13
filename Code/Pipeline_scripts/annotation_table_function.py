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
sys.path.append('/storage/cylin/home/cl6/pipeline/')

import pipeline_dfci
import utils
import string
import time
import datetime

from collections import defaultdict


#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'McKinsey'
pol2_dataFile = '/storage/cylin/grail/projects/%s/MCKINSEY_POL2_DATA_TABLE.txt' % (projectName)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names

#making folders
#folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

#for folder in folderList:
#    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#=======================LOADING DATA ANNOTATION============================
#==========================================================================

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK


#LOADING THE DATA TABLE
#dataDict = pipeline_dfci.loadDataTable(dataFile)

#print(dataDict.keys())

#pipeline_dfci.summary(dataFile)

#=============================================================
#=================FORMATTING GENIALIS TABLE===================
#=============================================================


def makeGenialisTable(dataFile,outFilePath,organism='',seqType='',paired=False,collection='',annotator='',source='',strain='',tissue='',age='',genotype='',molecule='',libraryStrategy='',exPro='',libraryConst='',other1='',other2=''):

    genialis_header = ['NAME','FASTQ_R1','FASTQ_R2','SEQ_TYPE','PAIRED','COLLECTION','ANNOTATOR','SOURCE','ORGANISM','STRAIN','TISSUE','AGE','GENOTYPE','MOLECULE',
                       'LIBRARY_STRATEGY','EXTRACTION_PROTOCOL','LIBRARY_CONSTRUCTION_PROTOCOL','OTHER_CHAR_1','OTHER_CHAR_2']

    genomeDict = {'HG38':'Homo sapiens',
                  'HG19':'Homo sapiens',
                  'HG19_ERCC':'Homo sapiens',
                  'MM9':'Mus musculus',
                  'MM10':'Mus musculus',
                  'RN6':'Rattus norvegicus',
                  'RN6_ERCC':'Rattus norvegicus',
                           }

    
    genialisTable=[]
    genialisTable.append(genialis_header)
    
    dataTable = utils.parseTable(dataFile,'\t')
    
    

    for line in dataTable[1:]:
        name = line[3]
        if paired==False:
            fastq1 = line[8]
            fastq2 = ''
            pair=0
        elif paired==True:
            foo = line[8].split('::')
            print(foo)
            fastq1 = foo[0]
            fastq2 = foo[1]
            pair=1
        if organism=='':
            organism=genomeDict[line[2].upper()]
        
        new_line = [name,fastq1,fastq2,seqType,pair,collection,annotator,source,organism,strain,tissue,age,
                    genotype,molecule,libraryStrategy,exPro,libraryConst,other1,other2]

        genialisTable.append(new_line)

    utils.unParseTable(genialisTable,outFilePath,'\t')



outFilePath1 = '%s%s_brd4_genialis_annotation_table.txt' % (projectFolder,projectName)
outFilePath2 = '%s%s_pol2_genialis_annotation_table.txt' % (projectFolder,projectName)

dataFile = pol2_dataFile

makeGenialisTable(dataFile,outFilePath2,organism='',seqType='CHiPseq',paired=False,collection='McKinsey',annotator='Rachel Hirsch',source='', tissue='',age='',genotype='',molecule='',libraryStrategy='',exPro='',libraryConst='',other1='',other2='')
