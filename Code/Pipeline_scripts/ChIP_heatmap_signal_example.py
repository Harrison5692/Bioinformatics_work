#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2018 Charles Lin

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


#Main method run script for processing of slam seq analysis from Muhar et al., 2018




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

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



projectName = 'McKinsey'
genome ='rn6'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
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

#mask Files


#genomeDirectory #select your genome
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/hg19/Sequence/Chromosomes/'

#making folders
#folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

#for folder in folderList:
#    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Seq
chip_data_file = '%sdata_tables/MCKINSEY_BRD4_DATA_TABLE.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(chip_data_file)


    print('\n\n')
    print('#======================================================================')
    print('#======================II. MAPPING FOR HEATMAP Enhancers===============')
    print('#======================================================================')
    print('\n\n')

    cellTypeList = ['brd4','tgfb','veh']
    gffList = ['/storage/cylin/grail/projects/McKinsey/gff/All_enhancers_20kb_regions.gff']
    
    mappedFolder_n200 = utils.formatFolder('%smappedFolder_n200' % (projectFolder),True)
    
    
    #pipeline_dfci.mapBams(chip_data_file,cellTypeList,gffList,mappedFolder_n200,nBin = 200,overWrite =False,rpm=True,nameList = [],extension=150)

    geneListFile =''
    heatFolder = utils.formatFolder('%sheatmaps_HS2/' % (projectFolder),True)
    namesList = ['brd4_tgfb','brd4_veh']
    orderByName = 'brd4_veh'

    pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[0],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)

    print('\n\n')
    print('#======================================================================')
    print('#======================II. MAPPING FOR HEATMAP Enhancers===============')
    print('#======================================================================')
    print('\n\n')

    cellTypeList = ['brd4','tgfb','veh']
    gffList = ['/storage/cylin/grail/projects/McKinsey/gff/All_enhancers_20kb_conserved.gff','/storage/cylin/grail/projects/McKinsey/gff/All_enhancers_20kb_veh_peaks.gff','/storage/cylin/grail/projects/McKinsey/gff/All_enhancers_20kb_tgfb_peaks.gff']

    mappedFolder_n200 = utils.formatFolder('%smappedFolder_n200' % (projectFolder),True)


#    pipeline_dfci.mapBams(chip_data_file,cellTypeList,gffList,mappedFolder_n200,nBin = 200,overWrite =False,rpm=True,nameList = [],extension=150)

    geneListFile =''
    heatFolder = utils.formatFolder('%sheatmaps_HS2/' % (projectFolder),True)
    namesList = ['brd4_tgfb','brd4_veh']
    orderByName = 'brd4_veh'

    pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[0],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)
    pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[1],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)
    pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[2],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)

    orderByName='brd4_tgfb'
    pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[2],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)



    print('\n\n')
    print('#======================================================================')
    print('#======================II. MAPPING FOR HEATMAP Promoters===============')
    print('#======================================================================')
    print('\n\n')

    cellTypeList = ['brd4','tgfb','veh']
    gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000_BRD4_conserved.gff','/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000_BRD4_veh.gff','/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000_BRD4_tgfb.gff']
    
    mappedFolder_n200 = utils.formatFolder('%smappedFolder_n200' % (projectFolder),True)
    
    
#    mapBams(chip_data_file,cellTypeList,gffList,mappedFolder_n200,nBin = 200,overWrite = False,rpm=True,nameList = [],extension=150)

    
    geneListFile =''
    heatFolder = utils.formatFolder('%sheatmaps_HS/' % (projectFolder),True)
    namesList = ['brd4_tgfb','brd4_veh']
    orderByName = 'brd4_veh'

    #pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[0],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)
    #pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[1],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)
    #pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[2],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)

    orderByName='brd4_tgfb'
    #pipeline_dfci.callHeatPlotOrdered(chip_data_file,gffList[2],namesList,orderByName,geneListFile,heatFolder,mappedFolder_n200,relative=False,useBackground=True)




############################################################################

def mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,rpm=True,nameList = [],extension=200):
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    if mappedFolder[-1] != '/':
        mappedFolder+='/'
    ticker = 0
    for gffFile in gffList:

        #check to make sure gff exists
        try:
            foo = open(gffFile,'r')
        except IOError:
            print('ERROR: GFF FILE %s DOES NOT EXIST' % (gffFile))

        gffName = gffFile.split('/')[-1].split('.')[0]

        #see if the parent directory exists, if not make it
        try:
            foo = os.listdir(mappedFolder)
        except OSError:
            print('%s directory not found for mapped bams. making it' % (mappedFolder))
            os.system('mkdir %s' % (mappedFolder))



        #make this directory specifically for the gff
        try:
            foo = os.listdir(mappedFolder+gffName)
        except OSError:
            print('%s directory not found for this gff: %s. making it' % (mappedFolder+gffName,gffName))
            os.system('mkdir %s%s' % (mappedFolder,gffName))

        outdir = mappedFolder+gffName+'/'

        if len(nameList) == 0:
            nameList = dataDict.keys()



        for name in nameList:
            print ('mapping %s to %s' % (name,gffFile))
            #filter based on celltype
            cellName = name.split('_')[0]
            if cellTypeList.count(cellName) != 1:
                print("this guy didn't get mapped %s" % (name))
                continue
            fullBamFile = dataDict[name]['bam']
            outFile = outdir+gffName+'_'+name+'.gff'



            if overWrite:
                cmd1 = "python2.7 %s/bamToGFF_turbo.py -e %s -m %s -b %s -i %s -o %s" % ('/storage/cylin/bin/pipeline',extension,nBin,fullBamFile,gffFile,outFile)
                if rpm:
                    cmd1 += ' -r'
                #cmd1 += ' &'

                print(cmd1)
                os.system(cmd1)

            else:
                try:
                    Foo = open(outFile,'r')
                    print('File %s Already Exists, not mapping' % (outFile))
                except IOError:
                    cmd1 = "python2.7 %s/bamToGFF_turbo.py -e %s -m %s -b %s -i %s -o %s" % ('/storage/cylin/bin/pipeline',extension,nBin,fullBamFile,gffFile,outFile)
                    if rpm:
                        cmd1 += ' -r'
                    #cmd1 += ' &'

                    print(cmd1)
                    os.system(cmd1)






#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()




