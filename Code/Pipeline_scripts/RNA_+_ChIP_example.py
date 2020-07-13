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
sys.path.append('/storage/cylin/bin/pipeline/')

import pipeline_dfci
import utils
import string
import time
import datetime
import subprocess

from collections import defaultdict


#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'McKinsey'
pol2_dataFile = '/storage/cylin/grail/projects/%s/MCKINSEY_POL2_DATA_TABLE.txt' % (projectName)
brd4_dataFile = '/storage/cylin/grail/projects/%s/MCKINSEY_BRD4_DATA_TABLE.txt' % (projectName)
rna_dataFile =  '/storage/cylin/grail/projects/%s/MCKINSEY_RNA_TABLE.txt' % (projectName)


genome ='rn6'
annotFile = '/storage/cylin/home/cl6/pipeline/annotation/%s_refseq.ucsc' % (genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#========================FORMATTING SAMPLE TABLE===========================
#==========================================================================

# #THIS SECTION CREATES A DATA TABLE FROM A WHITEHEAD ANNOTATION SPREADSHEET

# #give full path
# #sampleTableFile = 'YOUR_WIGTC_ANNOTATION.xls' #<- the .xls file in the seq data folder provided by WI

# dirpath = ''  <- provide full path of folder containing raw seq files
# #e.g. /ark/home/jr246/raw/130925_..../QualityScore/

# #bamPath <- where we store our bams.  Must have write access if you want to call bowtie
# #e.g. /ark/home/jr246/bam/
# bamPath = '/grail/bam/'

# pipeline_dfci.makePipelineTable(sampleTableFile,dirPath,bamPath,dataFile)

# dataDict = pipeline_dfci.loadDataTable(dataFile)

# namesList = dataDict.keys()

# print(namesList)


#==========================================================================
#=======================LOADING DATA ANNOTATION============================
#==========================================================================

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK


#LOADING THE DATA TABLE
#dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)
dataDict = pipeline_dfci.loadDataTable(pol2_dataFile)
#dataDict = pipeline_dfci.loadDataTable(rna_dataFile)

print(dataDict.keys())

pipeline_dfci.summary(pol2_dataFile)
#pipeline_dfci.summary(pol2_dataFile)
#pipeline_dfci.summary(rna_dataFile)


#==========================================================================
#==========================CALLING BOWTIE==================================
#==========================================================================

#THIS SECTION CALLS BOWTIE ON RAW READ FILES TO GENERATE SORTED AND INDEXED BAMS IN THE BAM FOLDER


#namesList = []  <- fill this in if you want to only map a subset of the data. otherwise leave blank
#namesList = []
#namesList = dataDict.keys()
#print(namesList)

##SET LAUNCH TO False to debug
#pipeline_dfci.makeBowtieBashJobs(brd4_dataFile,namesList,launch=True,overwrite=False,pCount=34)

# # # #THIS SECTION CALLS BOWTIE ON RAW READ FILES TO GENERATE SORTED AND INDEXED BAMS IN THE BAM FOLDER

#namesList = dataDict.keys()
#namesList = ['pol2_veh','pol2_tgfb']
# #SET LAUNCH TO False to debug
#pipeline_dfci.makeBowtieBashJobs(pol2_dataFile,namesList,launch=True,overwrite=False,pCount=34)

#==========================================================================
#=============================CALL MACS====================================
#==========================================================================

# THIS SECTION CALLS THE MACS ERROR MODEL

#dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)

#namesList = dataDict.keys()

#print(namesList)
#pipeline_dfci.callMacsSlurm(brd4_dataFile,macsFolder,namesList,overwrite=True,pvalue='1e-9')

################################################################################

#dataDict = pipeline_dfci.loadDataTable(pol2_dataFile)

#namesList = dataDict.keys()

#print(namesList)
#pipeline_dfci.callMacsSlurm(pol2_dataFile,macsFolder,namesList,overwrite=True,pvalue='1e-9')


#==========================================================================
#=======================FORMAT MACS OUTPUT=================================
#==========================================================================

#THIS SECTION FORMATS THE OUTPUT FROM MACS, CREATES THE MACSENRICHED FOLDER AND MOVES WIGGLES TO THE DESTINATION

#pipeline_dfci.formatMacsOutput(pol2_dataFile,macsFolder,macsEnrichedFolder,wiggleFolder, wigLink='')
#pipeline_dfci.formatMacsOutput(brd4_dataFile,macsFolder,macsEnrichedFolder,wiggleFolder, wigLink='')


#==========================================================================
#============================CALLING ROSE==================================
#==========================================================================

#dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)

#namesList = dataDict.keys()
#namesList = [name for name in dataDict.keys() if name.count('brd4') == 1]


#parentFolder = utils.formatFolder('%srose' % (projectFolder),True)
#bashFileName = '%sMCKINSEY_ENHANCER_BRD4_rose.sh' % (parentFolder)


#pipeline_dfci.callRose2Slurm(brd4_dataFile,macsEnrichedFolder,parentFolder,namesList,extraMap = [],inputFile='',tss=2500,stitch='',bashFileName =bashFileName,mask='')

##############################################################################

#dataDict = pipeline_dfci.loadDataTable(pol2_dataFile)

#namesList = dataDict.keys()
#namesList = [name for name in dataDict.keys() if name.count('pol2') == 1]


#parentFolder = utils.formatFolder('%srose' % (projectFolder),True)
#bashFileName = '%sMCKINSEY_ENHANCER_POL2_rose.sh' % (parentFolder)


#pipeline_dfci.callRose2(pol2_dataFile,macsEnrichedFolder,parentFolder,namesList,extraMap = [],inputFile='',tss=2500,stitch='',bashFileName =bashFileName,mask='')

###############################################################################

#dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)                                                                                                                                                                                                                          

#namesList = dataDict.keys()                                                                                                                                                                                                                                                    
#namesList = [name for name in dataDict.keys() if name.count('brd4') == 1]                                                                                                                                                                                                      


#parentFolder = utils.formatFolder('%srose_1000' % (projectFolder),True)                                                                                                                                                                                                            
#bashFileName = '%sMCKINSEY_ENHANCER_BRD4_rose_1000_stitch.sh' % (parentFolder)                                                                                                                                                                                                             


#pipeline_dfci.callRose2(brd4_dataFile,macsEnrichedFolder,parentFolder,namesList,extraMap = [],inputFile='',tss=2500,stitch=1000,bashFileName =bashFileName,mask='')                                                                                                              

##############################################################################                                                                                                                                                                                                  

#dataDict = pipeline_dfci.loadDataTable(pol2_dataFile)                                                                                                                                                                                                                          

#namesList = dataDict.keys()                                                                                                                                                                                                                                                    
#namesList = [name for name in dataDict.keys() if name.count('pol2') == 1]                                                                                                                                                                                                      


#parentFolder = utils.formatFolder('%srose_1000' % (projectFolder),True)                                                                                                                                                                                                             
#bashFileName = '%sMCKINSEY_ENHANCER_POL2_rose_1000_stitch.sh' % (parentFolder)                                                                                                                                                                                                             


#pipeline_dfci.callRose2(pol2_dataFile,macsEnrichedFolder,parentFolder,namesList,extraMap = [],inputFile='',tss=2500,stitch=1000,bashFileName =bashFileName,mask='')


#==========================================================================
#======================DYNAMIC COMPARISON==================================
#========================================================================== 


#set up the output
#name1 = 'brd4_veh'
#name2 = 'brd4_tgfb'

#analysisName = '%s_vs_%s_1000' % (name1,name2)

#outputFolder = '%sdynamic/%s/' % (projectFolder,analysisName)
#outputFolder = utils.formatFolder(outputFolder,True)

##get the rose folder
#roseParentFolder ='%srose_1000/' % (projectFolder)
#roseFolder1 = '%s%s_ROSE' % (roseParentFolder,name1)
#roseFolder2 = '%s%s_ROSE' % (roseParentFolder,name2)


#set up the bash file
#bashFileName = '%s%s_dynamic.sh' % (outputFolder,analysisName)
#bashFile = open(bashFileName,'w')

#bashFile.write('#!/usr/bin/bash\n')

#for now change into pipelinedir just to be safe
#bashFile.write('cd /storage/cylin/home/cl6/pipeline/\n')

#if you want to do w/ all do this
#this will run w/ default parameters
#see the documentation for available flags
#dynamicCmd = 'python /storage/cylin/home/cl6/pipeline/dynamicEnhancer.py -g %s #-d %s -n %s,%s -r %s,%s -o %s -e super -p' % (genome,brd4_dataFile,name1,name2,roseFolder1,roseFolder2,outputFolder)

#bashFile.write(dynamicCmd + '\n')
#print(dynamicCmd)

#bashFile.close()


#==========================================================================
#======================MAKING RN6 GFFS=====================================
#==========================================================================
#pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,species='RN6')

#==========================================================================
#======================MAPPING BAMS========================================
#========================================================================== 
'''

######## for pol2
#dataFile = pol2_dataFile
#gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_BODY_ALL_+300_+3000.gff','/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-300_+300.gff']

#pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList=['pol2_tgfb','pol2_veh','tgfb_control','veh_control'],extension=200)

#figureGFF = [['chr3','','',118106428,118548435,'','.','','Fgf7'],
#             ['chr13','','',111534031,111594413,'','.','Sertad4'],
#             ['chr1','','',200424698,200777694,'','-','','Fgfr2']]    

#figureGFF = [['chr10','','',69379587,69456654,'','.','','Ccl2']]


#figureGFFPath  ='%ssertad4_fgfr2_fgf7.gff' % (gffFolder)

#sertadGFF = '/storage/cylin/grail/projects/McKinsey/gff/Sertad4.gff'

sertadGFF = '/storage/cylin/grail/projects/McKinsey/gff/big_Sertad4.gff'

dataFile = brd4_dataFile

plotName = 'big_sertad4_brd4_plots'
outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
namesList = ['brd4_veh','brd4_tgfb']

#pipeline_dfci.callBatchPlot(dataFile,sertadGFF,plotName,outputFolder,namesList,uniform=True,bed ='')


dataFile = pol2_dataFile

plotName = 'big_sertad4_pol2_plots'
outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
namesList = ['pol2_veh','pol2_tgfb']

#pipeline_dfci.callBatchPlot(dataFile,sertadGFF,plotName,outputFolder,namesList,uniform=True,bed ='')



#dataFile = pol2_dataFile

#plotName = 'fgfr2_and_friends_pol2_ALL_TRACKS'
#outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
#namesList = []

#pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='')
'''
#------------


# harrison's plots (Revision 1)

#dataFile = brd4_dataFile
#dataFile = pol2_dataFile

#figure_gff = [['chr1','ACTA2','',252537613,252550394,'','-','','',''],
             # ['chr2','POSTN','',143656844,143688087,'','+','','','']
             
#             ]

#zoom out 50kb


#figure_gff = [['chr1','ACTA2','',252487613,252600394,'','-','','',''],
#              ['chr2','POSTN','',143606844,143738087,'','+','','','']

#            ]


#zoom out 100kb

#figure_gff = [['chr1','ACTA2','',252437613,252650394,'','-','','',''],
#             ['chr2','POSTN','',143556844,143788087,'','+','','','']

#             ]

#specific window -2kb +50kb
'''

figure_gff = [['chr1', 'ACTA2', '', 252477613, 252555394,'','-','',''],
              ['chr2', 'POSTN', '', 143646844, 143758087,'','+','','']

              ]




#figure_gff_path = '%srevision1/gff/brd4_plots.gff' % (projectFolder)
figure_gff_path = '%srevision1/gff/pol2_plots.gff' % (projectFolder)


utils.unParseTable(figure_gff,figure_gff_path, '\t')

#plot_name = 'BRD4_revisions'
#plot_name = 'BRD4_revisions_+50kb'
#plot_name = 'BRD4_revisions_+100kb'
#plot_name = 'BRD4_revisions_ACTA_-5kb_+60kb_&_postn_-10kb_+70kb'


#plot_name = 'Pol2_revision'
#plot_name = 'Pol2_revision_+50kb'
#plot_name = 'Pol2_revision_+100kb'
plot_name  = 'Pol2_revisions_ACTA_5kb_+60kb_POSTN_-10kb_+70kb'



#names_list = ['brd4_veh', 'brd4_tgfb']
names_list = ['pol2_veh', 'pol2_tgfb']

output_folder = '%srevision1/' % (projectFolder)

pipeline_dfci.callBatchPlot(dataFile,figure_gff_path, plot_name, output_folder, names_list, uniform=True, bed = '', plotType= 'MULTIPLE', extension=200, multiPage = False, debug=False, nameString = '', rpm=True, rxGenome = '', scaleFactorString = '') 
'''
#==========================================================================
#======================MAPPING ENRICHED====================================
#==========================================================================

#dataFile = pol2_dataFile
#gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000.gff']
#setName = 'NRVM_RNA_POL2'
#cellTypeList = ['pol2']
#enrichedFolder = macsEnrichedFolder
#namesList=['pol2_tgfb','pol2_veh']
#pipeline_dfci.mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)



#setList = [['pol2_tgfb'],['pol2_veh']]
#mappedEnrichedFile = '%smappedEnriched/RN6_TSS_ALL_-1000_+1000/RN6_TSS_ALL_-1000_+1000_NRVM_RNA_POL2.txt' % (projectFolder)
#output = '%stables/RN6_NRVM_POL2_TSS_ACTIVE_-1000_+1000.txt' % (projectFolder)
#pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)


'''

dataFile = brd4_dataFile
gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000.gff']
setName = 'BRD4'
cellTypeList = ['brd4']
enrichedFolder = macsEnrichedFolder
namesList=['brd4_tgfb','brd4_veh']
#pipeline_dfci.mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedEnrichedFolder,True,namesList,useBackground=True)
'''

####################################################################################
########################### make gffs ##############################################
####################################################################################
'''
mapped='%smappedEnriched/RN6_TSS_ALL_-1000_+1000/RN6_TSS_ALL_-1000_+1000_BRD4_conserved.txt' % (projectFolder)
table=utils.parseTable(mapped,'\t')
gff=[]
#chrom,name1,'',start,stop,'',sense,'',name1
for line in table:
    string=line[0]
    name=line[1]
    chrom=string.split('(')[0]
    sense=(string.split('(')[1]).split(')')[0]
    start_stop=string.split(':')[1]
    start=int(start_stop.split('-')[0])
    stop=int(start_stop.split('-')[1])
    new_line=[chrom,name,'',start,stop,'',sense,'',name]
    gff.append(new_line)

utils.unParseTable(gff,'/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000_BRD4_conserved.gff','\t')
'''

'''
mapped='%smappedEnriched/RN6_TSS_ALL_-1000_+1000/RN6_TSS_ALL_-1000_+1000_BRD4_tgfb.txt' % (projectFolder)
table=utils.parseTable(mapped,'\t')
gff=[]
#chrom,name1,'',start,stop,'',sense,'',name1
for line in table:
    string=line[0]
    name=line[1]
    chrom=string.split('(')[0]
    sense=(string.split('(')[1]).split(')')[0]
    start_stop=string.split(':')[1]
    start=int(start_stop.split('-')[0])
    stop=int(start_stop.split('-')[1])
    new_line=[chrom,name,'',start,stop,'',sense,'',name]
    gff.append(new_line)

utils.unParseTable(gff,'/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000_BRD4_tgfb.gff','\t')

'''
'''
mapped='%smappedEnriched/RN6_TSS_ALL_-1000_+1000/RN6_TSS_ALL_-1000_+1000_BRD4_veh.txt' % (projectFolder)
table=utils.parseTable(mapped,'\t')
gff=[]
#chrom,name1,'',start,stop,'',sense,'',name1
for line in table:
    string=line[0]
    name=line[1]
    chrom=string.split('(')[0]
    sense=(string.split('(')[1]).split(')')[0]
    start_stop=string.split(':')[1]
    start=int(start_stop.split('-')[0])
    stop=int(start_stop.split('-')[1])
    new_line=[chrom,name,'',start,stop,'',sense,'',name]
    gff.append(new_line)

utils.unParseTable(gff,'/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000_BRD4_veh.gff','\t')

'''


#==========================================================================
#======================MAP BAMS BATCH======================================
#==========================================================================


#dataFile = pol2_dataFile
#gffList =  ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000.gff']

#pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite=False,namesList = [],extension=200,rpm=False)



#==========================================================================
#======================SIGNAL TABLES=======================================
#========================================================================== 

#dataFile = pol2_dataFile
#gffFile = '/storage/cylin/grail/projects/McKinsey/gff/RN6_BODY_ALL_+300_+3000.gff'


#pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['pol2_tgfb','pol2_veh','tgfb_control','veh_control'],medianNorm=False,output ='/storage/cylin/grail/projects/McKinsey/signalTables/pol2_RN6_BODY_ALL_+300_+3000_Signal_Table.txt')



#dataFile = pol2_dataFile
#gffFile = '/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-300_+300.gff'


#pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['pol2_tgfb','pol2_veh','tgfb_control','veh_control'],medianNorm=False,output ='/storage/cylin/grail/projects/McKinsey/signalTables/pol2_RN6_TSS_ALL_-300_+300_Signal_Table.txt')


#dataFile = pol2_dataFile
#gffFile = '/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-1000_+1000.gff'

#pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['pol2_tgfb','pol2_veh','tgfb_control','veh_control'],medianNorm=False,output ='/storage/cylin/grail/projects/McKinsey/signalTables/pol2_RN6_TSS_ALL_-1000_+1000_Signal_Table.txt')

#==========================================================================
#==========================CALLING HISAT2==================================
#==========================================================================

#let's make a bash script to call hisat2 off of the SRA accession number

#start with the data file

#end w/ the bash script

#1 for each datafile we want to extract the SRA, the name, the genome etc...
#2 create a command and write it to the bash file
#3 add commands to convert sam to bam and to sort and index the bams

#def mapHisat(dataFile,bashFilePath,namesList=[],useSRA=True):

#     '''
#     maps using hisat2 if useSRA is flagged will try to extract an SRA ID from the fastq path and call directly
#     '''

#     #bashFilePath is where we will write the commands to

#     #first need to open
#    bashFile = open(bashFilePath,'w')

#    bashFile.write('#!/usr/bin/bash\n') #shebang line plus end of line characters
#    bashFile.write('\n\n\n')


#     #load the datasets
#    dataDict = pipeline_dfci.loadDataTable(dataFile)

#     #loading the datasets we want to process
#    if len(namesList) == 0:
#        namesList = dataDict.keys()

#     #want to write to the bam folder which we can deduce from the genomes of the datasets
#     #we want to make sure everyone has the same genome or else this is scary

#    genomeList = [dataDict[name]['genome'] for name in dataDict.keys()]
#    if len(utils.uniquify(genomeList)) > 1:
#        print('OH HECK NO YOU CANT ALIGN MULTIPLE GENOMES THAT WOULD BE STUPID')
#        sys.exit()

#    genome = string.upper(genomeList[0])


#    bamFolder = utils.formatFolder(dataDict[namesList[0]]['folder'])
#    hisatIndexDictionary = {'HG38':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Hisat2Index/hg38',
#                            'MM9': '/storage/cylin/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Hisat2Index/mm9',
#                            'RN6_ERCC': '/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Hisat2Index_ERCC/rn6_ercc',
#                            'RN6': '/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Hisat2Index/rn6'
#                           }

#    hisatIndex = hisatIndexDictionary[genome]


#     #first write a line to cd into the bam folder
#    bashFile.write('cd %s\n' %(bamFolder))

#     #now we can loop through the datasets
#    for name in namesList:
        
#        bashFileString =bashFilePath + name +'.sh'

#        bashFile = open(bashFileString,'w')

#        bashFile.write('#!/usr/bin/bash\n') #shebang line plus end of line characters
#        bashFile.write('\n\n\n')

#        bashFile.write('cd %s\n' %(bamFolder))        

#        uniqueID = dataDict[name]['uniqueID']
#        ts = time.time()
#        timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')
#        cmd = '#SBATCH --output=/storage/cylin/grail/slurm_out/hisat2_%s_%s' % (name,timestamp) + '_%j.out # Standard output and error log'
#        bashFile.write(cmd+'\n')
#        cmd = '#SBATCH -e /storage/cylin/grail/slurm_out/hisat2_%s_%s' % (name,timestamp) + '_%j.err # Standard output and error log'
#        bashFile.write(cmd+'\n')

#        cmd = 'pwd; hostname; date'
#        bashFile.write(cmd+'\n')
#        bashFile.write('\n\n\n')

#        bashFile.write('#===================\n')
#        bashFile.write('#PROCESSING %s\n' %(name))
#        bashFile.write('echo "processing %s"\n' % (name))

#        outputSam = '%s%s.%s.sam' % (bamFolder,uniqueID,genome)
#        outputBam = '%s%s.%s.bam' % (bamFolder,uniqueID,genome)
#        outputSortedBam = '%s%s.%s.sorted' % (bamFolder,uniqueID,genome)

#        if useSRA:
#            srrID = dataDict[name]['fastq'].split('/')[-2]
#            alignCmd = 'hisat2 -p 4 --no-unal -x %s --sra-acc %s -S %s' % (hisatIndex,srrID,outputSam)
#        else:
#            #check for paired end
#            fastqPath = dataDict[name]['fastq']
#            print('fastqPath')
#            print(fastqPath)
#            if fastqPath.count('::') == 1:
#                #this is paried end
#                [fastqPath_1,fastqPath_2] = fastqPath.split('::')
#                alignCmd = 'hisat2 -p 4 --no-unal -x %s -1 %s -2 %s -S %s' % (hisatIndex,fastqPath_1,fastqPath_2,outputSam)
#            else:
#                alignCmd = 'hisat2 -p 4 --no-unal -x %s -U %s -S %s' % (hisatIndex,fastqPath,outputSam)

#        bashFile.write(alignCmd)
#        bashFile.write('\n')

        #now convert the sam to a bam
#        generateBamCmd = 'samtools view -bS %s > %s' % (outputSam,outputBam)
#        bashFile.write(generateBamCmd)
#        bashFile.write('\n')

        #now we need to sort the bam
#        sortBamCmd = 'samtools sort %s %s' % (outputBam,outputSortedBam)
#        bashFile.write(sortBamCmd)
#        bashFile.write('\n')

        #now we need to index the bam
#        indexBamCmd = 'samtools index %s.bam' % (outputSortedBam)
#        bashFile.write(indexBamCmd)
#        bashFile.write('\n')

#         #now we need to delete the sam
#        deleteSamCmd = 'rm %s' % (outputSam)
#        bashFile.write(deleteSamCmd)
#        bashFile.write('\n')

#         #now we need to delete the unsorted bam
#        deleteBamCmd = 'rm %s' % (outputBam)
#        bashFile.write(deleteBamCmd)
#        bashFile.write('\n')
#        bashFile.write('\n\n\n')

#        bashFile.close()


#bashFilePath = '%srnaseq/' % (projectFolder)
#mapHisat(rna_dataFile,bashFilePath,namesList=[],useSRA=False)


#==========================================================================
#=======================RUNNING RNA-SEQ ANALYSIS===========================
#==========================================================================

'''
analysisName = 'mckinsey_rna'

gtfFile = '/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf'
cufflinksFolder = utils.formatFolder('%scufflinks' % (projectFolder),True)
bashFileName = '%s%s_cufflinks.sh' % (cufflinksFolder,analysisName)

groupList = [['MCKINSEY_RNA_DMSO_A','MCKINSEY_RNA_DMSO_B','MCKINSEY_RNA_DMSO_C'],['MCKINSEY_RNA_TGFB_DMSO_A','MCKINSEY_RNA_TGFB_DMSO_B','MCKINSEY_RNA_TGFB_DMSO_C'],['MCKINSEY_RNA_JQ1_A','MCKINSEY_RNA_JQ1_B','MCKINSEY_RNA_JQ1_C'],['MCKINSEY_RNA_TGFB_JQ1_A','MCKINSEY_RNA_TGFB_JQ1_B','MCKINSEY_RNA_TGFB_JQ1_C']]


print(groupList)
#pipeline_dfci.makeCuffTableSlurm(rna_dataFile,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)


# #flag useERCC to true

'''
#==========================================================================
#===================MAKING META GFFS OF TXN REGIONS========================
#==========================================================================

#def makeMetaGFFs(annotFile,gffFolder,genome,geneListFile =''):
#    '''
#    makes gffs of txn regions for meta genes
#    '''
        
#pipeline_dfci.makeMetaGFFs(annotFile,gffFolder,genome,geneListFile ='')

#==========================================================================
#===================MAPPING BAMS TO GFFS===================================
#==========================================================================
#dataFile = pol2_dataFile
#cellTypeList = ['pol2','tgfb','veh']
#gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-3000_+0.gff','/storage/cylin/grail/projects/McKinsey/gff/RN6_TTR_ALL_-0_+3000.gff'] #nBin = 60
#pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 60,overWrite =True,rpm=True,nameList = [])


#gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TXN_ALL_-0_+0.gff'] #nBin = 200
#pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200 ,overWrite =True,rpm=True,nameList = [])


###def mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,rpm=True,nameList = []):

#    '''
#    for each gff maps all of the data and writes to a specific folder named after the gff
#    can map either by cell type or by a specific name list
#    '''

#dataFile = pol2_dataFile
#dataFile = brd4_dataFile
#cellTypeList = ['brd4','tgfb','veh']

#gffList = ['/storage/cylin/grail/projects/McKinsey/gff/unique_genes_TSS_-1000_+1000.gff'] #nBin = 200
#pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200 ,overWrite =True,rpm=True,nameList = [])

#==================================================================
#==============Harrison's Edits for BRD4 meta======================
#==================================================================
'''
dataFile = brd4_dataFile
cellTypeList = ['brd4', 'tgfb', 'veh']
#gffList = ['/storage/cylin/grail/projects/McKinsey/gff/unique_genes_TSS_-1000_+1000.gff'] #nBin = 200
#gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-3000_+0.gff','/storage/cylin/grail/projects/McKinsey/gff/RN6_TTR_ALL_-0_+3000.gff', '/storage/cylin/grail/projects/McKinsey/gff/RN6_TXN_ALL_-0_+0.gff'] #nBin = 60 

gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TXN_ALL_-0_+0.gff'] #nBin = 200

pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200 ,overWrite =True,rpm=True,nameList = [])
'''

dataFile = pol2_dataFile
cellTypeList = ['pol2', 'tgfb', 'veh']
gffList = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TSS_ALL_-300_+300.gff']
#gfflist = ['/storage/cylin/grail/projects/McKinsey/gff/RN6_TXN_ALL_-0_+0.gff'] 
pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 1, overWrite = True, rpm = True, nameList = [])


















####################################################
################## GENIALIS TABLE ##################
####################################################

#outFilePath1 = '%s%s_brd4_genialis_annotation_table.txt' % (projectFolder,projectName)
#outFilePath2 = '%s%s_pol2_genialis_annotation_table.txt' % (projectFolder,projectName)
#outFilePath3 = '%s%s_rna_genialis_annotation_table.txt' % (projectFolder,projectName)


#pipeline_dfci.makeGenialisTable(rna_dataFile,outFilePath3,organism='',seqType='RNAseq',paired=False,collection='McKinsey',annotator='Rachel Hirsch',source='', tissue='',age='',genotype='',molecule='',libraryStrategy='',exPro='',libraryConst='',other1='',other2='')


####################################################
################## 0_STITCH GFF ##################
####################################################



#dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)
#namesList = ['brd4_tgfb','brd4_veh']

#allLoci = []
#for name in namesList:

#      collection = utils.importBoundRegion('/storage/cylin/grail/projects/McKinsey/macsEnriched/%s_peaks.bed' %(name),name)

#      allLoci += collection.getLoci()

# #do this for each one in the namesList
# #then make a giant collection

#giant_collection = utils.LocusCollection(allLoci,50)

#stitched_collection = giant_collection.stitchCollection()

#gff = utils.locusCollectionToGFF(stitched_collection)

#utils.unParseTable(gff,'/storage/cylin/grail/projects/McKinsey/gff/UNION_BRD4_STITCHED_-0_+0.gff','\t')

##################################################################
######################Bamplot example#############################
##################################################################


#testGFF = '/storage/cylin/grail/projects/McKinsey/gff/monika_test.gff'

#dataFile = brd4_dataFile

#plotName = 'Monika_test_plots'
#outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
#namesList = ['brd4_veh','brd4_tgfb']

#pipeline_dfci.callBatchPlot(dataFile,testGFF,plotName,outputFolder,namesList,uniform=True,bed ='',debug=False)

####################################################
### MAKE GFF TO FASTA FOR FIMO FROM DYNAMIC ENHANCER
####################################################

'''
gff = '/storage/cylin/grail/projects/McKinsey/dynamic/brd4_veh_vs_brd4_tgfb_all/RN6_brd4_veh_brd4_tgfb_merged_MERGED_REGIONS_-0_+0.gff'

fasta=utils.gffToFasta('RN6','/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Chromosomes/',gff,UCSC = True,useID=False)

utils.unParseTable(fasta,'/storage/cylin/grail/projects/McKinsey/FIMO/brd4_all.fasta','')

def makeMotifBackground(subpeakFasta,projectFolder,projectName):

'''
#    makes a 1st order markov background file for fimo
'''

    bgCmd = 'fasta-get-markov -m 1 < ' + subpeakFasta + '  > ' + projectFolder + projectName + '_bg.meme'
    bg_path = '%s%s_bg.meme' %(projectFolder,projectName)
    subprocess.call(bgCmd, shell=True)

    return bg_path

bg_path = makeMotifBackground('/storage/cylin/grail/projects/McKinsey/FIMO/brd4_all.fasta',projectFolder,'brd4_all')
'''

##################################################
### GENE LIST FOR GAINED AND LOST GENES FIMO
##################################################
'''
gl_tab_path='/storage/cylin/grail/projects/McKinsey/tables/gained_lost_genes.txt'
lost_tab_path='/storage/cylin/grail/projects/McKinsey/tables/lost_genes.txt'
gained_tab_path='/storage/cylin/grail/projects/McKinsey/tables/gained_genes.txt'
active_pol_path='/storage/cylin/grail/projects/McKinsey/tables/RN6_NRVM_POL2_TSS_ACTIVE_-1000_+1000.txt'
filtered_act_path='/storage/cylin/grail/projects/McKinsey/tables/filtered_active_table.txt'

gl_tab=utils.parseTable(gl_tab_path,'\t')
lost_tab=utils.parseTable(lost_tab_path,'\t')
gained_tab=utils.parseTable(gained_tab_path,'\t')
active_pol=utils.parseTable(active_pol_path,'\t')
filtered_act=utils.parseTable(filtered_act_path,'\t')


gl_act=[]
for line in gl_tab:
    tick=0
    tock=0
    gene=line[1]
    for act in active_pol:
        name=act[2]
        if gene==name:
            tick=1
            nm=act[1]
            new_line=[nm,name]
            gl_act.append(new_line)
    if tick == 0:
        for f_act in filtered_act:
            name=f_act[1]
            if gene==name:
                tock=1
                nm=f_act[0]
                new_line=[nm,name]
                gl_act.append(new_line)
    if tock == 0:
        print('%s has no match!') % (gene)

utils.unParseTable(gl_act,'/storage/cylin/grail/projects/McKinsey/tables/gl_active.txt','\t')


g_act=[]
for line in gained_tab:
    tick=0
    tock=0
    gene=line[1]
    for act in active_pol:
        name=act[2]
        if gene==name:
            tick=1
            nm=act[1]
            new_line=[nm,name]
            g_act.append(new_line)
    if tick == 0:
        for f_act in filtered_act:
            name=f_act[1]
            if gene==name:
                tock=1
                nm=f_act[0]
                new_line=[nm,name]
                g_act.append(new_line)
    if tock == 0:
        print('%s has no match!') % (gene)

utils.unParseTable(g_act,'/storage/cylin/grail/projects/McKinsey/tables/gained_active.txt','\t')


l_act=[]
for line in lost_tab:
    tick=0
    tock=0
    gene=line[1]
    for act in active_pol:
        name=act[2]
        if gene==name:
            tick=1
            nm=act[1]
            new_line=[nm,name]
            l_act.append(new_line)
    if tick == 0:
        for f_act in filtered_act:
            name=f_act[1]
            if gene==name:
                tock=1
                nm=f_act[0]
                new_line=[nm,name]
                l_act.append(new_line)
    if tock == 0:
        print('%s has no match!') % (gene)

utils.unParseTable(l_act,'/storage/cylin/grail/projects/McKinsey/tables/lost_active.txt','\t')
'''


'''
refseqTable = utils.parseTable(annotFile,'\t')

nm_list=[]
AGL=[]
for gene in gl_tab:
    print(gene[1])
    nm_id=[]
    for line in refseqTable[1:]:
        if line[12]==gene[1]:
            name=line[1]
            nm_id.append(name)
    print(nm_id)
    if len(nm_id)>0:
        new_line=[','.join(nm_id),gene[1]]
        nm_list.extend(nm_id)
        AGL.append(new_line)

outpath='/storage/cylin/grail/projects/McKinsey/tables/gl_activeGenesList.txt'
utils.unParseTable(AGL,outpath,'\t')

nm_list=[]
AGL=[]
for gene in gained_tab:
    print(gene[1])
    nm_id=[]
    for line in refseqTable[1:]:
        if line[12]==gene[1]:
            name=line[1]
            nm_id.append(name)
    print(nm_id)
    if len(nm_id)>0:
        new_line=[','.join(nm_id),gene[1]]
        nm_list.extend(nm_id)
        AGL.append(new_line)

outpath='/storage/cylin/grail/projects/McKinsey/tables/gained_activeGenesList.txt'
utils.unParseTable(AGL,outpath,'\t')

nm_list=[]
AGL=[]
for gene in lost_tab:
    print(gene[1])
    nm_id=[]
    for line in refseqTable[1:]:
        if line[12]==gene[1]:
            name=line[1]
            nm_id.append(name)
    print(nm_id)
    if len(nm_id)>0:
        new_line=[','.join(nm_id),gene[1]]
        nm_list.extend(nm_id)
        AGL.append(new_line)

outpath='/storage/cylin/grail/projects/McKinsey/tables/lost_activeGenesList.txt'
utils.unParseTable(AGL,outpath,'\t')
'''

'''
AGL=utils.parseTable('/storage/cylin/grail/projects/McKinsey/tables/gl_activeGenesList.txt','\t')

nm_list=[]
for line in AGL:
    nm_ids=line[0].split(',')
    nm_list.extend(nm_ids)

print(nm_list[0:100])

unique_nms=utils.uniquify(nm_list)


nm_genes_list=[]
for nm in unique_nms:
    genes=[]
    for line in AGL:
        nm_ids=line[0].split(',')
        if nm in nm_ids:
            genes.append(line[1])
    new_line=[nm,','.join(genes)]
    nm_genes_list.append(new_line)


outpath='/storage/cylin/grail/projects/McKinsey/tables/gl_NM_Genes.txt'
utils.unParseTable(nm_genes_list,outpath,'\t')



AGL=utils.parseTable('/storage/cylin/grail/projects/McKinsey/tables/lost_activeGenesList.txt','\t')

nm_list=[]
for line in AGL:
    nm_ids=line[0].split(',')
    nm_list.extend(nm_ids)

print(nm_list[0:100])

unique_nms=utils.uniquify(nm_list)


nm_genes_list=[]
for nm in unique_nms:
    genes=[]
    for line in AGL:
        nm_ids=line[0].split(',')
        if nm in nm_ids:
            genes.append(line[1])
    new_line=[nm,','.join(genes)]
    nm_genes_list.append(new_line)


outpath='/storage/cylin/grail/projects/McKinsey/tables/lost_NM_Genes.txt'
utils.unParseTable(nm_genes_list,outpath,'\t')


AGL=utils.parseTable('/storage/cylin/grail/projects/McKinsey/tables/gained_activeGenesList.txt','\t')

nm_list=[]
for line in AGL:
    nm_ids=line[0].split(',')
    nm_list.extend(nm_ids)

print(nm_list[0:100])

unique_nms=utils.uniquify(nm_list)


nm_genes_list=[]
for nm in unique_nms:
    genes=[]
    for line in AGL:
        nm_ids=line[0].split(',')
        if nm in nm_ids:
            genes.append(line[1])
    new_line=[nm,','.join(genes)]
    nm_genes_list.append(new_line)


outpath='/storage/cylin/grail/projects/McKinsey/tables/gained_NM_Genes.txt'
utils.unParseTable(nm_genes_list,outpath,'\t')

'''
#############################
# Enhancer region heatmaps ##
#############################

'''
gff_path = '/storage/cylin/grail/projects/McKinsey/gff/RN6_brd4_veh_brd4_tgfb_merged_MERGED_REGIONS_-0_+0.gff'

gff = utils.parseTable(gff_path,'\t')

new_gff = []

for line in gff:
    chrom=line[0]
    name = line[1]
    spacer = ''
    start = int(line[3])
    stop = int(line[4])
    sense = '.'
    mid_point = int(round((start+stop)/2))
    if mid_point >= 10000:
        new_start = mid_point - 10000
        new_stop = mid_point + 10000
    if mid_point < 10000:
        new_start = 0
        new_stop = 20000

    new_line = [chrom,name,spacer,new_start,new_stop,spacer,sense,spacer,name]
    new_gff.append(new_line)

utils.unParseTable(new_gff,'/storage/cylin/grail/projects/McKinsey/gff/All_enhancers_20kb_regions_no_neg.gff','\t')
'''

##########################
# means means means ######
##########################
'''
def getMeans(table):
    means=[]
    for i in range(2,202):
        print(i)
        sums=0
        for line in table[1:]:
            print(i)
            sums=sums+float(line[i])
        means.append(sums/7813)
    return(means)

veh = utils.parseTable('/storage/cylin/grail/projects/McKinsey/mappedFolder_n200/All_enhancers_20kb_regions/All_enhancers_20kb_regions_brd4_veh.gff','\t')
veh_control = utils.parseTable('/storage/cylin/grail/projects/McKinsey/mappedFolder_n200/All_enhancers_20kb_regions/All_enhancers_20kb_regions_veh_control.gff','\t')

tgfb = utils.parseTable('/storage/cylin/grail/projects/McKinsey/mappedFolder_n200/All_enhancers_20kb_regions/All_enhancers_20kb_regions_brd4_tgfb.gff','\t')
tgfb_control = utils.parseTable('/storage/cylin/grail/projects/McKinsey/mappedFolder_n200/All_enhancers_20kb_regions/All_enhancers_20kb_regions_tgfb_control.gff','\t')

sum_veh = getMeans(veh)
sum_veh_control = getMeans(veh_control)
sum_tgfb = getMeans(tgfb)
sum_tgfb_control = getMeans(tgfb_control)

utils.unParseTable(sum_veh,'/storage/cylin/grail/projects/McKinsey/tables/sum_veh.txt','\t')
utils.unParseTable(sum_veh_control,'/storage/cylin/grail/projects/McKinsey/tables/sum_veh_control.txt','\t')
utils.unParseTable(sum_tgfb,'/storage/cylin/grail/projects/McKinsey/tables/sum_tgfb.txt','\t')
utils.unParseTable(sum_tgfb_control,'/storage/cylin/grail/projects/McKinsey/tables/sum_tgfb_control.txt','\t')
'''
