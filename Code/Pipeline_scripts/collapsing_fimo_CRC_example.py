#!/usr/bin/env python

######################
#
# Core Regulatory Circuits
# Young and Bradner Labs
# Version 1.0
# 140724
#
######################

######################
# Dependencies
######################



import os
import sys
import string

pipeline_dir = '/storage/cylin/bin/pipeline/' 

sys.path.append(pipeline_dir)
print(pipeline_dir)
import utils




import numpy
import scipy
import scipy.stats

import subprocess
import os


from random import randrange
from collections import defaultdict

import pickle



#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

class Genome:
    __chrDict = dict()
    __featureDict = dict()

    #at its core, the genome has attributes of a build name, a fasta directory and an annotation file
    def __init__(self,name,genome_directory,annot_file):
        self._name = name
        self._directory = genome_directory
        self._annot = annot_file

    def name(self):
        return self._name

    def directory(self):
        return self._directory

    def annot(self):
        return self._annot

    def addFeature(self,feature,path):

        if feature in self.__featureDict:
            print('WARNING OVERRIDING %s PATH WITH %s' % (feature,path))
        self.__featureDict[feature] = path

    def returnFeature(self,feature):

        #tries to load the selected feature from the feature dictionary

        if feature not in self.__featureDict:
            print('ERROR: GENOME %s DOES NOT HAVE FEATURE %s' % (self.name(),feature))
            sys.exit()
        else:
            return self.__featureDict[feature]

    def hasFeature(self,feature):

        if feature in self.__featureDict:
            return True
        else:
            return False

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


def loadGenome(genome_build):

    '''
    loads annotation for a genome into a genome object
    '''

    #this nested dictionary has all of the useful information and likely will have to be
    #edited so it can be configured any time
    genome_build = string.upper(genome_build)
    genomeDict = {
        'HG19':{'annot_file':'%sannotation/hg19_refseq.ucsc' % (pipeline_dir),
                'genome_directory':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/',
                'tf_file':'%scrc/annotation/TFlist_NMid_hg19.txt' % (pipeline_dir),
                'mask_file':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Annotation/Masks/hg19_encode_blacklist.bed',
                'motif_convert':'%scrc/annotation/MotifDictionary.txt' % (pipeline_dir),
                'motif_database':'%scrc/annotation/VertebratePWMs.txt' % (pipeline_dir),
                },
        'RN6':{'annot_file':'%sannotation/rn6_refseq.ucsc' % (pipeline_dir),
                'genome_directory':'/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Chromosomes/',
                'tf_file':'%scrc/annotation/TFlist_NMid_rn6.txt' % (pipeline_dir),
                'motif_convert':'%scrc/annotation/MotifDictionary.txt' % (pipeline_dir),
                'motif_database':'%scrc/annotation/VertebratePWMs.txt' % (pipeline_dir),
                }

        }

    if genome_build not in genomeDict:
        print('ERROR: UNSUPPORTED GENOME BUILD %s. EXITING NOW' % (genome_build))
        sys.exit()

    #now attempt to load the genome
    
    genome = Genome(genome_build,genomeDict[genome_build]['genome_directory'],genomeDict[genome_build]['annot_file'])

    #adding additional optional features
    genome.addFeature('tf_file',genomeDict[genome_build]['tf_file'])
    if genome_build == 'HG19':
        genome.addFeature('mask',genomeDict[genome_build]['mask_file'])
    genome.addFeature('motif_convert',genomeDict[genome_build]['motif_convert'])
    genome.addFeature('motif_database',genomeDict[genome_build]['motif_database'])


    return genome



######################
# Functions
######################


def collapseFimo(fimo_output,gene_to_enhancer_dict,candidate_tf_list,output_folder,analysis_name,motifConvertFile):

    '''
    collapses motifs from fimo
    for each source node (TF) and each target node (gene enhancer regions), collapse motif instances
    then spit out a ginormous set of beds and a single crazy collapsed bed
    '''

    #first build up the motif name conversion database

    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = defaultdict(list)
    # The reverse of the other dict, from motif name to gene name
    # a motif can go to multiple genes
    for line in motifDatabase:
        motifDatabaseDict[line[0]].append(line[1])



    #make the folder to store motif beds
    utils.formatFolder('%smotif_beds/' % (output_folder),True)

    edgeDict = {}
    #first layer are source nodes
    for tf in candidate_tf_list:
        edgeDict[tf] = defaultdict(list) #next layer are target nodes which are derived from the fimo output


    fimoTable = utils.parseTable(fimo_output,'\t')
    print(fimo_output)

    #fimo sometimes puts the region in either the first or second column
    fimo_line = fimoTable[1]
    if fimo_line[1].count('|') >0:
        region_index = 1
    else:
        region_index = 2
    print('USING COLUMN %s OF FIMO OUTPUT FOR REGION' % (region_index))

    for line in fimoTable[1:]:
        source_tfs = motifDatabaseDict[line[0]]   #motifId
        for source in source_tfs:
            if candidate_tf_list.count(source) == 0:
                continue
            region = line[region_index].split('|')

            target = region[0]
            if region_index == 2:
                target_locus = utils.Locus(region[1],int(region[2]) + int(line[3]), int(region[2]) + int(line[4]),'.')
            else:
                target_locus = utils.Locus(region[1],int(region[2]) + int(line[2]), int(region[2]) + int(line[3]),'.')
            #what's missing here is the enhancer id of the target locus
            try:
                edgeDict[source][target].append(target_locus)
            except KeyError:
                print('this motif is not in the network')
                print(line)
                sys.exit()



    #now we actually want to collapse this down in a meaningful way
    #overlapping motifs count as a single binding site. This way a TF with tons of motifs
    #that finds the same site over and over again doesn't get over counted
    all_bed = []
    all_bed_path = '%s%s_all_motifs.bed' % (output_folder,analysis_name)
    for tf in candidate_tf_list:
        print(tf)
        target_nodes = edgeDict[tf].keys()
        bed_header = ['track name="%s" description="%s motifs in %s"' % (tf,tf,analysis_name)]
        all_bed.append(bed_header)
        target_bed = [bed_header]
        target_bed_path = '%smotif_beds/%s_motifs.bed' % (output_folder,tf)
        for target in target_nodes:
            edgeCollection = utils.LocusCollection(edgeDict[tf][target],50)
            edgeCollection = edgeCollection.stitchCollection()
            edgeLoci = edgeCollection.getLoci()
            edgeDict[tf][target] = edgeLoci
            for locus in edgeLoci:
                bed_line = [locus.chr(),locus.start(),locus.end(),target,'','+']
                target_bed.append(bed_line)
                all_bed.append(bed_line)

        utils.unParseTable(target_bed,target_bed_path,'\t')

    #now the loci are all stitched up
    utils.unParseTable(all_bed,all_bed_path,'\t')
    return edgeDict


def geneToEnhancerDict(genome, enhancer_file, activity_path):
    '''
    Assign each Super-Enhancer to the closest active TSS to its center
    Return a dictionary keyed by TF that points to a list of loci
    '''
    print('Identifying enhancers and target genes from %s' %(enhancer_file))
    #should this do gene assignment????
    #for now assume gene assignment has been done
    #can later toggle to do gene assignment

    #first load the TF lists

    tf_table = utils.parseTable(genome.returnFeature('tf_file'), '\t')

    motif_table = utils.parseTable(genome.returnFeature('motif_convert'),'\t')
    #this gives all tfs that have a motif
    motif_tfs = utils.uniquify([line[1] for line in motif_table])
    #intersect w/ the activity table
    if len(activity_path) > 0:
        activity_table = utils.parseTable(activity_path,'\t')
        #figure out the right column for actual gene names (basically not NM or NR and not a numeral)
        for i in range(len(activity_table[0])):
            try:
                foo = int(activity_table[0][i])
            except ValueError:
                continue
        if activity_table[0][i][0:2] != 'NM' and activity_table[0][i][0:2] != 'NR': #assumes refseq
            gene_col = i
        print('using column %s of %s gene activity table for common names' % (gene_col + 1, activity_path))


        active_gene_list = [string.upper(line[gene_col]) for line in activity_table]
        tf_list_refseq = [line[0] for line in tf_table if active_gene_list.count(line[1]) > 0 and motif_tfs.count(line[1]) > 0]
        tf_list_name = utils.uniquify([line[1] for line in tf_table if active_gene_list.count(line[1]) > 0 and motif_tfs.count(line[1]) > 0])
    else:
        tf_list_refseq = [line[0] for line in tf_table if motif_tfs.count(line[1]) >0]
        tf_list_name = [line[1] for line in tf_table if motif_tfs.count(line[1]) >0]

    print('Identified %s TFs from %s that have motifs' % (len(tf_list_name),genome.returnFeature('tf_file')))

    #keyed by gene with loci objects in the list
    gene_to_enhancer_dict = defaultdict(list)
    enhancer_to_gene_dict = defaultdict(list)


    #assuming id,chrom,start,stop w/ gene names in the last 3 columns per standard ROSE output
    enhancer_table = utils.parseTable(enhancer_file,'\t')
    print('Analyzing %s cis-regulatory regions' % (len(enhancer_table)))

    #now let's make the enhancer table by region and then by gene
    enhancerTable = [['ENHANCER_ID','CHROM','START','STOP','GENE_LIST']]
    enhancerTFTable = [['ENHANCER_ID','CHROM','START','STOP','GENE_LIST']]
    geneTable = [['GENE','TF','CHROM','START','STOP','ENHANCER_ID']]
    geneTFTable = [['GENE','CHROM','START','STOP','ENHANCER_ID']]
    geneSummaryTable = [['GENE','TF','ENHANCER_LIST']]

    #will need to track which ones are TFs
    candidate_tf_list = []
    #find the columns for gene assignment
    header = enhancer_table[0]
    header_length = len(enhancer_table[0])
    closest_index = header.index('CLOSEST_GENE')
    proximal_index = header.index('PROXIMAL_GENES')
    overlap_index = header.index('OVERLAP_GENES')
    for line in enhancer_table[1:]:
        if len(line) != header_length: #don't bother trying to figure out lines w/o target genes
            continue
        enhancer_locus = utils.Locus(line[1],line[2],line[3],'.',line[0])
        closest_gene_list = line[closest_index].split(',')
        proximal_gene_list = line[proximal_index].split(',')
        overlap_gene_list = line[overlap_index].split(',')
        all_gene_list = closest_gene_list + proximal_gene_list + overlap_gene_list
        all_gene_list = [string.upper(gene) for gene in all_gene_list]

        #print(all_gene_list)

        #print(activity_path)
        #print(active_gene_list)
        #gets a unique list of all tfs

        if len(activity_path) > 0:
            all_gene_list = utils.uniquify([gene for gene in all_gene_list if active_gene_list.count(gene) > 0])
        else:
            all_gene_list = utils.uniquify(all_gene_list)
        candidate_gene_list = utils.uniquify([gene for gene in all_gene_list if tf_list_name.count(gene) > 0])
        if len(all_gene_list) > 0:
            for gene in all_gene_list:

                gene_to_enhancer_dict[gene].append(enhancer_locus)
                enhancer_to_gene_dict[enhancer_locus].append(gene)
            newLine = line[0:4] + [','.join(all_gene_list)]
        else:
            newLine = line[0:4] + ['']
        enhancerTable.append(newLine)

        if len(candidate_gene_list) > 0:
            tfLine = line[0:4] + [','.join(candidate_gene_list)]
            enhancerTFTable.append(tfLine)

    #now iterate through each gene and list the enhancers
    gene_list = gene_to_enhancer_dict.keys()
    print(gene_list)
    gene_list.sort()
    for gene in gene_list:
        if tf_list_name.count(gene) > 0:
            tf_status = 1
            candidate_tf_list.append(gene)
        else:
            tf_status = 0
        enhancer_loci = gene_to_enhancer_dict[gene]
        enhancerString =','.join([enhancer.ID() for enhancer in enhancer_loci])
        geneSummaryTable.append([gene,tf_status,enhancerString])
        for enhancer in enhancer_loci:
            newLine = [gene,tf_status,enhancer.chr(),enhancer.start(),enhancer.end(),enhancer.ID()]
            geneTable.append(newLine)
            if tf_status == 1:
                newLine = [gene,enhancer.chr(),enhancer.start(),enhancer.end(),enhancer.ID()]
                geneTFTable.append(newLine)

    #return geneTable,geneTFTable,enhancerTable,enhancerTFTable,geneSummaryTable,candidate_tf_list,gene_to_enhancer_dict
    return candidate_tf_list,gene_to_enhancer_dict


######################
#
# Main Method
#
######################

def main():

    genome_build = 'RN6'

    #return genome type object
    #this is what needs to be passed as an argument to
    print('loaded genome')
    genome = loadGenome(genome_build)


    fimo_output = '/storage/cylin/grail/projects/McKinsey/FIMO/brd4_all_fimo.txt'

    #don't need to put motif_beds in the output_folder, adds it automatically
    output_folder = '/storage/cylin/grail/projects/McKinsey/FIMO/'
    
    #i would change this to SRF so we don't overwrite
    analysis_name = 'McKinsey_brd4_all'

    enhancer_file = '/storage/cylin/grail/projects/McKinsey/tables/gained_enhancers.txt'
    activity_path = '/storage/cylin/grail/projects/McKinsey/tables/gl_NM_Genes.txt'

    motifConvertFile = '/storage/cylin/bin/pipeline/crc/annotation/MotifDictionary.txt'

    print('loaded gene to enhancer dict')    
    candidate_tf_list,gene_to_enhancer_dict = geneToEnhancerDict(genome, enhancer_file, activity_path)

    #overwrite the candidate_tf_list since we only care about SRF
    candidate_tf_list = []

    
    collapseFimo(fimo_output,gene_to_enhancer_dict,candidate_tf_list,output_folder,analysis_name,motifConvertFile)



    print('done running! :D')



if __name__ == '__main__':
    main()
