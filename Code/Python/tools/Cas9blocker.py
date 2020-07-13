#This is the Cas9 Blocker! 

#Determine sequence match
#Locate pam site frame
#Alter the overlapping PAM sequence without changing the codon-AA

import re                             
import collections
from collections import defaultdict
import os
import sys
sys.path.append('/home/harrison/Downloads/pipeline-master')
import utils
import math
import pandas as pd
from openpyxl import load_workbook 
import csv
#from Cas9blocker3 import scrubber

#parses the tab delimited table
#
#
table=utils.parseTable('/home/harrison/PPscripts/HS_tools/codon_hg19.txt', '\t')




#make the first row the key and the next 2 rows the values for that key
#
#
hg19_bias = defaultdict()
for line in table[1:]:
    hg19_bias[line[0]]=[line[1],float(line[2])]



#This checks if the input is a text files, opens and reads it
#
#
if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2]):


    with open(sys.argv[1]) as f, open(sys.argv[2]) as g:

                        
        dna1 = f.read().rstrip()                                                
        gRNA = g.read().rstrip()


#This checks whether only one file is given and selects the cells in a table 
#to be used as input
#
#
elif sys.argv[1][0] == 't':

    
    xl = pd.ExcelFile('HelicaseCandidatesgRNAssequences3.xlsx')
    df = xl.parse('Sheet1')

    
    for index, row in df.iterrows():
        #scrubber(row['sequence1'], row['transcript_id'])
        
                
        gRNA = row['sequence1']
        dna1 = row['transcript_id']
        
        if gRNA > 1:
            print gRNA

        quit()
    
    

    #Here we check if the Accession number was given for the sequence
    #Search through a refGene file for coordinates to the gene
    #Find the sequence that pertains to the gene transcript or CDS
    #
    #
elif sys.argv[1][0] == "N":         # <--- might just be able to remove this block
    gRNA = sys.argv[2]
    dna1 = sys.argv[1]

else:                                                                                
    dna1 = sys.argv[1]
    gRNA = sys.argv[2]  
    pass

    
ref_gene_table = utils.parseTable('/home/harrison/PPscripts/hg38/refGene.txt', '\t')

for line in ref_gene_table[:]:
    if line[1] == dna1:
    
        if line[3] == '+':
            print line[1]
            line[6] = int(line[6])
            CDS_start = line[6] + 1
            line[7] = int(line[7])
            CDS_stop = line[7] + 1
            
            print ('( '+ line[2] + ', ' + line[3]), ', ', CDS_start, ', ', CDS_stop, ' )' 
            print '^ represents the CDS region'
            dna1 = utils.fetchSeq('/home/harrison/PPscripts/hg38/Chromosomes/', line[2], CDS_start, CDS_stop)
            dna1 = dna1.upper()
            
            
        elif line[3] == '-' : #and line[2] == 'chr6' :
             
            
            print line[1]
            line[6] = int(line[6])
            CDS_start = line[6] + 1
            line[7] = int(line[7])
            CDS_stop = line[7] + 1

            print ('( '+ line[2] + ', ' + line[3]), ', ', CDS_start, ', ', CDS_stop, ' )' 
            print '^ represents the CDS region'

            line[5] = int(line[5])
            
            #quit() 
            dna1 = utils.fetchSeq('/home/harrison/PPscripts/hg38/Chromosomes/', line[2], CDS_start, CDS_stop)
            dna1 = dna1.upper()
            dna1 = utils.revComp(dna1)
                

    

#If the input was not a file than continue with the code
#
#

            
#def scrubber (dna1,gRNA):


#This will replace all 'T's with 'U's if DNA was inputted otherwise    
#
#
dna1_Ts = dna1
gRNA_Ts = gRNA
dna1 = dna1.upper()
gRNA = gRNA.upper()
dna1 = dna1.replace('T', 'U')                                                
gRNA = gRNA.replace('T', 'U')
 
#create reverse complement function
#
#
def reverse_complement(dna):
                complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A',}
                return ''.join([complement[base] for base in dna[::-1]])

#this gives you a reverse version of the target and guide
#
#
#rev_gRNA = gRNA[::-1]

revcom_dna1 = reverse_complement(dna1)


#Checks if the match is unique in the original sequence
#
#
if dna1.count(gRNA) == 1 and revcom_dna1.count(gRNA) == 0:
    
    print 'This is the target:\n'+ dna1_Ts + '(Coding region)'
    print 'This is the guide:\n' + gRNA_Ts  

    #Finds the location of the match
    #
    #
    for a in re.finditer(gRNA, dna1):    

        #This makes sure the adjacent 3 bases make a PAM site
        #
        #
        if dna1[a.end() + 1] == 'G' and dna1[a.end() + 2] == 'G':
            
            #This makes the PAM site a new PAM variable
            #
            #
            PAM = dna1[a.end()], dna1[a.end() + 1], dna1[a.end() + 2]

            PAM = ''.join(PAM)

            print 'The guide matches the sequence once at: '
            print (a.start(), a.end())

            #This print's out the same codon/PAM in DNA if DNA was given or leaves it as RNA if not.
            #
            #
            if sys.argv[1].find('U') == -1:
                t = PAM.replace('U', 'T')
            elif sys.argv[1].find('T') == -1:
                t = PAM    
            print t + ' is the PAM sequence'
                
            #remainder indirectly gives postion for base 1 in PAM telling us the frame its in
            #
            #
            if (a.end()%3) == 0:
                print'The PAM is in a 123 frame'
                
                #if 123 frame, PAM is the same as the codon that needs to be changed
                #
                #
                print t + ' Is the overlapping codon'
                

                #Incase the overlapping codon is special and does not have replacements
                #The code will alter what's known as the core sequence
                #Changing the core should have a similar effect as altering the PAM
                #
                #
                if hg19_bias[PAM][0] == 'W' or hg19_bias[PAM][0] == 'M':
                    print 'Cannot alter Tryptophan or Methionine due to only one codon'
                    print 'Now attempting to alter a codon in the core seqeunce (4 bases upstream from PAM site)'
                    Core_Seq = dna1[a.end() - 6] + dna1[a.end() - 5] + dna1[a.end() - 4]
                                            
                    Core_Seq = ''.join(Core_Seq)
                    
                    print Core_Seq + ' produces ' + hg19_bias[Core_Seq][0]
                    

                    if hg19_bias[Core_Seq][0] == 'W':

                        print 'Cannot alter a portion of the core sequence codon due to Tryptophan'

                        Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                        
                        i = Core_Seq_2

                        if hg19_bias[Core_Seq_2][0] == 'W':



                            #for finding the swap easily in the sequence
                            #
                            #
                            i = i.lower()
                         
                            dna2 = list(dna1)
                        
                            #The first codon near the core sequence is changed but not the second
                            #
                            #
                            #dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i
                            dna2 = ''.join(dna2)
                            print 'cannot alter adjacent codon because it also produces tryptophan'
                            print 'Sorry, but the target sequence cannot be adjusted to block cas9 at this time.'
                            
                            
                            quit()

                        else:

                            print 'Now swapping adjacent codon for efficacy'
                            #So it is easier to read in the seqeunce
                            

                            i = i.lower()
                             

                            #this makes a copy to work with to keep the original for whatever reason
                            #
                            # 
                            dna2 = list(dna1)
                            



                            #replaces old codon overlapping PAM with new similar codon
                            #
                            #
                            #dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i


                            #Now the neighboring codon will be altered since for efficacy 
                            #and since its not a special codon
                            #
                            #
                            Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                            Core_Seq_2 = ''.join(Core_Seq_2)      
                            

                            print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2][0]
                            
                            g = Core_Seq_2
                            
                            
                            codon_list_2 = []
                            for h in hg19_bias:
                                if hg19_bias[h][0]==hg19_bias[Core_Seq_2][0]:
                                    codon_list_2.append(h)


                            for j in codon_list_2:
                                if j == g:
                                    codon_list_2.remove(j) 

                                    k = hg19_bias[Core_Seq_2]
                                    del hg19_bias[Core_Seq_2]   


                                    freq_diff_2 = 0
                                    freq_diff_list_2 = [] 
                                
                                    for l in codon_list_2:
                                        

                                        freq_diff_2 = k[1] - hg19_bias[l][1]
                                        freq_diff_2 = math.fabs(freq_diff_2)
                                        freq_diff_list_2.append(freq_diff_2)
                                                  
                                        
                                    freq_diff_list_2.sort()
                                    z = freq_diff_list_2[0] 
                                    for l in codon_list_2:
                                        if math.fabs(k[1] - hg19_bias[l][1]) == z: # and l != 'AGG':
                           

                                            if sys.argv[1].find('U') == -1:
                                                l = l.replace('U', 'T')
                                            elif sys.argv[1].find('T') == -1:
                                                l = l

                                            print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                            
                                            #So it is easier to read in the seqeunce
                                            l = l.lower()
                                             
                                            
                                            
                                            # #replaces old codon overlapping PAM with new similar codon
                                            dna2[a.end() - 9], dna2[a.end() - 8], dna2[a.end() - 7] = l

                                            #list becomes a string again for legibility 
                                            dna2 = ''.join(dna2)
                            
                                            # print out the new DNA string
                                            if sys.argv[1].find('U') == -1:
                                                dna2 = dna2.replace('U', 'T')
                                                print dna2
                                                quit()
                                            else:
                                                print dna2   
                                                quit()

                        
                    
                    t = Core_Seq
                    
                    #make a list of other codons that have the same amino acid
                    #
                    #
                    codon_list = []
                    for i in hg19_bias:
                        if hg19_bias[i][0]==hg19_bias[Core_Seq][0]:
                            codon_list.append(i)
        
                    
                    hg19_bias_2 = hg19_bias

                    #finds a different codon with the highest codon bias       
                    #
                    #
                    for b in codon_list:
                        if b == t:
                            codon_list.remove(b)
                                           
                            f = hg19_bias[Core_Seq]
                            del hg19_bias[Core_Seq]
                            

                            freq_diff = 0
                            freq_diff_list = [] 
                        
                            #Here we created a list that holds frequencie differences between
                            #a new codon and the original
                            #the list is sorted for each codon and the codon with the smallest difference
                            #is selected to replace the original
                            #
                            #
                            for i in codon_list:
                                

                                freq_diff = f[1] - hg19_bias[i][1]
                                freq_diff = math.fabs(freq_diff)
                                freq_diff_list.append(freq_diff)
                                          
                            #sort the list from smallest to greatest
                            #the first on the list has a frequency closest to that of original
                            #
                            #    
                            freq_diff_list.sort()
                            x = freq_diff_list[0] 
                            for i in codon_list:
                                

                                #Select the first codon in the list
                                #
                                #
                                if math.fabs(f[1] - hg19_bias[i][1]) == x:
                                   

                                    if sys.argv[1].find('U') == -1:
                                        i = i.replace('U', 'T')
                                    elif sys.argv[1].find('T') == -1:
                                        i = i

                                    print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                    
                                    i = i.lower()

                                    dna2 = list(dna1)
                                    
                                    dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i
  
                                        

                                    #So now 1 codon near the core sequence has been changed
                                    #This is good but changes in all the codons overlapping the core is more beneficial
                                    #Now we will change the neighboring codon unless its a special codon that cannot be
                                    #
                                    #
                                    Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                                    if Core_Seq_2 == Core_Seq:

                                        print 'The adjecent codon is the same and will also be changed for efficacy'
                                        dna2[a.end() - 9], dna2[a.end() - 8], dna2[a.end() - 7] = i
                                        
                                        dna2 = ''.join(dna2)
                                        if sys.argv[1].find('U') == -1:
                                            dna2 = dna2.replace('U', 'T')

                                            print dna2
                                            quit()

                                    if hg19_bias_2[Core_Seq_2][0] == 'W':



                                        print 'cannot alter adjacent codon because it produces tryptophan'
                                        
                                        dna2 = ''.join(dna2)
                                        if sys.argv[1].find('U') == -1:
                                            dna2 = dna2.replace('U', 'T')

                                        print dna2
                                        quit()

                                    else:

                                        print 'Also swapping adjacent codon for efficacy'
                                        #So it is easier to read in the seqeunce
                                        

                                        i = i.lower()
                                         

                                        #this makes a copy to work with to keep the original for whatever reason
                                        #
                                        # 
                                        dna2 = list(dna1)
                                        



                                        #replaces old codon overlapping PAM with new similar codon
                                        #
                                        #
                                        dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i


                                        #Now the neighboring codon will be altered since for efficacy 
                                        #and since its not a special codon
                                        #
                                        #
                                        Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                                        Core_Seq_2 = ''.join(Core_Seq_2)      
                                        

                                        print Core_Seq_2 + ' produces ' + hg19_bias_2[Core_Seq_2][0]
                                        
                                        g = Core_Seq_2
                                        
                                        
                                        codon_list_2 = []
                                        for h in hg19_bias_2:
                                            if hg19_bias_2[h][0]==hg19_bias_2[Core_Seq_2][0]:
                                                codon_list_2.append(h)


                                        for j in codon_list_2:
                                            if j == g:
                                                codon_list_2.remove(j) 

                                                k = hg19_bias_2[Core_Seq_2]
                                                del hg19_bias_2[Core_Seq_2]   


                                                freq_diff_2 = 0
                                                freq_diff_list_2 = [] 
                                            
                                                for l in codon_list_2:
                                                    

                                                    freq_diff_2 = f[1] - hg19_bias_2[l][1]
                                                    freq_diff_2 = math.fabs(freq_diff_2)
                                                    freq_diff_list_2.append(freq_diff_2)
                                                              
                                                    
                                                freq_diff_list_2.sort()
                                                z = freq_diff_list_2[0] 
                                                for l in codon_list_2:
                                                    if math.fabs(f[1] - hg19_bias_2[l][1]) == z: # and l != 'AGG':
                                       

                                                        if sys.argv[1].find('U') == -1:
                                                            l = l.replace('U', 'T')
                                                        elif sys.argv[1].find('T') == -1:
                                                            l = l

                                                        print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                                        
                                                        #So it is easier to read in the seqeunce
                                                        l = l.lower()
                                                         
                                                        
                                                        
                                                        # #replaces old codon overlapping PAM with new similar codon
                                                        dna2[a.end() - 9], dna2[a.end() - 8], dna2[a.end() - 7] = l

                                        #list becomes a string again for legibility 
                                        dna2 = ''.join(dna2)
                                        
                                        # print out the new DNA string
                                        if sys.argv[1].find('U') == -1:
                                            dna2 = dna2.replace('U', 'T')
                                            print dna2
                                            quit()
                                        else:
                                            print dna2   
                                            quit()

                            


                #Incase the codon is not special
                #
                #
                print t + ' produces ' + hg19_bias[PAM][0]
                    
             

                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[PAM][0]:
                        codon_list.append(i)

                 
    
                
                for b in codon_list:
                    if b == t:
                        codon_list.remove(b)
                                       
                        f = hg19_bias[PAM]
                        del hg19_bias[PAM]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x and i != 'AGG':
                               

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i

                                print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                
                                
                                i = i.lower()

                                dna2 = list(dna1)
                                
                                
                                dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = i
                                    
                                
                                dna2 = ''.join(dna2)
                                
                                
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                    print dna2
                                    quit()
                                else:
                                    print dna2   
                                    quit()
                   
            #This asks whether the end position is in a different frame based on the positions remainder
            #is it in a 231?
            #
            #
            elif (a.end()%3) == 1:
                        
                #everything after this point is the same as above for a different frame or checking the reverse complement of the target
                #
                #
                codon_PAM = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

                codon_PAM = ''.join(codon_PAM)
                        
                print'The PAM is in a 231 frame'
                
                if sys.argv[1].find('T') == -1:
                    t = codon_PAM.replace('T', 'U')
                elif sys.argv[1].find('U') == -1:
                    t = codon_PAM
                      
                print t + ' Is the overlapping codon'

                print t +' produces ' + hg19_bias[codon_PAM][0]
                

                if hg19_bias[codon_PAM][0] == 'W' or hg19_bias[codon_PAM][0] == 'M':
                    print 'Cannot alter Tryptophan or Methionine due to only one codon'
                    print 'Now attempting to alter a codon in the core seqeunce (4 bases upstream from PAM site)'
                    Core_Seq = dna1[a.end() - 4] + dna1[a.end() - 3] + dna1[a.end() - 2]
                                            
                    Core_Seq = ''.join(Core_Seq)
                    
                    print Core_Seq + ' produces ' + hg19_bias[Core_Seq][0]
                    
                    if hg19_bias[Core_Seq][0] == 'W' or hg19_bias[Core_Seq][0] == 'M':
                        print 'cannot swap tryptophan or methionine, now attempting to swap adjacent codon'

                        print 'Now swapping adjacent codon for efficacy'
                        #So it is easier to read in the seqeunce
                        #i = i.lower()
                         
                        dna2 = list(dna1)
                        
                        # #replaces old codon overlapping PAM with new similar codon
                        #dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = i

                        Core_Seq_2 = dna1[a.end() - 7] + dna1[a.end() - 6] + dna1[a.end() - 5]
                        Core_Seq_2 = ''.join(Core_Seq_2)      
                        

                        print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2][0]
                        
                        if hg19_bias[Core_Seq_2][0] == 'M' or hg19_bias[Core_Seq_2][0] == 'W':
                            print 'Cannot alter pam-site or both codons in coding sequence. You are out of luck :/'
                            quit()
                        
                        else:
                            

                            g = Core_Seq_2
                            #make a list of other codons that have the same amino acid
                            
                            codon_list_2 = []
                            for h in hg19_bias:
                                if hg19_bias[h][0]==hg19_bias[Core_Seq_2][0]:
                                    codon_list_2.append(h)


                            for j in codon_list_2:
                                if j == g:
                                    codon_list_2.remove(j) 

                                    k = hg19_bias[Core_Seq_2]
                                    del hg19_bias[Core_Seq_2]   


                                    freq_diff_2 = 0
                                    freq_diff_list_2 = [] 
                                
                                    for l in codon_list_2:
                                        

                                        freq_diff_2 = k[1] - hg19_bias[l][1]
                                        freq_diff_2 = math.fabs(freq_diff_2)
                                        freq_diff_list_2.append(freq_diff_2)
                                                  
                                        
                                    freq_diff_list_2.sort()
                                    z = freq_diff_list_2[0] 
                                    for l in codon_list_2:
                                        if math.fabs(k[1] - hg19_bias[l][1]) == z:
                           

                                            if sys.argv[1].find('U') == -1:
                                                l = l.replace('U', 'T')
                                            elif sys.argv[1].find('T') == -1:
                                                l = l

                                            print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                            
                                            #So it is easier to read in the seqeunce
                                            l = l.lower()
                                             
                                            
                                            
                                            # #replaces old codon overlapping PAM with new similar codon
                                            dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = l

                            #list becomes a string again for legibility 
                            dna2 = ''.join(dna2)
                            
                            # print out the new DNA string
                            if sys.argv[1].find('U') == -1:
                                dna2 = dna2.replace('U', 'T')
                                print dna2
                                quit()
                            else:
                                print dna2   
                                quit()

                    else:

                        t = Core_Seq
                        
                        
                        
                        codon_list = []
                        for i in hg19_bias:
                            if hg19_bias[i][0]==hg19_bias[Core_Seq][0]:
                                codon_list.append(i)
            
                        
                              
            
                        
                        for b in codon_list:
                            if b == t:
                                codon_list.remove(b)
                                               
                                f = hg19_bias[Core_Seq]
                                del hg19_bias[Core_Seq]
                                

                                freq_diff = 0
                                freq_diff_list = [] 
                            
                                for i in codon_list:
                                    

                                    freq_diff = f[1] - hg19_bias[i][1]
                                    freq_diff = math.fabs(freq_diff)
                                    freq_diff_list.append(freq_diff)
                                              
                                    
                                freq_diff_list.sort()
                                x = freq_diff_list[0] 
                                for i in codon_list:
                                    
                                    if math.fabs(f[1] - hg19_bias[i][1]) == x:
                                       

                                        if sys.argv[1].find('U') == -1:
                                            i = i.replace('U', 'T')
                                        elif sys.argv[1].find('T') == -1:
                                            i = i

                                        print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]


                                        i = i.lower()
                                         
                                        dna2 = list(dna1)
                                        

                                        dna2[a.end() - 4], dna2[a.end() - 3], dna2[a.end() - 2] = i
                                        dna2 = ''.join(dna2)
                                        
                                        Core_Seq_2 = dna1[a.end() - 7] + dna1[a.end() - 6] + dna1[a.end() - 5]
                                            
                                        if hg19_bias[Core_Seq_2][0] == 'W':    

                                            print 'cannot alter adjacent codon because it produces tryptophan'
                                            
                                            if sys.argv[1].find('U') == -1:
                                                dna2 = dna2.replace('U', 'T')
                                                print dna2
                                                quit()
                                            else:
                                                print dna2   
                                                quit()

                                           

                                        else:
                                            
                                            print 'Also swapping adjacent codon for efficacy'
                                            #So it is easier to read in the seqeunce
                                            i = i.lower()
                                             
                                            dna2 = list(dna2)
                                            
                                            # #replaces old codon overlapping PAM with new similar codon
                                            dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = i

                                            Core_Seq_2 = dna1[a.end() - 7] + dna1[a.end() - 6] + dna1[a.end() - 5]
                                            Core_Seq_2 = ''.join(Core_Seq_2)      
                                            

                                            print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2][0]
                                            
                                            g = Core_Seq_2
                                            #make a list of other codons that have the same amino acid
                                            
                                            codon_list_2 = []
                                            for h in hg19_bias:
                                                if hg19_bias[h][0]==hg19_bias[Core_Seq_2][0]:
                                                    codon_list_2.append(h)


                                            for j in codon_list_2:
                                                if j == g:
                                                    codon_list_2.remove(j) 

                                                    k = hg19_bias[Core_Seq_2]
                                                    del hg19_bias[Core_Seq_2]   


                                                    freq_diff_2 = 0
                                                    freq_diff_list_2 = [] 
                                                
                                                    for l in codon_list_2:
                                                        

                                                        freq_diff_2 = f[1] - hg19_bias[l][1]
                                                        freq_diff_2 = math.fabs(freq_diff_2)
                                                        freq_diff_list_2.append(freq_diff_2)
                                                                  
                                                        
                                                    freq_diff_list_2.sort()
                                                    z = freq_diff_list_2[0] 
                                                    for l in codon_list_2:
                                                        if math.fabs(f[1] - hg19_bias[l][1]) == z:
                                           

                                                            if sys.argv[1].find('U') == -1:
                                                                l = l.replace('U', 'T')
                                                            elif sys.argv[1].find('T') == -1:
                                                                l = l

                                                            print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                                            
                                                            #So it is easier to read in the seqeunce
                                                            l = l.lower()
                                                             
                                                            
                                                            
                                                            # #replaces old codon overlapping PAM with new similar codon
                                                            dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = l

                                            #list becomes a string again for legibility 
                                            dna2 = ''.join(dna2)
                                            
                                            # print out the new DNA string
                                            if sys.argv[1].find('U') == -1:
                                                dna2 = dna2.replace('U', 'T')
                                                print dna2
                                                quit()
                                            else:
                                                print dna2   
                                                quit()





                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[codon_PAM][0]:
                        codon_list.append(i)

                
                for b in codon_list:
                    if b == t:
                        codon_list.remove(b)
              
                        f = hg19_bias[codon_PAM]
                        del hg19_bias[codon_PAM]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x:

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i
                            
                                print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                
                                i = i.lower()

                                dna2 = list(dna1)

                                dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = i
                                
                                dna2 = ''.join(dna2)
                            
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                print dna2
                                quit()



                
            else:
                print 'PAM is in 312 frame. You cannot alter the PAM without altering the future peptide'
                
                   
                print 'Now attempting to alter a codon in the core seqeunce (4 bases upstream from PAM site)'
                Core_Seq = dna1[a.end() - 5] + dna1[a.end() - 4] + dna1[a.end() - 3]
                 
                

                Core_Seq = ''.join(Core_Seq)
                
                CS = Core_Seq

                if sys.argv[1].find('U') == -1:
                    CS = CS.replace('U', 'T')
                elif sys.argv[1].find('T') == -1:
                    CS = CS


                print CS + ' produces ' + hg19_bias[Core_Seq][0]
                
                if hg19_bias[Core_Seq][0] == 'W':
                    print 'Only one codon produces tryptophan and the codon cannot be altered'
                   
                    
                    print 'Now swapping adjacent codon for efficacy'
                                
                    

                                              

                    Core_Seq_2 = dna1[a.end() - 8] + dna1[a.end() - 7] + dna1[a.end() - 6]
                    Core_Seq_2 = ''.join(Core_Seq_2)      
                    
                    Core_Seq_2_2 = Core_Seq_2

                    if sys.argv[1].find('U') == -1:
                        Core_Seq_2 = Core_Seq_2.replace('U', 'T')
                    elif sys.argv[1].find('T') == -1:
                        Core_Seq_2 = Core_Seq_2

                    print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2_2][0]
                    
                    if hg19_bias[Core_Seq_2_2][0] == 'W':
                        print 'Both codons in core sequence are Tryptophan and the Target sequence cannot be altered at this time'
                        
                        quit()


                    g = Core_Seq_2_2
                        #make a list of other codons that have the same amino acid
                        
                    codon_list_2 = []
                    for h in hg19_bias:
                        if hg19_bias[h][0]==hg19_bias[Core_Seq_2_2][0]:
                            codon_list_2.append(h)


                    for j in codon_list_2:
                        if j == g:
                            codon_list_2.remove(j) 


                            k = hg19_bias[Core_Seq_2_2]
                            del hg19_bias[Core_Seq_2_2]   


                            freq_diff_2 = 0
                            freq_diff_list_2 = [] 
                        
                            for l in codon_list_2:
                                

                                freq_diff_2 = k[1] - hg19_bias[l][1]
                                freq_diff_2 = math.fabs(freq_diff_2)
                                freq_diff_list_2.append(freq_diff_2)
                                          
                                
                            freq_diff_list_2.sort()
                            z = freq_diff_list_2[0] 
                            for l in codon_list_2:
                                if math.fabs(k[1] - hg19_bias[l][1]) == z:
                   

                                    if sys.argv[1].find('U') == -1:
                                        l = l.replace('U', 'T')
                                    elif sys.argv[1].find('T') == -1:
                                        l = l

                                    print 'Now swapping ' + Core_Seq_2 + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                    
                                    #So it is easier to read in the seqeunce
                                    l = l.lower()
                                    
                                    dna2 = dna1

                                    if dna2.find('U') == -1:
                                        dna2 = dna2.replace('U', 'T')
                                    elif dna2.find('T') == -1:
                                        dna2 = dna2

                                    dna2 = list(dna1)

                                    # #replaces old codon overlapping PAM with new similar codon
                                    dna2[a.end() - 8], dna2[a.end() - 7], dna2[a.end() - 6] = l

                    #list becomes a string again for legibility 
                    dna2 = ''.join(dna2)
                    
                    # print out the new DNA string
                    if sys.argv[1].find('U') == -1:
                        dna1 = dna2.replace('U', 'T')
                        print dna1
                        quit()
                    else:
                        print dna2  
                        quit()

                t = Core_Seq
                # # g = Core_Seq_2
                # #make a list of other codons that have the same amino acid
                hg19_bias_2 = hg19_bias
                

                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[Core_Seq][0]:
                        codon_list.append(i)
    
                
                # #finds a different codon with the highest codon bias       
                codon_list_temp = codon_list
                
                for b in codon_list_temp:
                    if b == t:
                        codon_list_temp.remove(b)
                                       
                        f = hg19_bias[Core_Seq]
                        del hg19_bias[Core_Seq]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list_temp:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list_temp:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x:
                               

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i

                                print 'Now swapping ' + CS + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                

                                print 'Swapping adjacent codon for efficacy'
                                
                                i = i.lower()
                                     
                                dna2 = dna1             
                                dna2 = list(dna2)


                                dna2[a.end() - 5], dna2[a.end() - 4], dna2[a.end() - 3] = i

                                                             

                                Core_Seq_2 = dna1[a.end() - 8] + dna1[a.end() - 7] + dna1[a.end() - 6]
                                #Core_Seq_2 = ''.join(Core_Seq_2)      
                                                                
                                
                                if Core_Seq == Core_Seq_2:
                                    print 'Adjecent codon is the same and will be swapped similarly'

                                    dna1 = list(dna1)
                                    dna1[a.end() - 5], dna1[a.end() - 4], dna1[a.end() - 3] = i
                                    dna1[a.end() - 8], dna1[a.end() - 7], dna1[a.end() - 6] = i
                                    dna1 = ''.join(dna1)
                                    print dna1
                                    quit()

                                if hg19_bias_2[Core_Seq_2][0] == 'W':
                                    print 'Cannot alter adject codon. Target DNA is complete'
                                    dna2 = ''.join(dna2)
                                

                                    if sys.argv[1].find('U') == -1:
                                        #dna2 = ''.join(dna2)
                                        dna2 = dna2.replace('U', 'T')
                                        print dna2
                                        quit()
                                    else:
                                        print dna2  
                                        quit()



                                    
                                Core_Seq_2_2 = Core_Seq_2

                                if sys.argv[1].find('U') == -1:
                                    Core_Seq_2 = Core_Seq_2.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    Core_Seq_2 = Core_Seq_2


                                print Core_Seq_2 + ' produces ' + hg19_bias_2[Core_Seq_2_2][0]
                                

                                g = Core_Seq_2_2
                                #make a list of other codons that have the same amino acid
                                    
                                codon_list_2 = []
                                for h in hg19_bias_2:
                                    if hg19_bias_2[h][0]==hg19_bias_2[Core_Seq_2_2][0]:
                                        codon_list_2.append(h)


                                for j in codon_list_2:
                                    if j == g:
                                        codon_list_2.remove(j) 


                                        k = hg19_bias_2[Core_Seq_2_2]
                                        del hg19_bias_2[Core_Seq_2_2]   


                                        freq_diff_2 = 0
                                        freq_diff_list_2 = [] 
                                    
                                        for l in codon_list_2:
                                            

                                            freq_diff_2 = k[1] - hg19_bias_2[l][1]
                                            freq_diff_2 = math.fabs(freq_diff_2)
                                            freq_diff_list_2.append(freq_diff_2)
                                                      
                                            
                                        freq_diff_list_2.sort()
                                        z = freq_diff_list_2[0] 
                                        for l in codon_list_2:
                                            if math.fabs(k[1] - hg19_bias_2[l][1]) == z:
                               

                                                if sys.argv[1].find('U') == -1:
                                                    l = l.replace('U', 'T')
                                                elif sys.argv[1].find('T') == -1:
                                                    l = l

                                                print 'Now swapping ' + Core_Seq_2 + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                                
                                                #So it is easier to read in the seqeunce
                                                l = l.lower()
                                                 
                                                dna2 = list(dna2)
                                                
                                                # #replaces old codon overlapping PAM with new similar codon
                                                dna2[a.end() - 8], dna2[a.end() - 7], dna2[a.end() - 6] = l
                                                dna2[a.end() - 5], dna2[a.end() - 4], dna2[a.end() - 3] = i

                                #list becomes a string again for legibility 
                                dna2 = ''.join(dna2)
                                
                                # print out the new DNA string
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                    print dna2
                                    quit()
                                else:
                                    print dna2   
                                    quit()
                                
                                
###########################################################################
###########################################################################
###########################################################################
###########################################################################################################

        else:

            print 'Target not adjacent to a PAM site'

elif dna1.count(gRNA) == 0 and revcom_dna1.count(gRNA) == 1:
    
    
    print 'Reverse-complementing the target :'
    
    if sys.argv[1].find('U') == -1:
        rD = revcom_dna1.replace('U', 'T')
        rG = gRNA.replace('U', 'T')
    elif sys.argv[1].find('T') == -1:
        rD = revcom_dna1
        rG = gRNA
        
    #print 'This is the target:\n'+ rD
    
    print 'This is the guide:\n' + rG    

    dna1 = revcom_dna1
    gRNA = gRNA
    for a in re.finditer(gRNA, dna1):    

    #This makes sure the adjacent 3 bases make a PAM site
    #
    #
        if dna1[a.end() + 1] == 'G' and dna1[a.end() + 2] == 'G':
            
            #This makes the PAM site a new PAM variable
            #
            #
            PAM = dna1[a.end()], dna1[a.end() + 1], dna1[a.end() + 2]

            PAM = ''.join(PAM)

            print 'The guide matches the sequence once at: '
            print (a.start(), a.end())

            #This print's out the same codon/PAM in DNA if DNA was given or leaves it as RNA if not.
            #
            #
            if sys.argv[1].find('U') == -1:
                t = PAM.replace('U', 'T')
            elif sys.argv[1].find('T') == -1:
                t = PAM    
            print t + ' is the PAM sequence'
                
            #remainder indirectly gives postion for base 1 in PAM telling us the frame its in
            #
            #
            if (a.end()%3) == 0:
                print'The PAM is in a 123 frame'
                
                #if 123 frame, PAM is the same as the codon that needs to be changed
                #
                #
                print t + ' Is the overlapping codon'
                

                #Incase the overlapping codon is special and does not have replacements
                #The code will alter what's known as the core sequence
                #Changing the core should have a similar effect as altering the PAM
                #
                #
                if hg19_bias[PAM][0] == 'W' or hg19_bias[PAM][0] == 'M':
                    print 'Cannot alter Tryptophan or Methionine due to only one codon'
                    print 'Now attempting to alter a codon in the core seqeunce (4 bases upstream from PAM site)'
                    Core_Seq = dna1[a.end() - 6] + dna1[a.end() - 5] + dna1[a.end() - 4]
                                            
                    Core_Seq = ''.join(Core_Seq)
                    
                    print Core_Seq + ' produces ' + hg19_bias[Core_Seq][0]
                    

                    if hg19_bias[Core_Seq][0] == 'W':

                        print 'Cannot alter a portion of the core sequence codon due to Tryptophan'

                        Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                        
                        i = Core_Seq_2

                        if hg19_bias[Core_Seq_2][0] == 'W':



                            #for finding the swap easily in the sequence
                            #
                            #
                            i = i.lower()
                         
                            dna2 = list(dna1)
                        
                            #The first codon near the core sequence is changed but not the second
                            #
                            #
                            #dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i
                            dna2 = ''.join(dna2)
                            print 'cannot alter adjacent codon because it also produces tryptophan'
                            print 'Sorry, but the target sequence cannot be adjusted to block cas9 at this time.'
                            
                            
                            quit()

                        else:

                            print 'Now swapping adjacent codon for efficacy'
                            #So it is easier to read in the seqeunce
                            

                            i = i.lower()
                             

                            #this makes a copy to work with to keep the original for whatever reason
                            #
                            # 
                            dna2 = list(dna1)
                            



                            #replaces old codon overlapping PAM with new similar codon
                            #
                            #
                            #dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i


                            #Now the neighboring codon will be altered since for efficacy 
                            #and since its not a special codon
                            #
                            #
                            Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                            Core_Seq_2 = ''.join(Core_Seq_2)      
                            

                            print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2][0]
                            
                            g = Core_Seq_2
                            
                            
                            codon_list_2 = []
                            for h in hg19_bias:
                                if hg19_bias[h][0]==hg19_bias[Core_Seq_2][0]:
                                    codon_list_2.append(h)


                            for j in codon_list_2:
                                if j == g:
                                    codon_list_2.remove(j) 

                                    k = hg19_bias[Core_Seq_2]
                                    del hg19_bias[Core_Seq_2]   


                                    freq_diff_2 = 0
                                    freq_diff_list_2 = [] 
                                
                                    for l in codon_list_2:
                                        

                                        freq_diff_2 = k[1] - hg19_bias[l][1]
                                        freq_diff_2 = math.fabs(freq_diff_2)
                                        freq_diff_list_2.append(freq_diff_2)
                                                  
                                        
                                    freq_diff_list_2.sort()
                                    z = freq_diff_list_2[0] 
                                    for l in codon_list_2:
                                        if math.fabs(k[1] - hg19_bias[l][1]) == z: # and l != 'AGG':
                           

                                            if sys.argv[1].find('U') == -1:
                                                l = l.replace('U', 'T')
                                            elif sys.argv[1].find('T') == -1:
                                                l = l

                                            print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                            
                                            #So it is easier to read in the seqeunce
                                            l = l.lower()
                                             
                                            
                                            
                                            # #replaces old codon overlapping PAM with new similar codon
                                            dna2[a.end() - 9], dna2[a.end() - 8], dna2[a.end() - 7] = l

                                            #list becomes a string again for legibility 
                                            dna2 = ''.join(dna2)
                            
                                            # print out the new DNA string
                                            if sys.argv[1].find('U') == -1:
                                                dna2 = dna2.replace('U', 'T')
                                                print dna2
                                                quit()
                                            else:
                                                print dna2   
                                                quit()

                        
                    
                    t = Core_Seq
                    
                    #make a list of other codons that have the same amino acid
                    #
                    #
                    codon_list = []
                    for i in hg19_bias:
                        if hg19_bias[i][0]==hg19_bias[Core_Seq][0]:
                            codon_list.append(i)
        
                    
                    hg19_bias_2 = hg19_bias

                    #finds a different codon with the highest codon bias       
                    #
                    #
                    for b in codon_list:
                        if b == t:
                            codon_list.remove(b)
                                           
                            f = hg19_bias[Core_Seq]
                            del hg19_bias[Core_Seq]
                            

                            freq_diff = 0
                            freq_diff_list = [] 
                        
                            #Here we created a list that holds frequencie differences between
                            #a new codon and the original
                            #the list is sorted for each codon and the codon with the smallest difference
                            #is selected to replace the original
                            #
                            #
                            for i in codon_list:
                                

                                freq_diff = f[1] - hg19_bias[i][1]
                                freq_diff = math.fabs(freq_diff)
                                freq_diff_list.append(freq_diff)
                                          
                            #sort the list from smallest to greatest
                            #the first on the list has a frequency closest to that of original
                            #
                            #    
                            freq_diff_list.sort()
                            x = freq_diff_list[0] 
                            for i in codon_list:
                                

                                #Select the first codon in the list
                                #
                                #
                                if math.fabs(f[1] - hg19_bias[i][1]) == x:
                                   

                                    if sys.argv[1].find('U') == -1:
                                        i = i.replace('U', 'T')
                                    elif sys.argv[1].find('T') == -1:
                                        i = i

                                    print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                    
                                    i = i.lower()

                                    dna2 = list(dna1)
                                    
                                    dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i

                                        

                                    #So now 1 codon near the core sequence has been changed
                                    #This is good but changes in all the codons overlapping the core is more beneficial
                                    #Now we will change the neighboring codon unless its a special codon that cannot be
                                    #
                                    #
                                    Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                                    if Core_Seq_2 == Core_Seq:

                                        print 'The adjecent codon is the same and will also be changed for efficacy'
                                        dna2[a.end() - 9], dna2[a.end() - 8], dna2[a.end() - 7] = i
                                        
                                        dna2 = ''.join(dna2)
                                        if sys.argv[1].find('U') == -1:
                                            dna2 = dna2.replace('U', 'T')

                                            print dna2
                                            quit()

                                    if hg19_bias_2[Core_Seq_2][0] == 'W':



                                        print 'cannot alter adjacent codon because it produces tryptophan'
                                        
                                        dna2 = ''.join(dna2)
                                        if sys.argv[1].find('U') == -1:
                                            dna2 = dna2.replace('U', 'T')

                                        print dna2
                                        quit()

                                    else:

                                        print 'Also swapping adjacent codon for efficacy'
                                        #So it is easier to read in the seqeunce
                                        

                                        i = i.lower()
                                         

                                        #this makes a copy to work with to keep the original for whatever reason
                                        #
                                        # 
                                        dna2 = list(dna1)
                                        



                                        #replaces old codon overlapping PAM with new similar codon
                                        #
                                        #
                                        dna2[a.end() - 6], dna2[a.end() - 5], dna2[a.end() - 4] = i


                                        #Now the neighboring codon will be altered since for efficacy 
                                        #and since its not a special codon
                                        #
                                        #
                                        Core_Seq_2 = dna1[a.end() - 9] + dna1[a.end() - 8] + dna1[a.end() - 7]
                                        Core_Seq_2 = ''.join(Core_Seq_2)      
                                        

                                        print Core_Seq_2 + ' produces ' + hg19_bias_2[Core_Seq_2][0]
                                        
                                        g = Core_Seq_2
                                        
                                        
                                        codon_list_2 = []
                                        for h in hg19_bias_2:
                                            if hg19_bias_2[h][0]==hg19_bias_2[Core_Seq_2][0]:
                                                codon_list_2.append(h)


                                        for j in codon_list_2:
                                            if j == g:
                                                codon_list_2.remove(j) 

                                                k = hg19_bias_2[Core_Seq_2]
                                                del hg19_bias_2[Core_Seq_2]   


                                                freq_diff_2 = 0
                                                freq_diff_list_2 = [] 
                                            
                                                for l in codon_list_2:
                                                    

                                                    freq_diff_2 = f[1] - hg19_bias_2[l][1]
                                                    freq_diff_2 = math.fabs(freq_diff_2)
                                                    freq_diff_list_2.append(freq_diff_2)
                                                              
                                                    
                                                freq_diff_list_2.sort()
                                                z = freq_diff_list_2[0] 
                                                for l in codon_list_2:
                                                    if math.fabs(f[1] - hg19_bias_2[l][1]) == z: # and l != 'AGG':
                                       

                                                        if sys.argv[1].find('U') == -1:
                                                            l = l.replace('U', 'T')
                                                        elif sys.argv[1].find('T') == -1:
                                                            l = l

                                                        print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                                        
                                                        #So it is easier to read in the seqeunce
                                                        l = l.lower()
                                                         
                                                        
                                                        
                                                        # #replaces old codon overlapping PAM with new similar codon
                                                        dna2[a.end() - 9], dna2[a.end() - 8], dna2[a.end() - 7] = l

                                        #list becomes a string again for legibility 
                                        dna2 = ''.join(dna2)
                                        
                                        # print out the new DNA string
                                        if sys.argv[1].find('U') == -1:
                                            dna2 = dna2.replace('U', 'T')
                                            print dna2
                                            quit()
                                        else:
                                            print dna2   
                                            quit()

                            


                #Incase the codon is not special
                #
                #
                print t + ' produces ' + hg19_bias[PAM][0]
                    
             

                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[PAM][0]:
                        codon_list.append(i)

                 

                
                for b in codon_list:
                    if b == t:
                        codon_list.remove(b)
                                       
                        f = hg19_bias[PAM]
                        del hg19_bias[PAM]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x and i != 'AGG':
                               

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i

                                print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                
                                
                                i = i.lower()

                                dna2 = list(dna1)
                                
                                
                                dna2[a.end()], dna2[a.end() + 1], dna2[a.end() + 2] = i
                                    
                                
                                dna2 = ''.join(dna2)
                                
                                
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                    print dna2
                                    quit()
                                else:
                                    print dna2   
                                    quit()
                   
            #This asks whether the end position is in a different frame based on the positions remainder
            #is it in a 231?
            #
            #
            elif (a.end()%3) == 1:
                        
                #everything after this point is the same as above for a different frame or checking the reverse complement of the target
                #
                #
                codon_PAM = dna1[a.end() - 1], dna1[a.end()], dna1[a.end() + 1]

                codon_PAM = ''.join(codon_PAM)
                        
                print'The PAM is in a 231 frame'
                
                if sys.argv[1].find('T') == -1:
                    t = codon_PAM.replace('T', 'U')
                elif sys.argv[1].find('U') == -1:
                    t = codon_PAM
                      
                print t + ' Is the overlapping codon'

                print t +' produces ' + hg19_bias[codon_PAM][0]
                

                if hg19_bias[codon_PAM][0] == 'W' or hg19_bias[codon_PAM][0] == 'M':
                    print 'Cannot alter Tryptophan or Methionine due to only one codon'
                    print 'Now attempting to alter a codon in the core seqeunce (4 bases upstream from PAM site)'
                    Core_Seq = dna1[a.end() - 4] + dna1[a.end() - 3] + dna1[a.end() - 2]
                                            
                    Core_Seq = ''.join(Core_Seq)
                    
                    print Core_Seq + ' produces ' + hg19_bias[Core_Seq][0]
                    
                    if hg19_bias[Core_Seq][0] == 'W' or hg19_bias[Core_Seq][0] == 'M':
                        print 'cannot swap tryptophan or methionine, now attempting to swap adjacent codon'

                        print 'Now swapping adjacent codon for efficacy'
                        #So it is easier to read in the seqeunce
                        #i = i.lower()
                         
                        dna2 = list(dna1)
                        
                        # #replaces old codon overlapping PAM with new similar codon
                        #dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = i

                        Core_Seq_2 = dna1[a.end() - 7] + dna1[a.end() - 6] + dna1[a.end() - 5]
                        Core_Seq_2 = ''.join(Core_Seq_2)      
                        

                        print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2][0]
                        
                        if hg19_bias[Core_Seq_2][0] == 'M' or hg19_bias[Core_Seq_2][0] == 'W':
                            print 'Cannot alter pam-site or both codons in coding sequence. You are out of luck :/'
                            quit()
                        
                        else:
                            

                            g = Core_Seq_2
                            #make a list of other codons that have the same amino acid
                            
                            codon_list_2 = []
                            for h in hg19_bias:
                                if hg19_bias[h][0]==hg19_bias[Core_Seq_2][0]:
                                    codon_list_2.append(h)


                            for j in codon_list_2:
                                if j == g:
                                    codon_list_2.remove(j) 

                                    k = hg19_bias[Core_Seq_2]
                                    del hg19_bias[Core_Seq_2]   


                                    freq_diff_2 = 0
                                    freq_diff_list_2 = [] 
                                
                                    for l in codon_list_2:
                                        

                                        freq_diff_2 = k[1] - hg19_bias[l][1]
                                        freq_diff_2 = math.fabs(freq_diff_2)
                                        freq_diff_list_2.append(freq_diff_2)
                                                  
                                        
                                    freq_diff_list_2.sort()
                                    z = freq_diff_list_2[0] 
                                    for l in codon_list_2:
                                        if math.fabs(k[1] - hg19_bias[l][1]) == z:
                           

                                            if sys.argv[1].find('U') == -1:
                                                l = l.replace('U', 'T')
                                            elif sys.argv[1].find('T') == -1:
                                                l = l

                                            print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                            
                                            #So it is easier to read in the seqeunce
                                            l = l.lower()
                                             
                                            
                                            
                                            # #replaces old codon overlapping PAM with new similar codon
                                            dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = l

                            #list becomes a string again for legibility 
                            dna2 = ''.join(dna2)
                            
                            # print out the new DNA string
                            if sys.argv[1].find('U') == -1:
                                dna2 = dna2.replace('U', 'T')
                                print dna2
                                quit()
                            else:
                                print dna2   
                                quit()

                    else:

                        t = Core_Seq
                        
                        
                        
                        codon_list = []
                        for i in hg19_bias:
                            if hg19_bias[i][0]==hg19_bias[Core_Seq][0]:
                                codon_list.append(i)
            
                        
                              
            
                        
                        for b in codon_list:
                            if b == t:
                                codon_list.remove(b)
                                               
                                f = hg19_bias[Core_Seq]
                                del hg19_bias[Core_Seq]
                                

                                freq_diff = 0
                                freq_diff_list = [] 
                            
                                for i in codon_list:
                                    

                                    freq_diff = f[1] - hg19_bias[i][1]
                                    freq_diff = math.fabs(freq_diff)
                                    freq_diff_list.append(freq_diff)
                                              
                                    
                                freq_diff_list.sort()
                                x = freq_diff_list[0] 
                                for i in codon_list:
                                    
                                    if math.fabs(f[1] - hg19_bias[i][1]) == x:
                                       

                                        if sys.argv[1].find('U') == -1:
                                            i = i.replace('U', 'T')
                                        elif sys.argv[1].find('T') == -1:
                                            i = i

                                        print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]


                                        i = i.lower()
                                         
                                        dna2 = list(dna1)
                                        

                                        dna2[a.end() - 4], dna2[a.end() - 3], dna2[a.end() - 2] = i
                                        dna2 = ''.join(dna2)
                                        
                                        Core_Seq_2 = dna1[a.end() - 7] + dna1[a.end() - 6] + dna1[a.end() - 5]
                                            
                                        if hg19_bias[Core_Seq_2][0] == 'W':    

                                            print 'cannot alter adjacent codon because it produces tryptophan'
                                            
                                            if sys.argv[1].find('U') == -1:
                                                dna2 = dna2.replace('U', 'T')
                                                print dna2
                                                quit()
                                            else:
                                                print dna2   
                                                quit()

                                           

                                        else:
                                            
                                            print 'Also swapping adjacent codon for efficacy'
                                            #So it is easier to read in the seqeunce
                                            i = i.lower()
                                             
                                            dna2 = list(dna2)
                                            
                                            # #replaces old codon overlapping PAM with new similar codon
                                            dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = i

                                            Core_Seq_2 = dna1[a.end() - 7] + dna1[a.end() - 6] + dna1[a.end() - 5]
                                            Core_Seq_2 = ''.join(Core_Seq_2)      
                                            

                                            print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2][0]
                                            
                                            g = Core_Seq_2
                                            #make a list of other codons that have the same amino acid
                                            
                                            codon_list_2 = []
                                            for h in hg19_bias:
                                                if hg19_bias[h][0]==hg19_bias[Core_Seq_2][0]:
                                                    codon_list_2.append(h)


                                            for j in codon_list_2:
                                                if j == g:
                                                    codon_list_2.remove(j) 

                                                    k = hg19_bias[Core_Seq_2]
                                                    del hg19_bias[Core_Seq_2]   


                                                    freq_diff_2 = 0
                                                    freq_diff_list_2 = [] 
                                                
                                                    for l in codon_list_2:
                                                        

                                                        freq_diff_2 = f[1] - hg19_bias[l][1]
                                                        freq_diff_2 = math.fabs(freq_diff_2)
                                                        freq_diff_list_2.append(freq_diff_2)
                                                                  
                                                        
                                                    freq_diff_list_2.sort()
                                                    z = freq_diff_list_2[0] 
                                                    for l in codon_list_2:
                                                        if math.fabs(f[1] - hg19_bias[l][1]) == z:
                                           

                                                            if sys.argv[1].find('U') == -1:
                                                                l = l.replace('U', 'T')
                                                            elif sys.argv[1].find('T') == -1:
                                                                l = l

                                                            print 'Now swapping ' + g + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                                            
                                                            #So it is easier to read in the seqeunce
                                                            l = l.lower()
                                                             
                                                            
                                                            
                                                            # #replaces old codon overlapping PAM with new similar codon
                                                            dna2[a.end() - 7], dna2[a.end() - 6], dna2[a.end() - 5] = l

                                            #list becomes a string again for legibility 
                                            dna2 = ''.join(dna2)
                                            
                                            # print out the new DNA string
                                            if sys.argv[1].find('U') == -1:
                                                dna2 = dna2.replace('U', 'T')
                                                print dna2
                                                quit()
                                            else:
                                                print dna2   
                                                quit()





                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[codon_PAM][0]:
                        codon_list.append(i)

                
                for b in codon_list:
                    if b == t:
                        codon_list.remove(b)
              
                        f = hg19_bias[codon_PAM]
                        del hg19_bias[codon_PAM]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x:

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i
                            
                                print 'Now swapping ' + t + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                
                                i = i.lower()

                                dna2 = list(dna1)

                                dna2[a.end() - 1], dna2[a.end()], dna2[a.end() + 1] = i
                                
                                dna2 = ''.join(dna2)
                            
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                print dna2
                                quit()



                
            else:
                print 'PAM is in 312 frame. You cannot alter the PAM without altering the future peptide'
                
                   
                print 'Now attempting to alter a codon in the core seqeunce (4 bases upstream from PAM site)'
                Core_Seq = dna1[a.end() - 5] + dna1[a.end() - 4] + dna1[a.end() - 3]
                 
                

                Core_Seq = ''.join(Core_Seq)
                
                CS = Core_Seq

                if sys.argv[1].find('U') == -1:
                    CS = CS.replace('U', 'T')
                elif sys.argv[1].find('T') == -1:
                    CS = CS


                print CS + ' produces ' + hg19_bias[Core_Seq][0]
                
                if hg19_bias[Core_Seq][0] == 'W':
                    print 'Only one codon produces tryptophan and the codon cannot be altered'
                   
                    
                    print 'Now swapping adjacent codon for efficacy'
                                
                    

                                              

                    Core_Seq_2 = dna1[a.end() - 8] + dna1[a.end() - 7] + dna1[a.end() - 6]
                    Core_Seq_2 = ''.join(Core_Seq_2)      
                    
                    Core_Seq_2_2 = Core_Seq_2

                    if sys.argv[1].find('U') == -1:
                        Core_Seq_2 = Core_Seq_2.replace('U', 'T')
                    elif sys.argv[1].find('T') == -1:
                        Core_Seq_2 = Core_Seq_2

                    print Core_Seq_2 + ' produces ' + hg19_bias[Core_Seq_2_2][0]
                    
                    if hg19_bias[Core_Seq_2_2][0] == 'W':
                        print 'Both codons in core sequence are Tryptophan and the Target sequence cannot be altered at this time'
                        
                        quit()


                    g = Core_Seq_2_2
                        #make a list of other codons that have the same amino acid
                        
                    codon_list_2 = []
                    for h in hg19_bias:
                        if hg19_bias[h][0]==hg19_bias[Core_Seq_2_2][0]:
                            codon_list_2.append(h)


                    for j in codon_list_2:
                        if j == g:
                            codon_list_2.remove(j) 


                            k = hg19_bias[Core_Seq_2_2]
                            del hg19_bias[Core_Seq_2_2]   


                            freq_diff_2 = 0
                            freq_diff_list_2 = [] 
                        
                            for l in codon_list_2:
                                

                                freq_diff_2 = k[1] - hg19_bias[l][1]
                                freq_diff_2 = math.fabs(freq_diff_2)
                                freq_diff_list_2.append(freq_diff_2)
                                          
                                
                            freq_diff_list_2.sort()
                            z = freq_diff_list_2[0] 
                            for l in codon_list_2:
                                if math.fabs(k[1] - hg19_bias[l][1]) == z:
                   

                                    if sys.argv[1].find('U') == -1:
                                        l = l.replace('U', 'T')
                                    elif sys.argv[1].find('T') == -1:
                                        l = l

                                    print 'Now swapping ' + Core_Seq_2 + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                    
                                    #So it is easier to read in the seqeunce
                                    l = l.lower()
                                    
                                    dna2 = dna1

                                    if dna2.find('U') == -1:
                                        dna2 = dna2.replace('U', 'T')
                                    elif dna2.find('T') == -1:
                                        dna2 = dna2

                                    dna2 = list(dna1)

                                    # #replaces old codon overlapping PAM with new similar codon
                                    dna2[a.end() - 8], dna2[a.end() - 7], dna2[a.end() - 6] = l

                    #list becomes a string again for legibility 
                    dna2 = ''.join(dna2)
                    
                    # print out the new DNA string
                    if sys.argv[1].find('U') == -1:
                        dna1 = dna2.replace('U', 'T')
                        print dna1
                        quit()
                    else:
                        print dna2  
                        quit()

                t = Core_Seq
                # # g = Core_Seq_2
                # #make a list of other codons that have the same amino acid
                hg19_bias_2 = hg19_bias
                

                codon_list = []
                for i in hg19_bias:
                    if hg19_bias[i][0]==hg19_bias[Core_Seq][0]:
                        codon_list.append(i)

                
                # #finds a different codon with the highest codon bias       
                codon_list_temp = codon_list
                
                for b in codon_list_temp:
                    if b == t:
                        codon_list_temp.remove(b)
                                       
                        f = hg19_bias[Core_Seq]
                        del hg19_bias[Core_Seq]
                        

                        freq_diff = 0
                        freq_diff_list = [] 
                    
                        for i in codon_list_temp:
                            

                            freq_diff = f[1] - hg19_bias[i][1]
                            freq_diff = math.fabs(freq_diff)
                            freq_diff_list.append(freq_diff)
                                      
                            
                        freq_diff_list.sort()
                        x = freq_diff_list[0] 
                        for i in codon_list_temp:
                            
                            if math.fabs(f[1] - hg19_bias[i][1]) == x:
                               

                                if sys.argv[1].find('U') == -1:
                                    i = i.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    i = i

                                print 'Now swapping ' + CS + ' for ' + i + ' because its frequencies are close and both produce ' + f[0]
                                

                                print 'Swapping adjacent codon for efficacy'
                                
                                i = i.lower()
                                     
                                dna2 = dna1             
                                dna2 = list(dna2)


                                dna2[a.end() - 5], dna2[a.end() - 4], dna2[a.end() - 3] = i

                                                             

                                Core_Seq_2 = dna1[a.end() - 8] + dna1[a.end() - 7] + dna1[a.end() - 6]
                                #Core_Seq_2 = ''.join(Core_Seq_2)      
                                                                
                                
                                if Core_Seq == Core_Seq_2:
                                    print 'Adjecent codon is the same and will be swapped similarly'

                                    dna1 = list(dna1)
                                    dna1[a.end() - 5], dna1[a.end() - 4], dna1[a.end() - 3] = i
                                    dna1[a.end() - 8], dna1[a.end() - 7], dna1[a.end() - 6] = i
                                    dna1 = ''.join(dna1)
                                    print dna1
                                    quit()

                                if hg19_bias_2[Core_Seq_2][0] == 'W':
                                    print 'Cannot alter adject codon. Target DNA is complete'
                                    dna2 = ''.join(dna2)
                                

                                    if sys.argv[1].find('U') == -1:
                                        #dna2 = ''.join(dna2)
                                        dna2 = dna2.replace('U', 'T')
                                        print dna2
                                        quit()
                                    else:
                                        print dna2  
                                        quit()



                                    
                                Core_Seq_2_2 = Core_Seq_2

                                if sys.argv[1].find('U') == -1:
                                    Core_Seq_2 = Core_Seq_2.replace('U', 'T')
                                elif sys.argv[1].find('T') == -1:
                                    Core_Seq_2 = Core_Seq_2


                                print Core_Seq_2 + ' produces ' + hg19_bias_2[Core_Seq_2_2][0]
                                

                                g = Core_Seq_2_2
                                #make a list of other codons that have the same amino acid
                                    
                                codon_list_2 = []
                                for h in hg19_bias_2:
                                    if hg19_bias_2[h][0]==hg19_bias_2[Core_Seq_2_2][0]:
                                        codon_list_2.append(h)


                                for j in codon_list_2:
                                    if j == g:
                                        codon_list_2.remove(j) 


                                        k = hg19_bias_2[Core_Seq_2_2]
                                        del hg19_bias_2[Core_Seq_2_2]   


                                        freq_diff_2 = 0
                                        freq_diff_list_2 = [] 
                                    
                                        for l in codon_list_2:
                                            

                                            freq_diff_2 = k[1] - hg19_bias_2[l][1]
                                            freq_diff_2 = math.fabs(freq_diff_2)
                                            freq_diff_list_2.append(freq_diff_2)
                                                      
                                            
                                        freq_diff_list_2.sort()
                                        z = freq_diff_list_2[0] 
                                        for l in codon_list_2:
                                            if math.fabs(k[1] - hg19_bias_2[l][1]) == z:
                               

                                                if sys.argv[1].find('U') == -1:
                                                    l = l.replace('U', 'T')
                                                elif sys.argv[1].find('T') == -1:
                                                    l = l

                                                print 'Now swapping ' + Core_Seq_2 + ' for ' + l + ' because its frequencies are close and both produce ' + k[0]
                                                
                                                #So it is easier to read in the seqeunce
                                                l = l.lower()
                                                 
                                                dna2 = list(dna2)
                                                
                                                # #replaces old codon overlapping PAM with new similar codon
                                                dna2[a.end() - 8], dna2[a.end() - 7], dna2[a.end() - 6] = l
                                                dna2[a.end() - 5], dna2[a.end() - 4], dna2[a.end() - 3] = i

                                #list becomes a string again for legibility 
                                dna2 = ''.join(dna2)
                                
                                # print out the new DNA string
                                if sys.argv[1].find('U') == -1:
                                    dna2 = dna2.replace('U', 'T')
                                    print dna2
                                    quit()
                                else:
                                    print dna2   
                                    quit()
    
        else:
            print 'Target not adjacent to a PAM site'
else:
    print 'Guide does not match or is not unique'
