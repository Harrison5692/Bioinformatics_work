#Overlap test
import sys

with open('Toy_bed.csv','r') as bed_table1, open('Toy_bed2.csv','r') as bed_table2:
    for line in bed_table1:

        chrm_1 = line.split('\t')[0]
        strand_1 = line.split('\t')[5]
        Pstart_1 = line.split('\t')[1]
        Pstop_1 = line.split('\t')[2]

        
        for line in bed_table2:
            chrm_2 = line.split('\t')[0]
            strand_2 = line.split('\t')[5]
            Pstart_2 = line.split('\t')[1]
            Pstop_2 = line.split('\t')[2]

            
            if chrm_2 == chrm_1 and strand_2 == strand_1:
                
                Pstart_1 = int(Pstart_1)
                Pstart_2 = int(Pstart_2)
                Pstop_1 = int(Pstop_1)
                Pstop_2 = int(Pstop_2)

                
                #Here the T1 peak starts first
                if Pstart_1 <= Pstart_2 and Pstart_2 <= Pstop_1 or Pstart_1 <= Pstart_2 and Pstart_2 == (Pstop_1+1):
                    
                    New_start = Pstart_1
                    
                    if Pstop_2 <= Pstop_1:
                    
                        New_stop = Pstop_1

                        break
                    
                    elif Pstop_2 >= Pstop_1:
                    
                        New_stop = Pstop_2

                        break       


                #Here the T2 peak starts first 
                elif Pstart_1 >= Pstart_2 and Pstop_2 >= Pstart_1 or Pstart_1 >= Pstart_2 and Pstop_2 == (Pstart_1-1):

                    New_start = Pstart_2

                    if Pstop_2 <= Pstop_1:

                        New_stop = Pstop_1

                        break

                    elif Pstop_2 >= Pstop_1:

                        New_stop = Pstop_2

                        break

                                             
                #no overlap and the T1 peak is first
                if Pstop_1 < (1+Pstart_2):
                    New_start = Pstart_1
                    New_stop = Pstop_1
                    
                    New_start2 = Pstart_2
                    New_stop2 = Pstop_2

                    break

                #no overlap and the T2 peak is first
                elif Pstop_2 < (Pstart_1-1):
                    New_start = Pstart_2
                    New_stop = Pstop_2
                    
                    New_start2 = Pstart_1 
                    New_stop2 = Pstop_1
                    break




        #in the event of no overlap we want to write both original beds                            
        if int(Pstop_1) < (1+int(Pstart_2)) or int(Pstop_2) < (int(Pstart_1)-1):
            list_of_peaks = []
            list_of_peaks.append(New_start)
            list_of_peaks.append(New_stop)
            list_of_peaks.append(line.split('\t')[0])
            list_of_peaks.append(line.split('\t')[5])

            list_of_peaks2 = []
            list_of_peaks2.append(New_start2)
            list_of_peaks2.append(New_stop2)
            list_of_peaks2.append(line.split('\t')[0])
            list_of_peaks2.append(line.split('\t')[5])

            #if multiple peaks merge then...
            # if New_start <= Pstart_1 and New_stop <= Pstop_1:
            #     New_start = New_start
            #     New_stop = Pstop_1
            #     pass
            # elif New_start <= Pstart_2 and New_stop <= Pstop_2:
            #     New_start = New_start
            #     New_stop = Pstop_2
            #     pass

            print(list_of_peaks)
            print(list_of_peaks2)



        #in the event of overlap again we want the new peak
        else: 
                

            list_of_peaks = []
            list_of_peaks.append(New_start)
            list_of_peaks.append(New_stop)
            list_of_peaks.append(line.split('\t')[0])
            list_of_peaks.append(line.split('\t')[5])  




            print(list_of_peaks)


     

       

'''
#Sample information of two Chromosomes: name, start, and stop location.
dna1 = ['chr1', 10, 30]
dna2 = ['chr1', 18, 19]

#If sequence is similar and any overlap exists, state whether this is overlap
if dna1[0] == dna2[0]:
	if ((dna1[1] >= dna2[1] and dna1[1] <= dna2[2]) or (dna2[1] >= dna1[1] and dna2[1] <= dna1[2])) : #Checks if dna1 start is in range/touching second sequence and vice versa for dna2 start
		print ("overlap exists") 
	else:
		print('regions do not overlap')
else:
	print("Regions are not on the same chromosome")    #If the starting chromosomes arent the same the program doesn't bother checking for overlap. 
'''