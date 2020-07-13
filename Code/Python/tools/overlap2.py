from collections import defaultdict

region_dict = defaultdict(list)
#new_regions = []

#Opens the toy beds that have been created
with open('Toy_bed.csv','r') as bed_table1, open('Toy_bed2.csv','r') as bed_table2:


    #for each line in both of the table give variables to each column
    for line in bed_table1:
        
        #this lets us validate looking at the correct files
        print(line)

        #labels for each column item
        chrs = line.split('\t')[0]
        Pstart = int(line.split('\t')[1])
        Pstop = int(line.split('\t')[2])
        strand = line.split('\t')[5]
        region_dict[chrs].append((Pstart, Pstop))

    #do the same for the other bed file
    for line in bed_table2:
        print(line)
        chrs = line.split('\t')[0]
        Pstart = int(line.split('\t')[1])
        Pstop = int(line.split('\t')[2])
        strand = line.split('\t')[5]
        region_dict[chrs].append((Pstart, Pstop))



new_regions = []

# For all the chromsomes in the dictionary we created. In this case its only the one but this
# is what would need to be done for multple keys (multple chromosomes)
for chrs in region_dict:
    
    
    new_stop = 0
    #sorting will order the regions making them easier to work with
    region_dict[chrs].sort()
    while len(region_dict[chrs]) > 1:


        region_1 = region_dict[chrs][0]
        start_1 = region_1[0]
        stop_1 = region_1[1]

        region_2 = region_dict[chrs][1]
        start_2 = region_2[0]
        stop_2 = region_2[1]

    

        if stop_1 >= start_2:

            start = start_1
            stop = stop_2
            region_new = (start, stop)

            region_dict[chrs][0] = region_new
            
            del region_dict[chrs][1]

            print(region_dict)

            #new_regions.append(region_new)
        else:
            
            new_regions.append(region_1)
            del region_dict[chrs][0]
            
    if len(region_dict[chrs]) == 1:

        new_regions.append(region_dict[chrs][0])

            

    print(region_dict)
    print(new_regions)





#------------------------------------------------------------------------

   #remove the first region from the original file and loop back to see if this subsequent
    #peak merges with the next
        # else:

        #     if stop == new_stop:
        #         del region_dict[chrs][0]
        #     else:    
        #         new_regions.append(region_1)
        #         del region_dict[chrs][0]

             

        # print(region_dict)
        # print(new_regions)







	#Start with empty dictionary
	#I want this dictonary to have one key for each chromosome number
	#Each key(chromosome) can have multiple regions

	#for every line in the file that has the same chromosome I want to
	#get their respective regions and append them as the value for their respective key



	#If chromsomes are the same in both files then write one Key for that chromosome to a dictionary
	# if chrs_b1 == chrs_b2:

	# 	dictionary = {chrs_b1:[Pstart_b1,Pstop_b1]}
	# 	print(dictionary)
#-------------------------------------------


####heres where jost is helping me 

        # i = 0
        # region_1 = region_dict[chrs][i]
        # start_1 = region_1[0]
        # stop_1 = region_1[1]
        # merge = True
        # while merge:
        #     print(len(region_dict[chrs]))
        #     print(i)
        #     if len(region_dict[chrs]) > i:
        #         region_2 = region_dict[chrs][i+1]
        #         start_2 = region_2[0]
        #         stop_2 = region_2[1]
        #     else:
        #         region_new = (start_1, stop_1)
        #     #Check if the first list region overlaps with the second list region
        #     if start_2 <= stop_1:

        #         #this is the new region
        #         region_new = (start_1, stop_2)

        #         #if so, then write the new merged coordinates
        #         stop_1 = stop_2
        #         i += 1
        #         merge = True

        #     else:
        #         region_new = (start_1, stop_1)
        #         merge = False

        # #print region_dict[chrs]
        # new_regions.append(region_new)
        # del region_dict[chrs][:i+1]
        # #print region_dict[chrs]

        # print new_regions