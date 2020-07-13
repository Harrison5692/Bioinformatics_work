import sys
import pandas as pd

if sys.argv[1][0] == 't':

    xl = pd.ExcelFile('HelicaseCandidatesgRNAssequences2.xlsx')
    df = xl.parse('Sheet1')

    for i in df['sequence1']:
    	i = str(i)
    	i = i.lower()
    	sys.stdout = open("file1.txt", 'a+')   
	print i 
    quit()



   
