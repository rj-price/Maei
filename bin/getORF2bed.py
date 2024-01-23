#Jamie Pike - 14/12/2022
#GetORF to bed format. 

###############################################
#Set up and import correct modules.
from Bio import SeqIO
import pandas as pd
import numpy as np
import re, sys, os

###############################################
#Establish input file (BED).
infile = sys.argv[1]
outfile = sys.argv[2]

#PREPARE EMPTY DATA FRAME AND ESTABLISH VARIABLES. 
###############################################
Bed = { 
   'Chrom':[],
   'ChromStart':[],
   'ChromEnd':[],
   'sense':[],
} #Create a Dataframe which follows the bed format, using the output from getORF FASTA headers and convert to Pandas dataframe. 
df = pd.DataFrame(Bed)

#PREPARE FASTA HEADERS
###############################################
for record in SeqIO.parse(infile, "fasta"): #Loop through the input FASTA file headers 
   record = record.description #Replace the double underscore for ",_" in the contig names.
   FastaHeaders = re.split(" \[|\]| \- | \(*\)",record)  #Separate the contig name, start, stop and sense locations in FASTA headers.
   df.loc[len(df)] = FastaHeaders #Add the lists created by splitting the FASTA headers to the dataframe. 

#INSERT SENSE INFORMATION
###############################################
df['strand'] = np.where(df['sense']!= " (REVERSE SENSE)", "-", "+") #Add a column where the sense location reported by getORF is printed as a "+" or "-" 

#MOVE GETORF UID TO NAME COLUMN
###############################################
df['Names'] = df.loc[:, 'Chrom'] #Duplicate the chrom the column, which contains the getORF uid. 
column_to_move = df.pop("Names") #Set variable for the column to move.
df.insert(4, "Names", column_to_move) #Move the Names column which contains getORF UIDs to the NAMES column position in the bedfile.

#TIDY CONTIG NAMES AND COLUMNS TO FIT BED FORMAT
###############################################
df['Chrom'] = df['Chrom'].str.split('_').str[:-1].str.join('_') #Remove the unique identified from getORF which appears at the end of each contig name.  
df.insert(5, "Score", ".") #Insert Score column which is standard in bedfiles with strand information, and infill rows with a "." as we do not have this data.

#SWAP THE START AND END POS IF END IS > START
###############################################
m = df['ChromStart'] > df['ChromEnd']
df.loc[m, ['ChromStart', 'ChromEnd']] = (
    df.loc[m, ['ChromEnd', 'ChromStart']].values)


col = "sense" #Set a vraibale for the sense column. 
df1 = df.loc[:, df.columns != col] #Remove the sense column created, as we now have the stand column to indicate direction.

#WRITE OUT TO BED FILE.
###############################################
df1.to_csv(outfile, header=False, index=False, sep="\t") #Save the output to a tab seprated file, removing the row index and column headers. 
