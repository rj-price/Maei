#!/usr/bin/env python3

#Jamie Pike - 14/12/2022
#Covert BED file from Maei pipeline to GFF 

###########################################################
#Set up and import correct modules.
import re, sys, os

###############################################
#Establish input file (BED).
infile = sys.argv[1]
type = sys.argv[2]
SeqID =sys.argv[3]

###############################################
#Set count variable
count = 0 

#Loop through the BED file by line.
for line in open(infile):
   count = count + 1 
   #Split the BED file into columns using  \t
   Column = line.strip().split("\t")
   #Add a value of 1 to the start of the candiate effector location as GFFs the BED format developed at UCSC uses a zero based indexing and an open end interval whereas the GFF format developed at Sanger assumes a 1 based coordinate system that includes both start and end coordinates.
   Start = str(int(Column[1]) + 1)
   #Trim the candidate effector name by removing the information following the bracket. 
   Separator = '('
   TrimName = Column[3]
   #Print the output to screen (tab separated).  
   print(Column[0], "bed2gff", type, Start, Column[2], Column[4], Column[5], ".", "ID=%s;Name=%s;" % (SeqID, (TrimName)), sep = '\t')