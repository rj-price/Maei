#!/bin/python

#### ORF Finder ####
#Jamie Pike 26.01.2021

#Search through a FASTA file from getorf and identify all open reading frames within the fasta.

#Set up and import correct modules.

from datetime import datetime
import re, sys

from Bio import SeqIO
# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

############################################
#Set colours

CRED = '\033[91m'
CEND = '\033[0m'

############################################
#Establish input file (FASTA) and start time.

startTime = datetime.now()

if len(sys.argv) <= 1: #if there is no fasta file provided, print this error message and exit.
    print('\nERROR: No input fasta provided.\nPlease follow the command with a fasta file, e.g. <python ./find_all_ORF_from_getorf.py Input_file.fsa_aa>.\n')
    exit(1)
else:
    infile = sys.argv[1]  #takes a getorf fasat file as input.

#############################################
#Create an output fasta file.

    name = infile[:-7] #removes ".fsa_aa" at the end of the input fatsa.

    outfile = f'{name}_all_ORF.fsa_aa'  #to rename the output file.

    print(f'\n\nRunning find_all_ORF_from_getorf.py on {infile} at {startTime}.') #Print statement to indicate the script has started.

    #Open a .fsa_aa file using the name of the infile for the orfs to be deposited with.

    with open(outfile, "a") as file_object:

#############################################
#Search for smaller ORFs and deposit them in the new fasta file.

    #Parse and loop through the FASTA input file. Search through each scaffold/contig of the FASTA looking for each methionine.
    #For each methonine write that methonine to the end of the sequence into the new output fasta file.

        for sequence in SeqIO.parse(infile, "fasta"): #Loop through each sequence with the Fasta.
            print("="*80 + "\n")
            print(f'Searching {sequence.description} for shorter ORFs.\n') #print the description of the sequence that is being searched.
            hit = re.finditer(r"M", str(sequence.seq), re.IGNORECASE) #find all "M" (methonine) within the sequence.
            M_pos = 1
            for h in hit: #for each instance there is an "M"
                if h:
                    h_start = h.start() #record the location of the "M"
                    sequence_length = len(sequence.seq) - h_start #How long is the sequence from the current "M"
                    if h_start == 0: #if h_start = 0 then this is the origanl sequence produced by getorf and will need to be included in the output fasta.
                        print(f'--Original ORF: "{h.group()}" at position {h.start()}') #print the location of the "M" to the screen.
                        print(CRED + f'{sequence.seq[h_start:]}' + CEND ) #then print the remaining sequence to the screen.
                        file_object.write(f'>{sequence.description}_ORF_from_position_{h.start()},\n{sequence.seq[h_start:]}\n') #print original sequence found by getorf in the output fatsa.
                    else:
                        if sequence_length > 30:  #if the length of the sequence after M is greater than 30.
                            if h_start > 1: #if the location of the "M" is greater that one (not the first M in the sequence, the one getorf identified as the start of the ORF)
                                print(f'--Smaller ORF found >30aa: "{h.group()}" at position {h.start()}') #print the location of the "M" to the screen.
                                print(CRED + f'{sequence.seq[h_start:]}' + CEND ) #then print the remaining sequence to the screen.
                                file_object.write(f'>{sequence.description}_ORF_from_position_{h.start()},\n{sequence.seq[h_start:]}\n') #Within the fasta created using the "with" statement, write the ">"followed by the sequence description and the location of the "M". Enter a new line, print the sequence from the current "M" identified, then print a new line.

###########################################
#End.
        print("="*80 + "\n")

        print(f'\nfind_all_ORF_from_getorf.py complete. \n\nFile searched: {infile}\n\nFinished at: {startTime}') #print statement to indicate that the search has ended.
