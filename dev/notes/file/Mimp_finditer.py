#!/bin/python

#############################################
#### SEARCHES FOR ALL MIMP IN SEQUENCE ####
#############################################
from datetime import datetime
startTime = datetime.now()
import re, sys

from Bio import SeqIO
# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

############################################

infile = sys.argv[1]

print(f'\n\nRunning Mimp_finder.py on {infile} at {startTime}.')


with open(f'{infile}_mimp_hits.bed'.format(), "a") as file_object:
    for sequence in SeqIO.parse(infile, "fasta"):
        print("="*80 + "\n")
        print(f'Searching {sequence.description} for mimps.\n')
        hit = re.finditer(r"CAGTGGG..GCAA[TA]AA", str(sequence.seq), re.IGNORECASE)
        mimp_length = 400
        mimp_count = 0
        for h in hit:
            if h:
                h_start = h.start()
                print(f'--Mimp TIR found: "{h.group()}" at position {h.start()} to {h.end()}--')
                hit_rc = re.finditer(r"TT[TA]TTGC..CCCACTG", str(sequence.seq), re.IGNORECASE)
                for h_rc in hit_rc:
                    print('Looking for reverse complement sequence(s)...')
                    h_rc_end = h_rc.end()
                    if h_rc:
                        print(f'Mimp reverse complement TIR found: "{h_rc.group()}" at position {h_rc.start()} to {h_rc.end()}.  Is this complete mimp?')
                        length = h_rc_end - h_start
                        if length > 0:
                            if length < mimp_length:
                                print(f'Yes! Full mimp found at position {h.start()} to {h_rc.end()}.')
                                print("The mimp length is:")
                                print(f'{length}bp\n')
                                mimp_count = mimp_count+1
                                print(sequence.description, h.start(), h_rc.end(), file=file_object, sep="\t")
                            else:
                                print("No, TIRs are at a greater distance than 400bp. \nOnly a partial mimp.\n")
                        else:
                            print("No, the TIR's reverse complement is from the previous mimp as it comes before the start location of this TIR\n")
                print(f'\n{sequence.description} contains {mimp_count} mimp(s)')

    print("="*80 + "\n")

    print(f'\nMimp_finder.py complete. \n\nFile searched: {infile}\n\nFinished at: {startTime}')
