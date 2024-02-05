#!/bin/bash

#Jamie Pike - Adapted 04/11/2022
#Version 5. Overhaul of pipeline structure, generating GTFs, masking non-mimp regions and performing effector analysis (EffectorP) before clustering.
#Identify MIMP sequences in a list of Foc genomes, and prepare a fasta of these hits, as well as table which lists the Subject ID, Query ID, location, and sequence.
#Note, the -b flag can be altered in bedtools to produce a fasta file which is inlcludes hit plus x many bases either side. This will alter the resulting location_and_seq table however, as the hit plus the number of bases either side will be prinited in location and sequence section.
#all genomes must be in the same directory unless Paths are created.


#PATHS TO MAEI ASSOCIATED SCRIPTS. 
ProcessingCDHIT=/Volumes/Jamie_EXT/Projects/Maei/bin/Processingcdhit.py

##############################
###PLEASE DO NOT EDIT BELOW###
##############################

#Insert name of file containing a list of genomes
Cl=${1?Error: No sequence identity threshold set. \n Please provide a sequence identity threshold for cd-hit, e.g. 0.8.
 	Calculated as:
 	Number of identical amino acids in alignment divided by the full length of the shorter sequence.}

#Print Summary to the screen
python -c "print('=' * 75)"
echo "Mimp-Associated Effector Clustering"
echo "-----------------------------------"
echo $(date)
echo "usage: ./Maei_v5.sh <List_of_Genome_FASTAs> <Sequence identity threshold for cd-hit, e.g. 0.8>"
python -c "print('=' * 75)"

echo "Clustering final effector set..."
mkdir cdhit
cd-hit -i ./AllCandidateEffectors-AfterBLASTandFilter.fasta -d 0 -o ./AllCandidateEffectorSet-AfterBLASTandFilter_cd-hit.out.${Cl} -c ${Cl} -n 5  -G 1 -g 0 -b 20 -l 10 -s 0.0 -aL 0.0 -aS 0.0 1> ./cdhit/cd-hit.log #Cluster the combined effector set based on sequnce identitiy.
mv ./AllCandidateEffectorSet-AfterBLASTandFilter_cd-hit.out.${Cl} ./cdhit/./AllCandidateEffectorSet_CulsterRepresentativeSeqs.fasta
echo "Processing Cd-Hit output..."
python ${ProcessingCDHIT}  ./AllCandidateEffectorSet-AfterBLASTandFilter_cd-hit.out.${Cl}.clstr 
echo "Done."

python -c "print('=' * 75)"
echo "Mimp-Associated Effector Identification Complete"
echo "------------------------------------------------"
echo $(date)
python -c "print('=' * 75)"

