#!/bin/bash

#Jamie Pike - Adapted 04/11/2022
#Version 4 to include aa size filtering and cluster at 80% identity, not 90.
#Identify MIMP sequences in a list of Foc genomes, and prepare a fasta of these hits, as well as table which lists the Subject ID, Query ID, location, and sequence.
#Note, the -b flag can be altered in bedtools to produce a fasta file which is inlcludes thit plus x many bases either side. This will alter the resulting location_and_seq table however, as the hit plus the number of bases either side will be prinited in location and sequence section.
#all genomes must be in the same directory unless Paths are created.


#############################
###PLEASE EDIT BELOW PATHS###
#############################

#PLEASE SET PATHS TO THE FOLLOWING DEPENDENCIES:
Run_Bio_python=/home/u1983390/miniconda3/envs/MaeiEnv/bin/python3           #path to biopython
Mimp_finditer=/home/u1983390/Fusarium_data/MimpAssociatedEffectorIdent/Scripts/Mimp_finditer.py                     #path to mimp_finditer.py script
Run_Augustus=/home/u1983390/miniconda3/envs/MaeiEnv/bin/augustus             #path to Augustus executable
Run_getAnnoFasta=/home/u1983390/miniconda3/envs/MaeiEnv/bin/getAnnoFasta.pl  #path to Augustus getAnnoFasta.pl
Emboss_get_orf=/home/u1983390/miniconda3/bin/getorf               #path to EMBOSS getorf executable
Run_SignalP=/home/u1983390/apps/FoEC/signalp-4.1/signalp                         #path to SignalP executable
Run_EffectorP=/home/u1983390/apps/EffectorP-2.0-2.0.1/Scripts/EffectorP.py       #path to EffectorP.py
Run_find_small_ORFs=//home/u1983390/Fusarium_data/MimpAssociatedEffectorIdent/Scripts/find_all_ORF_from_getorf.py           #Path to find_all_ORF_from_getorf.py script
Run_nhmmer=/home/u1983390/miniconda3/envs/hmmerEnv/bin/nhmmer                    #path to HMMER


###############################
###PLEASE DO NOT EDIT BELOW ###
###############################

#Insert name of file containing a list of genomes
gL=${1?Error: no list of FASTAs given. \n Please provide a .txt file containing a list of FASTAs for analysis}


#Loop to identify mimps in each genome
for i in $(cat ${gL}); do

  #Create an index file for bedtools
  samtools faidx ${i}; cut -f 1,2 ${i}.fai > ${i}.chrom.sizes ;

  #Run the HMMER mimp finder
  ${Run_nhmmer} --tblout "${i}_mimp_hits" --dna realn_All.aln.trimmed.hmm ${i} ;
  awk 'BEGIN { OFS="\t" }  {print $1,$7,$8,$12}' "${i}_mimp_hits" > ${i}_mimp_hmm_hits.bed
  #Sort the output from hmmer search
  awk '{if ($2>$3)print $1,$3,$2,".",".","-";else print $1,$2,$3,".",".","+";}' OFS='\t' ${i}_mimp_hmm_hits.bed | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t' >> ${i}_mimp_hits.bed.ordered

  #Run the mimp_finder script.
  ${Run_Bio_python} ${Mimp_finditer} ${i} ;
  #Sort the output from mimp_finditer.py
  awk '{if ($2>$3)print $1,$3,$2,".",".","-";else print $1,$2,$3,".",".","+";}' OFS='\t' ${i}_mimp_hits.bed | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t' >> ${i}_mimp_hits.bed.ordered;

  #Sort the bed files.
  bedtools sort -i ${i}_mimp_hits.bed.ordered > ${i}_mimp_hits.sorted.bed ;
  #Create a fasta containing the mimp sequence
  bedtools slop -i ${i}_mimp_hits.sorted.bed -g ${i}.chrom.sizes -b 0 > ${i}_mimp_hits.sorted.0.bed;
  bedtools getfasta -s -fi ${i} -bed ${i}_mimp_hits.sorted.0.bed -fo ${i}_mimp_hits.0.fasta;
  #create a fasta by expanding 2.5kb upstream and downstream of the mimp bases.
  bedtools slop -i ${i}_mimp_hits.sorted.bed -g ${i}.chrom.sizes -b 2500 > ${i}_mimp_hits.sorted.2500.bed;
  bedtools getfasta -s -fi ${i} -bed ${i}_mimp_hits.sorted.2500.bed -fo ${i}_mimp_hits.2500.fasta;
done

python -c "print('=' * 130)"
python -c "print('=' * 130)"
python -c "print('\n')"

echo "Creating directories"
#taking all of the previous outputs, create a folder for each genome in gL placing the all of the ouputs for the genome into the folder.
for i in $(cat ${gL});
  do mkdir ${i}_MIMPS; mv ./*${i}* ./${i}_MIMPS/; done
echo "//Ignore the mv: rename error message"
#This does produce an error message saying that the folder cannot be moved into itself. This can be ignored as the folder should not be moved into itself

for i in $(cat ${gL});
  do
  mkdir ./${i}_MIMPS/${i}_Predictions;
  mkdir ./${i}_MIMPS/${i}_SignalP_Raw; done

python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "Mimp searching complete."
python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "Predicting effectors..."
python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "Running AUGUSTUS on mimps"
python -c "print('=' * 80)"

for i in $(cat ${gL});
  do
  echo "//Running predictions for ${i}";
  ${Run_Augustus} --species=fusarium ${i}_MIMPS/${i}_mimp_hits.2500.fasta > ${i}_MIMPS/${i}_mimp_hits.2500.augustus.gff;
  echo "//Creating fasta for ${i}";
  ${Run_getAnnoFasta} ${i}_MIMPS/${i}_mimp_hits.2500.augustus.gff; done

for i in $(cat ${gL});
  do
  sed "s/>/>${i}_Candidate_Sequence_/" ${i}_MIMPS/${i}_mimp_hits.2500.augustus.aa > ${i}_MIMPS/${i}_mimp_hits.2500.augustus.fsa_aa; done


python -c "print('=' * 80)"
echo "Running getorf on mimps"
python -c "print('=' * 80)"

for i in $(cat ${gL});
  do echo "//Running predictions for ${i}";
  ${Emboss_get_orf} -find 1 ${i}_MIMPS/${i}_mimp_hits.2500.fasta -outseq ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf ; done

for i in $(cat ${gL});
  do
  sed "s/>/>${i}_Candidate_Sequence_/" ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf > ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.fsa_aa ; ${Run_Bio_python} ${Run_find_small_ORFs} ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.fsa_aa ;  done


python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "Running SignalP on predictions..."
python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "AUGUSTUS predictions"
python -c "print('=' * 80)"

#AUGUSTUS SIGNALP ANALYSIS AND FILTERING

for i in $(cat ${gL});
  do echo "//Running SignalP for ${i} Augustus predictions";
  ${Run_SignalP} -f summary  ./${i}_MIMPS/${i}_mimp_hits.2500.augustus.fsa_aa > ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.augustus.SignalP_raw.out;
  echo "//Filtering SignalP results for ${i} Augustus predictions";
  grep "^Name*" ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.augustus.SignalP_raw.out | awk '$2 ~ /YES/ {print $0}' > ./${i}_MIMPS/${i}_mimp_hits.2500.augustus.SignalP.out ; done

for i in $(cat ${gL});
  do echo "//Generating fasta of Augustus predictions for ${i} which contain a signal peptide";
  awk  '{gsub("=","\t",$0); print;}' ${i}_MIMPS/${i}_mimp_hits.2500.augustus.SignalP.out | awk '{print $2}' > ${i}_MIMPS/${i}_SignalPep_sequences_augustus_list.txt;
    for j in $(cat ${i}_MIMPS/${i}_SignalPep_sequences_augustus_list.txt;); do samtools faidx ${i}_MIMPS/${i}_mimp_hits.2500.augustus.fsa_aa ${j}; done > ${i}_MIMPS/${i}_mimp_hits.2500.augustus.SignalP.fsa_aa; done

python -c "print('=' * 80)"
echo "Getorf predictions"
python -c "print('=' * 80)"
#EMBOSS getorf SIGNALP ANALYSIS AND FILTERING

for i in $(cat ${gL});
  do echo "//Running SignalP for ${i} getorf output";
  ${Run_SignalP} -f summary  ./${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf_all_ORF.fsa_aa > ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.emboss.all_ORF.SignalP_raw.out;
  echo "//Filtering SignalP results for ${i} getorf output";
  grep "^Name*" ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.emboss.all_ORF.SignalP_raw.out | awk '$2 ~ /YES/ {print $0}' > ./${i}_MIMPS/${i}_mimp_hits.2500.emboss.all_ORF.SignalP.out ; done

for i in $(cat ${gL});
  do echo "//Generating fasta of getorf output for ${i} which contain a signal peptide";
  awk  '{gsub("=","\t",$0); print;}' ${i}_MIMPS/${i}_mimp_hits.2500.emboss.all_ORF.SignalP.out | awk '{print $2}' > ${i}_MIMPS/${i}_SignalPep_sequences_emboss.orf_list.txt;
    for j in $(cat ${i}_MIMPS/${i}_SignalPep_sequences_emboss.orf_list.txt;); do samtools faidx ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf_all_ORF.fsa_aa ${j}; done > ${i}_MIMPS/${i}_mimp_hits.2500.emboss.all_ORF.SignalP.fsa_aa ; done

python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "Tidying up"
python -c "print('=' * 130)"
python -c "print('=' * 130)"

for i in $(cat ${gL});
  do
  mv ./${i}_MIMPS/${i}_SignalPep_* ./${i}_MIMPS/${i}_SignalP_Raw/ ;
  mv ./${i}_MIMPS/${i}*.SignalP.out ./${i}_MIMPS/${i}_SignalP_Raw/ ;
  mv ./${i}_MIMPS/${i}*.2500.augustus.aa ./${i}_MIMPS/${i}_Predictions/ ;
  mv ./${i}_MIMPS/${i}*.2500.emboss.orf ./${i}_MIMPS/${i}_Predictions/ ;
  mv ./${i}_MIMPS/${i}*.augustus.fsa_aa.fai  ./${i}_MIMPS/${i}_Predictions/ ;
  mv ./${i}_MIMPS/${i}*.emboss.orf.fsa_aa.fai ./${i}_MIMPS/${i}_Predictions/ ;
  mv ./${i}_MIMPS/${i}*.2500.augustus.gff ./${i}_MIMPS/${i}_Predictions/ ; done

python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "Creating non-redundant set of predicted protiens with signal peptide..."
python -c "print('=' * 130)"
python -c "print('=' * 130)"


#Create fasta of all the augustus and getorf results and place in one fasta file
for i in $(cat ${gL});
 do
   echo "//Creating fasta from the augutsus and getorf predictions for ${i}" ;
   cat ${i}_MIMPS/${i}*.SignalP.fsa_aa > ${i}_MIMPS/${i}.SignalP.fsa_aa ; done


#Create a non-redundant protein set from the grouped fasta created above.
 for i in $(cat ${gL});
  do
    echo "//Creating non-redundant fasta from the augutsus and getorf predictions for ${i}" ;
    cd-hit -i ${i}_MIMPS/${i}.SignalP.fsa_aa -d 0 -o ${i}_MIMPS/${i}.SignalP.non_redundant.fasta -c 1 -n 5  -G 1 -g 1 -b 20 -l 10 -s 0.0 -aL 1 -aS 0.0 ; done

#Create a master fasta file containing all of the candidate effectors from all of the genomes.
echo "//Creating master fasta file for candiate sequences "
 for i in $(cat ${gL});
  do
    #Filter all sequences by size. <300 aa accepted.
    samtools faidx ${i}_MIMPS/${i}.SignalP.non_redundant.fasta ;
    awk '{if($2 < 300) print $1 "\t0\t" $2 "\t" $1}' ${i}_MIMPS/${i}.SignalP.non_redundant.fasta.fai > ${i}_MIMPS/${i}.SignalP.non_redundant.filtered.bed
    #Sort the bed files.
    bedtools sort -i ${i}_MIMPS/${i}.SignalP.non_redundant.filtered.bed > ${i}_MIMPS/${i}.SignalP.non_redundant.filtered.sorted.bed ;
    #Create a fasta containing the mimp sequence
    bedtools slop -i ${i}_MIMPS/${i}.SignalP.non_redundant.filtered.sorted.bed -g ${i}_MIMPS/${i}.SignalP.non_redundant.fasta.fai -b 0 > ${i}_MIMPS/${i}.SignalP.non_redundant.sorted.0.bed;
    bedtools getfasta -s -fi ${i}_MIMPS/${i}.SignalP.non_redundant.fasta -bed ${i}_MIMPS/${i}.SignalP.non_redundant.sorted.0.bed -fo ${i}_MIMPS/${i}.SignalP.non_redundant.filtered.fasta ;
    cat ${i}_MIMPS/${i}.SignalP.non_redundant.filtered.fasta ; done  > ./all_mimp_associated_effectors_unclustered.fasta

python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo "Clustering the non-redundant protien set..."
python -c "print('=' * 130)"
python -c "print('=' * 130)"

#Cluster the combined effector set at 80% identity
cd-hit -i all_mimp_associated_effectors_unclustered.fasta -d 0 -o all_mimp_associated_effectors_CLUSTERED.fasta -c 0.8 -n 5  -G 1 -g 1 -b 20 -l 10 -s 0.0 -aL 0.0 -aS 0.0


python -c "print('=' * 130)"
python -c "print('=' * 130)"
echo 'Running EffectorP on all the  clustered candidate effectors...'
python -c "print('=' * 130)"
python -c "print('=' * 130)"


#run EffectorP for effector prediction on all filtered hits with a signal peptide and produce a fasta with the output.
python ${Run_EffectorP} -i all_mimp_associated_effectors_CLUSTERED.fasta -E all_mimp_associated_effectors_CLUSTERED.EffectorP.fasta

python -c "print('=' * 130)"
python -c "print('=' * 130)"
