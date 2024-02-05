#!/bin/bash 
#Filter BLAST Hit results. 

set -e 


###################
gL=${1?Error: no list of FASTAs given. \n Please provide a .txt file containing a list of FASTAs for analysis}

python -c "print('=' * 75)"
echo "Filtering BigBlast results for effector candidates"
echo "--------------------------------------------------"
echo $(date)
python -c "print('=' * 75)"


for i in $(cat ${gL}); do
#SignalP Analysis of each BLAST Hit Set
    echo "Running for: ${i}"
    echo "---"
    echo "Merging hits in the same location..."
    sed 's/:/_/' ./${i}/*hits_within_threshold.bed > ./${i}/${i}-hits_within_threshold.bed ; #Ensure that headers are not trimmed by transeq due to species characters. 
    bedtools sort -i ./${i}/${i}-hits_within_threshold.bed | bedtools merge -header -d 10 -c 4,6 -o distinct  > ./${i}/${i}-mergedhits.bed ; #Merge Hits in the same location
    bed2fasta.sh -b ./${i}/${i}-mergedhits.bed -f ./${i}/${i}.fna -e 0 ; #Extract the hits from the genome for processing
    echo "Running Transeq..."
    transeq -sformat pearson -frame 6 -sequence ./${i}/${i}-mergedhits.bed.expanded_by_0.fasta -outseq  ./${i}/${i}-Translated.fasta  #Translate to protein sequnce. 
    echo "Running SignalP for ${i}..." ;
    signalp -f summary ./${i}/${i}-Translated.fasta 1>./${i}/${i}-SignalP.out 2>./${i}/${i}-SignalP.log #Identify sequences with signal peptide.
    echo "Generating index files..." ;
    grep "^Name*" ./${i}/${i}-SignalP.out | awk '$2 ~ /YES/ {print $0}' > ./${i}/${i}-SignalP.filtered.out ;
    awk  '{gsub("=","\t",$0); print;}' ./${i}/${i}-SignalP.filtered.out | awk '{print $2}' > ./${i}/${i}-SignalP.filtered.list.txt ; #Filter by presence of signal peptide.
    for j in $(cat ./${i}/${i}-SignalP.filtered.list.txt); 
        do samtools faidx ./${i}/${i}-Translated.fasta ${j} ; 
    done > ./${i}/${i}-SignalP.filtered.fasta ; 
    echo "Running EffectorP..." ;
    EffectorP.py -i ./${i}/${i}-SignalP.filtered.fasta -E ./${i}/${i}-EffectorP.filtered.fasta > ./${i}/${i}-EffectorP.filtered.log ; 
    echo "Generating output bed files..." 
    grep ">"  ./${i}/${i}-EffectorP.filtered.fasta  | sed 's/>//g' | sed 's/\:0.*//' | sed 's/\.*t1//' | sed 's/ .*//' |  #Create an effector set list.
    awk '{gsub("::","\t",$1); gsub(":","\t",$1) ; print $0}' | 
    awk -v OFS="\t" '{ gsub("-","\t",$3); gsub(/\(/,"\t",$3) ; gsub(/\)/,"\t",$3); gsub("_","",$3);  print $0}'| 
    awk  -v OFS="\t" '{print $2, $3, $4, $1, $5, $6}'it | bedtools sort -i > ./${i}/${i}-EffectorP.filtered.bed
    echo "Done."
    echo "----------------------"
done 

python -c "print('=' * 75)"


