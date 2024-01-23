# Effector Set Analysis

Using the set of effectors from the run of Maei on 15/02/2021. 

A non-redundant candidate effector set was generated using MAFFT. 

The Script Effector_Seq_blastback.sh was written and then used to BLAST the set of candidate effectors (All_effector_sets_combined.non_redundant.fsa_aa)  against all the genomes, and determine presence or absence. The output can be found on Vettel here:

/home/u1983390/Fusarium_data/EFFECTOR_SET_ANALYSIS/EFFECTOR_SEQ_BLASTBACK_10-03-2021

The Effector_seq_blastback.sh command: 

```
./Effector_seq_blastback.sh -b -q All_effector_sets_combined.non_redundant.fsa_aa -g Fasta_list.txt
```

The Effector_seq_blastback.sh code:

```
#!/bin/bash

#Jamie Pike
#Using TBLASTN, from a list of genomes identify which genomes the in which the sequence can be found.

function main(){
  check $@

python -c "print('=' * 80)" ;
#Set Variables

date=$(date '+%d-%m-%Y')
echo "Date: $date"
query_filename="${query##*/}"
query_extension="${query_filename##*.}"
query_filename="${query_filename%.*}"

red='\033[01;31m'
none='\033[0m'

#Indicate the bedfile option was selected.

  if [ $b ]
  then
    echo "Bed file option selected. A bed file will be generated for each hit which meet the threshold (>90% identity and an evalue of < 1e-6)." ;
  else
    echo "No bed files being generated. Please use the -b flag if you would like a bed file for each genome (sequences included in bedfile will have >90% identity and an evalue of < 1e-6)."
  fi

#Indicate the genomes being searched
  echo "Databases searched:"
  for i in $(cat ${genomes}) ;
    do echo "${i}" ;
  done
  python -c "print('=' * 80)" ;

#Make Blast Databases.
  #For the query protiens
  echo "Making blast database for ${query}..." ;
  makeblastdb -dbtype prot -in ${query} ;
  python -c "print('=' * 80)" ;

  for i in $(cat ${genomes}) ;
    do
    echo "Making blast database for ${i}..." ;
    makeblastdb -dbtype nucl -in ${i} ;
    python -c "print('=' * 80)" ;
    done

    #Run tblastn
      mkdir EFFECTOR_SEQ_BLASTBACK_${date}

  if [ $b ]
  then
      for i in $(cat ${genomes}) ;
        do
        python -c "print('=' * 80)" ;
        subject_filename="${i##*/}" ;
        subject_extension="${subject_filename##*.}" ;
        subject_filename="${subject_filename%.*}" ;
        mkdir ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename} ;
        mkdir ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/BLAST_DB ;
        echo "Running TBLASTN for ${i}..." ;
        echo "query: ${query}"
        echo "db: ${i}"
        echo "out: ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out"
        tblastn -query ${query} -db ${i}  -outfmt 6 -evalue 1e-6  -out ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out ;
        if [[ -s ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out ]] ;
        then
          awk '{ if ($3 >= 90) print $0}' ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out | awk '{ if ($11 < 1e-6) print $1,"\t",$2}' > ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out ;
          cut -f 1,2,9,10 ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out | awk '{print $2,$3,$4,$1}' OFS='\t' | awk '{if ($2>$3)print $1,$3,$2,$4,".","-";else print $1,$2,$3,$4,".","+";}' OFS='\t' |awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.all_hits.bed ;
          awk '{ if ($3 >= 90) print $0}' OFS='\t' ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out | awk '{ if ($11 < 1e-6) print $0}' OFS='\t' | cut -f 1,2,9,10 | awk '{print $2,$3,$4,$1}' OFS='\t' | awk '{if ($2>$3)print $1,$3,$2,$4,".","-";else print $1,$2,$3,$4,".","+";}' OFS='\t' |awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.hits_within_threshold.bed
          echo " "
          echo -e "\e[1;31mQuery hits in ${i}: \e[0m" ;
          echo -e "\e[1;31mQUERY SEQUENCE  NUMBER OF HITS \e[0m" ;
          echo -e "\e[1;31m--------------  -------------- \e[0m" ;
          awk -v r="$red" -v n="$none" '{a[$1]++}END{for(i in a){printf "%s%s\t%s%s\n", r, i, a[i], n}}' ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out | column -t
          awk '{a[$1]++}END{for(i in a){print i,"\t",a[i]}}' ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out > ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/Occurances_of_${query_filename}_sequences_in_${subject_filename}.txt
          echo "\nBed file generated: ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.bed"
          cp ${i} ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${i} ;
          mv ${i}.n* ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/BLAST_DB ;
        else
          echo "
          No hits found.
          "
          echo "##No hits found." >> ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out
        fi
        done
    else
      for i in $(cat ${genomes}) ;
        do
        python -c "print('=' * 80)" ;
        subject_filename="${i##*/}" ;
        subject_extension="${subject_filename##*.}" ;
        subject_filename="${subject_filename%.*}" ;
        mkdir ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename} ;
        mkdir ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/BLAST_DB ;
        echo "Running TBLASTN for ${i}..." ;
        echo "query: ${query}"
        echo "db: ${i}"
        echo "out: ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out"
        tblastn -query ${query} -db ${i}  -outfmt 6 -evalue 1e-6  -out ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out ;
        awk '{ if ($3 >= 90) print $0}' ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.tblastn.outfmt6.evalue_used_1e-6.out | awk '{ if ($11 < 1e-6) print $1,"\t",$2}' > ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out ;
        echo " "
        echo -e "\e[1;31mQuery hits in ${i}: \e[0m" ;
        echo -e "\e[1;31mQUERY SEQUENCE  NUMBER OF HITS \e[0m" ;
        echo -e "\e[1;31m--------------  -------------- \e[0m" ;
        awk -v r="$red" -v n="$none" '{a[$1]++}END{for(i in a){printf "%s%s\t%s%s\n", r, i, a[i], n}}' ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out | column -t
        awk '{a[$1]++}END{for(i in a){print i,"\t",a[i]}}' ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out > ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/Occurances_of_${query_filename}_sequences_in_${subject_filename}.txt
        cp ${i} ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/${i} ;
        mv ${i}.n* ./EFFECTOR_SEQ_BLASTBACK_${date}/${subject_filename}/BLAST_DB ;
        done
      fi

  mkdir ./EFFECTOR_SEQ_BLASTBACK_${date}/$query_filename
  cp $query ./EFFECTOR_SEQ_BLASTBACK_${date}/$query_filename
  mv $query.p* ./EFFECTOR_SEQ_BLASTBACK_${date}/$query_filename
}

function help(){
  python -c "print('=' * 80)"
  echo "
  ##############################################################################
  Effector_seq_blastback.sh

  Usage: cmd [OPTIONAL ARGUMENTS] -q <FASTA> -g <TEXT FILE>

  -q  Provide a FASTA file as input which contains a set of query sequences (sequences you're searching for).
  -g  Provide a TEXT file listing the genomes to search through.

  OPTIONAL ARGUMENT:

  -b  Generate a bedfile from each hit with >90% identity and an evalue of < 1e-6.

  Note the genome files should be in the directory you wish to run the script or the text files should contain the full path to the genome files.

  DEPENDENCIES
  BLAST: Protein Query-Translated Subject BLAST 2.2.31+
  BEDTOOLS: v2.25.0
  ##############################################################################
  "
  python -c "print('=' * 80)"
}

function help2(){
  echo "
  ##############################################################################
  Effector_seq_blastback.sh

  Usage: cmd [OPTIONAL ARGUMENTS] -q <FASTA> -g <TEXT FILE>

  -q  Provide a FASTA file as input which contains a set of query sequences (sequences you're searching for).
  -g  Provide a TEXT file listing the genomes to search through.

  OPTIONAL ARGUMENT:

  -b  Generate a bedfile from each hit with >90% identity and an evalue of < 1e-6.

  Note the genome files should be in the directory you wish to run the script or the text files should contain the full path to the genome files.

  DEPENDENCIES
  BLAST: Protein Query-Translated Subject BLAST 2.2.31+
  BEDTOOLS: v2.25.0
  ##############################################################################
  "
  python -c "print('=' * 80)"
}

function check(){
#  local OPTIND opt i
  while getopts ":q:g:b" opt; do
    case $opt in
      q) query="$OPTARG" ;;
      g) genomes="$OPTARG" ;;
      b) b=true ;;
      \?) help; exit 1
      ;;
    esac
  done
#  shift $((OPTINF -1))

  if [ "$query" = "" ]
  then
    python -c "print('=' * 80)"
    echo "
    //Please provide a FASTA file containing query sequences." ;
    help2 ;
    exit 2
  fi

  if [ "$genomes" = "" ]
  then
    python -c "print('=' * 80)"
    echo "
    //Please provide a genome to be searched." ;
    help2 ;
    exit 3
  fi
}

main $@

#  echo "Invalid option: $OPTARG requires an argument" ; exit 1
```

Effector presence or absence was then identified for each genome in an excel document, and an overall effector repertoire for races was determined. 

[List_of_effectors_and_the_genomes_they_belong_to.xlsx](./file/List_of_effectors_and_the_genomes_they_belong_to.xlsx)

![image.png](image/image.png)
