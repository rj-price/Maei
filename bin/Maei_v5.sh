#!/bin/bash

#Jamie Pike - Adapted 04/11/2022
#Version 5. Overhaul of pipeline structure, generating GTFs, masking non-mimp regions and performing effector analysis (EffectorP) before clustering.
#Identify MIMP sequences in a list of Foc genomes, and prepare a fasta of these hits, as well as table which lists the Subject ID, Query ID, location, and sequence.
#Note, the -b flag can be altered in bedtools to produce a fasta file which is inlcludes hit plus x many bases either side. This will alter the resulting location_and_seq table however, as the hit plus the number of bases either side will be prinited in location and sequence section.
#all genomes must be in the same directory unless Paths are created.


#############################
###PLEASE EDIT BELOW PATHS###
#############################

#PLEASE SET PATHS TO THE FOLLOWING DEPENDENCIES:
Run_Bio_python=/home/u1983390/miniconda3/envs/biopythonEnv/bin/python3 #path to biopython
Run_Augustus=/home/u1983390/miniconda3/envs/MaeiEnv/bin/augustus  #path to Augustus executable
Run_getAnnoFasta=/home/u1983390/miniconda3/envs/MaeiEnv/bin/getAnnoFasta.pl #path to Augustus getAnnoFasta.pl
Emboss_get_orf=/home/u1983390/miniconda3/bin/getorf #path to EMBOSS getorf executable
Run_SignalP=/home/u1983390/apps/FoEC/signalp-4.1/signalp  #path to SignalP executable
Run_EffectorP=/home/u1983390/apps/EffectorP-2.0-2.0.1/Scripts/EffectorP.py  #path to EffectorP.py
Run_nhmmer=/home/u1983390/miniconda3/envs/hmmerEnv/bin/nhmmer #path to HMMER
Run_AGAT_manage=/home/u1983390/miniconda3/bin/agat_sp_manage_attributes.pl
Run_AGAT_extract=/home/u1983390/miniconda3/bin/agat_sp_extract_sequences.pl #Path to AGAT for extarcting cds from mimp-associated gene/ORF gff. 
Run_cdhit=/usr/bin/cdhit #Path to cd-hit.

#PATHS TO MAEI ASSOCIATED SCRIPTS. 
Mimp_finditer=/home/u1983390/Fusarium_data/MimpAssociatedEffectorIdent/Fo._of_Lettuce_Analysis/Maei_v5_Scripts/bin/Mimp_finditer.py #path to mimp_finditer.py script
bed2gff=/home/u1983390/Fusarium_data/MimpAssociatedEffectorIdent/Fo._of_Lettuce_Analysis/Maei_v5_Scripts/bin/bed2gff.py #Bed to GFF file converter. 
getORF2bed=/home/u1983390/Fusarium_data/MimpAssociatedEffectorIdent/Fo._of_Lettuce_Analysis/Maei_v5_Scripts/bin/getORF2bed.py #Convert getORF output to bed file. 
ProcessingCDHIT=/home/u1983390/Fusarium_data/MimpAssociatedEffectorIdent/Fo._of_Lettuce_Analysis/Maei_v5_Scripts/bin/Processingcdhit.py #Script for processing the cd-hit output.

##############################
###PLEASE DO NOT EDIT BELOW###
##############################

#Insert name of file containing a list of genomes
gL=${1?Error: No list of FASTAs given. \n Please provide a .txt file containing a list of FASTAs for analysis}
Cl=${2?Error: No sequence identity threshold set. \n Please provide a sequence identity threshold for cd-hit, e.g. 0.8.
 	Calculated as:
 	Number of identical amino acids in alignment divided by the full length of the shorter sequence.}

#Print Summary to the screen
python -c "print('=' * 75)"
echo "Mimp-Associated Effector Identification"
echo "---------------------------------------"
echo $(date)
echo "usage: ./Maei_v5.sh <List_of_Genome_FASTAs> <Sequence identity threshold for cd-hit, e.g. 0.8>"
python -c "print('=' * 75)"

#Loop to identify mimps in each genome
for i in $(cat ${gL}); do
  echo "Building ouptput directories for ${i}..." ; 
  #Make Directories to save output. 
  mkdir ${i}_Results ;
  mkdir ${i}_Results/Mimps ;
  mkdir ${i}_Results/Index ;
  mkdir ${i}_Results/Augustus ;
  mkdir ${i}_Results/ORFfinding ;
  mkdir ${i}_Results/EffectorFiltering ;
  mkdir ${i}_Results/EffectorFiltering/SignalP ;
  mkdir ${i}_Results/EffectorFiltering/EffectorP ;

  #Create an index file for bedtools
  echo "Generating Index files..." ; 
  samtools faidx ${i}; 
  cut -f 1,2 ${i}.fai > ${i}_Results/Index/${i}.chrom.sizes ;
  awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ${i}.fai > ${i}_Results/Index/${i}.bed ;
  mv ${i}.fai ${i}_Results/Index/ ;


  ####################
  ###MIMP SEARCHING###
  ####################
  echo "Running mimp identification..." ; 
  #Run the HMMER mimp finder
  ${Run_nhmmer} --tblout "${i}_Results/Mimps/${i}_mimp_hmm_hits.tblout" -o "${i}_Results/Mimps/${i}_mimp_hmm_search.log"  --dna Mimp.hmm ${i} ;
  awk 'BEGIN { OFS="\t" }  {print $1,$7,$8,$12}' "${i}_Results/Mimps/${i}_mimp_hmm_hits.tblout" > ${i}_Results/Mimps/${i}_mimp_hmm_hits.bed
  #Sort the output from hmmer search
  awk '{if ($2>$3)print $1,$3,$2,".",".","-";else print $1,$2,$3,".",".","+";}' OFS='\t' ${i}_Results/Mimps/${i}_mimp_hmm_hits.bed | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t' >> ${i}_Results/Mimps/${i}_mimp_hits.ordered.bed
  #Run the mimp_finder script.
  ${Run_Bio_python} ${Mimp_finditer} ${i} > ${i}_Results/Mimps/${i}_mimp_regex_search.log;
  #Sort the output from mimp_finditer.py
  mv ${i}_mimp_regex_hits.bed ${i}_Results/Mimps/ ;
  awk '{if ($2>$3)print $1,$3,$2,".",".","-";else print $1,$2,$3,".",".","+";}' OFS='\t' ${i}_Results/Mimps/${i}_mimp_regex_hits.bed | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t' >> ${i}_Results/Mimps/${i}_mimp_hits.ordered.bed;

  ##########################################
  ## MIMP REGION PROCESSING AND EXPANSION ##
  ##########################################
  bedtools sort -i ${i}_Results/Mimps/${i}_mimp_hits.ordered.bed > ${i}_Results/Mimps/${i}_mimp_hits.sorted.bed ; #Sort and Merge the outputs from both mimp searching approaches.
  bedtools merge -c 4,5,6 -o first,first,first -i ${i}_Results/Mimps/${i}_mimp_hits.sorted.bed > ${i}_Results/Mimps/${i}_mimp_hits.merged.bed ;  #Sort and Merge the outputs from both mimp searching approaches.
  bedtools slop -i ${i}_Results/Mimps/${i}_mimp_hits.merged.bed -g ${i}_Results/Index/${i}.chrom.sizes -b 2500 > ${i}_Results/Mimps/${i}_mimp_hits.sorted.2500.bed ; #Extract the 2.5kb around the mimps as MIMP REGIONS
  bedtools merge -s -c 4,5,6 -o first,first,first -i ${i}_Results/Mimps/${i}_mimp_hits.sorted.2500.bed > ${i}_Results/Mimps/${i}_mimp_hits.sorted.2500.merged.bed ; #Merge again, as there will now be intesections of yje mimp regions, and this will limit duplications in the gff. 
  ${Run_Bio_python} ${bed2gff} ${i}_Results/Mimps/${i}_mimp_hits.sorted.2500.bed mimp_region mimp_region > ${i}_Results/Mimps/${i}_mimp_hits.gff ;    #Run bed2gff python script to generate mimp region GFF. 
  bedtools sort -faidx ${i}_Results/Index/${i}.fai -i ${i}_Results/Mimps/${i}_mimp_hits.gff > ${i}_Results/Mimps/${i}_mimp_regions.gff ; #Sort the mimp search region GFFs. 
  cp ${i}_Results/Mimps/${i}_mimp_regions.gff ${i}_Results/Index/${i}_mimp_regions.gff #Ensure that the mimp region gff can also be found in the index directory. 
  bedtools complement -i ${i}_Results/Mimps/${i}_mimp_regions.gff -g ${i}_Results/Index/${i}.bed | awk '{if ($3 <=0) {next} print $0}' > ${i}_Results/Index/${i}_non-mimp_regions.gff   #Use bedtools complement to generate non-mimp regions of the genome. If there are regions which go from 0 - <0, they are removed by awk as they will generate an error for bedtools
  bedtools maskfasta -fi ${i} -bed ${i}_Results/Index/${i}_non-mimp_regions.gff -fo ${i}_Results/Mimps/${i}.MIMP-REGIONS-ONLY.fasta ; #Now generate a fasta where the non-mimp regions have been masked, this will prevent Augustus and GetORF from annotating the non-mimps, but will also mean that the candidate effector's location in the genome fasta is preserved. 

  #Extract regions for Augustus to search 
  bedtools slop -i ${i}_Results/Mimps/${i}_mimp_hits.merged.bed -g ${i}_Results/Index/${i}.chrom.sizes -b 10000 > ${i}_Results/Augustus/${i}.Augustus_search_regions.bed ;   #Extract 20kb for Augutsus to search so that genes in mimp regions are not trucated.
  bedtools sort -i ${i}_Results/Augustus/${i}.Augustus_search_regions.bed > ${i}_Results/Augustus/${i}.Augustus_search_regions.sorted.bed ; #Sort the bedfile.
  ${Run_Bio_python} ${bed2gff} ${i}_Results/Augustus/${i}.Augustus_search_regions.sorted.bed Augustus_region Augustus_region > ${i}_Results/Augustus/${i}.Augustus_search_regions.gff ; #convert the bed file to GFF format 
  bedtools sort -faidx ${i}_Results/Index/${i}.fai -i ${i}_Results/Augustus/${i}.Augustus_search_regions.gff  > ${i}_Results/Augustus/${i}.Augustus_search_regions.sorted.gff  ; #Sort the augustus search region GFFs according to the input genome. 
  bedtools complement -i ${i}_Results/Augustus/${i}.Augustus_search_regions.sorted.gff -g ${i}_Results/Index/${i}.bed | awk '{if ($3 <=0) {next} print $0}' > ${i}_Results/Augustus/${i}.non-Augustus_search_regions.sorted.bed #Use bedtools complement to generate non-augustus search regions of the genome, if there are regions which go from 0 - 0 (Augustus search loaction goes from start of the contig), they are removed by awk as they will generate an error for bedtools. 
  bedtools maskfasta -fi ${i} -bed ${i}_Results/Augustus/${i}.non-Augustus_search_regions.sorted.bed -fo ${i}_Results/Augustus/${i}.Augustus_search-REGIONS-ONLY.fasta ; #Now generate a fasta where the non-augustus search regions have been masked, this will prevent Augustus and GetORF from annotating the non-mimps, but will also mean that the candidate effector's location in the genome fasta is preserved. 
  echo "Done."
  python -c "print('*' * 20)"
  done

python -c "print('=' * 75)"

  ##########################################
  ######### GENE AND ORF PREDICTION ########
  ##########################################
for i in $(cat ${gL});
  do
  echo "Running Augustus predictions for augustus search regions in ${i}...";
  ${Run_Augustus} --species=fusarium ${i}_Results/Augustus/${i}.Augustus_search-REGIONS-ONLY.fasta > ${i}_Results/Augustus/${i}.augustus.gff; #Run Augustus on the Augustus search regions using Fusarium as the species model.
  echo "Identifying intersections with mimp regions and generating FASTA...";
  bedtools intersect -wa -a ${i}_Results/Augustus/${i}.augustus.gff -b ${i}_Results/Mimps/${i}_mimp_hits.gff > ${i}_Results/Augustus/${i}.augustus-mimp_intersections.gff #identifying gene models that intersect with mimp regions. This prevents gene models which extend beyond the mimp regions from being truncated.
  ${Run_AGAT_manage} -gff ${i}_Results/Augustus/${i}.augustus-mimp_intersections.gff -att gene_id/Name -cp -o ${i}_Results/Augustus/${i}.augustus-mimp_intersections.agat.gff 1>${i}_Results/Augustus/${i}.augustus-mimp_intersections.agat_manage.log
  ${Run_AGAT_extract} -g ${i}_Results/Augustus/${i}.augustus-mimp_intersections.agat.gff -f ${i} -t cds -p --clean_final_stop -o ${i}_Results/Augustus/${i}.augustus-mimp_intersections.fasta 1>${i}_Results/Augustus/${i}.augustus-mimp_intersections.agat_extract.log #Creating an output FASTA from the gene models in mimp regions. 
  mv ${i}.augustus-mimp_intersections.agat* ./${i}_Results/Augustus/ 
  mv ${i}.index ${i}_Results/Index/
  echo "Done."
  echo "Running mimp-region ORF search...";
  ${Emboss_get_orf}  -minsize 90 -maxsize 2400 -find 1 ${i}_Results/Mimps/${i}.MIMP-REGIONS-ONLY.fasta -outseq ${i}_Results/ORFfinding/${i}.emboss.orf.fasta ;   #Run get orf on the mimp regions to find ORFs associated with mimps. Min and Max size set to 30aa and 800aa, respectively. 
  ${Run_Bio_python} ${getORF2bed} ${i}_Results/ORFfinding/${i}.emboss.orf.fasta ${i}_Results/ORFfinding/${i}.emboss.orf.bed #Convert the getORF output to a bed file. 
  ${Run_Bio_python} ${bed2gff} ${i}_Results/ORFfinding/${i}.emboss.orf.bed ORF mimp-regionORF > ${i}_Results/ORFfinding/${i}.emboss.orf.gff #Create a gff file from the getORF .bed file. 
  echo "Done."
  echo "Generating mimp-associated gene model and ORF gff and fasta..."
  touch ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.gff #Create an empty file which will contain the augutus and getORF gffs combined. 
  cat ${i}_Results/Augustus/${i}.augustus-mimp_intersections.agat.gff >> ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.gff #Add the augutsus GFF to the new augutus and getORF gffs combined GFF.
  cat ${i}_Results/ORFfinding/${i}.emboss.orf.gff >> ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.gff #Add the getORF GFF to the new augutus and getORF gffs combined GFF.

  ##I think I need to use Bedtools merge here again to make sure that no ORFs and gene models are in the same place and prioritise the gene models too. 
  ##I'll need to then compare the merged gff to the fasta files from getorf and Augustus to pull out the correct sequences into the ${i}_mimp-region_GeneModelsAndORFs.fasta to ensure there are no duplicates. 

  echo "${i}_mimp-region_GeneModelsAndORFs.gff."
  touch ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.fasta #Create an empty file which will contain the mimp-region augutus gene models and getORF fasta combined.
  cat ${i}_Results/Augustus/${i}.augustus-mimp_intersections.fasta >> ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.fasta #Add the mimp-region augutus gene models to the new combined fasta.
  cat ${i}_Results/ORFfinding/${i}.emboss.orf.fasta >> ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.fasta #Add the getORF mimp-region ORFs to the new combined fasta.
  echo "${i}_mimp-region_GeneModelsAndORFs.fasta."
  python -c "print('*' * 20)"
done

python -c "print('=' * 75)"


for i in $(cat ${gL});
  do echo "Generating Effector predictions for ${i}...";
  ##########################################
  ######### SEQUENCE SIZE FILTERING ########
  ##########################################
  echo "Filtering by size (<450aa && >30aa)..."
  samtools faidx ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.fasta ; #Generate index file for size filtering.
  awk '{if($2 < 450 && $2 > 30) print $1 "\t0\t" $2 "\t"}' ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.fasta.fai > ${i}_Results/EffectorFiltering/SignalP/${i}_.SizeFiltered.bed     #Filter all sequences by size. <300 aa accepted.
  bedtools sort -i ${i}_Results/EffectorFiltering/SignalP/${i}_.SizeFiltered.bed > ${i}_Results/EffectorFiltering/SignalP/${i}_.SizeFiltered.sorted.bed ; #Sort the bed files for getfasta
  bedtools getfasta -s -fi ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.fasta -bed ${i}_Results/EffectorFiltering/SignalP/${i}_.SizeFiltered.sorted.bed  -fo ${i}_Results/EffectorFiltering/SignalP/${i}_.SizeFiltered.sorted.fasta ;     #Create a fasta containing the size filtered sequences.
  echo "Done." ;
  ##########################################
  ############ SIGNALP FILTERING ###########
  ##########################################
  echo "Running SignalP...";
  ${Run_SignalP} -f summary ${i}_Results/EffectorFiltering/SignalP/${i}_.SizeFiltered.sorted.fasta > ${i}_Results/EffectorFiltering/SignalP/${i}_SignalP.log;
  echo "Filtering SignalP results...";
  grep "^Name*" ${i}_Results/EffectorFiltering/SignalP/${i}_SignalP.log | awk '$2 ~ /YES/ {print $0}' > ${i}_Results/EffectorFiltering/SignalP/${i}_SignalPepPositives.out ;
  echo "Generating fasta...";
  awk  '{gsub("Name=","\t",$0); print;}' ${i}_Results/EffectorFiltering/SignalP/${i}_SignalPepPositives.out | awk '{print $1}' > ${i}_Results/EffectorFiltering/SignalP/${i}_SignalPepPositives_list.txt;
    for j in $(cat ${i}_Results/EffectorFiltering/SignalP/${i}_SignalPepPositives_list.txt;); do samtools faidx ${i}_Results/EffectorFiltering/SignalP/${i}_.SizeFiltered.sorted.fasta ${j}; done > ${i}_Results/EffectorFiltering/SignalP/${i}_SignalPepPositives.fasta;
  echo "Done." ;
  
  ##########################################
  ######### EFFECTORP PREDICTIONS ##########
  ##########################################
  echo "Running EffectorP predictions for ${i}"; 
  python ${Run_EffectorP} -i ${i}_Results/EffectorFiltering/SignalP/${i}_SignalPepPositives.fasta -E ${i}_Results/EffectorFiltering/EffectorP/${i}_EffectorP.fasta  > ${i}_Results/EffectorFiltering/EffectorP/${i}_EffectorP.log ;
  echo "Done."
  echo "Creating final candidate effector gff..."; 
  grep ">" ${i}_Results/EffectorFiltering/EffectorP/${i}_EffectorP.fasta  | sed 's/>//g' | sed 's/\:0.*//' | sed 's/\.*t1//' > ${i}_Results/Index/${i}_EffectorList.txt ; #Split the fasta header into just sequence names that can be found in the gene and ORF gff and save as a file.
  grep  -f ${i}_Results/Index/${i}_EffectorList.txt  ${i}_Results/Index/${i}_mimp-region_GeneModelsAndORFs.gff > ${i}_Results/${i}_CandidateEffector.gff #Using the file which lists the fasta headers from EffectorP, grep the line containing those headers from the genes and ORF gff and save as a new gff.
  sed 's/\:0.*//' ${i}_Results/EffectorFiltering/EffectorP/${i}_EffectorP.fasta | sed 's/\.*t1//' > ${i}_Results/${i}_CandidateEffector.fasta
  echo "Done."
  python -c "print('*' * 20)"
done

python -c "print('=' * 75)"

  ##########################################
  ############ CD-HIT CLUSTERING ###########
  ##########################################
touch ./AllCandidateEffectorSets.fasta #Create an empty file to hold all of the candidate effectors.
echo "Clustering final effector sets..."
for i in $(cat ${gL}); 
  do
  sed -i "s/^>/>${i};;_/" "$i" ${i}_Results/${i}_CandidateEffector.fasta ; #Add genome filename to start of Fasta Headers so we know which isolate this came from. 
  cat ${i}_Results/${i}_CandidateEffector.fasta >> ./AllCandidateEffectorSets.fasta ; #Combine all of the individual candidate effector sets by adding them to the empty candidate effector FASTA.
done 

mkdir cdhit
${Run_cdhit} -i ./AllCandidateEffectorSets.fasta -d 0 -o ./AllCandidateEffectorSets -c ${Cl} -n 5  -G 1 -g 0 -b 20 -l 10 -s 0.0 -aL 0.0 -aS 0.0 1> ./cdhit/cd-hit.log #Cluster the combined effector set based on sequnce identitiy.
mv ./AllCandidateEffectorSets ./cdhit/./AllCandidateEffectorSet_CulsterRepresentativeSeqs.fasta
echo "Processing Cd-Hit output..."
${Run_Bio_python} ${ProcessingCDHIT} ./AllCandidateEffectorSets.clstr 
echo "Done."

python -c "print('=' * 75)"
echo "Mimp-Associated Effector Identification Complete"
echo "------------------------------------------------"
echo $(date)
python -c "print('=' * 75)"
