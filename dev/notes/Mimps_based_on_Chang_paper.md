# Mimps based on Chang paper

Paper bookmarked on chrome under /Papers/Pathogenicity genes/MITES and mimps/ 

![1eca1295-b142-4591-a7eb-0ff467ca06c1.png](image/1eca1295-b142-4591-a7eb-0ff467ca06c1.png)

Received fasta sequences from Dr Chang, following email (subject: 13 reference mimp sequences help regarding your recent publication). 

3/3/2020

Merged the independent sequences given into a fasta containing all 14 mimps:

Focub_all_Chang_mimps.faa

Also created a folder contains all mimps with number 12 removed, as 12 and 8 are identical. 

Focub_mimp12rm_Chang_mimps.faa

**Blastn** 

Following fasta generation, BLASTn was performed using the 13 reference mimp sequences Chang sourced from Van Dam. These sequences were blasted against the TR4_Warmington genome, and the R1_Aisa genome. 

**Warmingtion as subject�**�

All mimps as query 

blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_all_Chang_mimps.faa -db ../../TR4_Warmington_Raw_Data/VMNF01.1.fsa_nt -perc_identity 90.00 -out cubense_all_mimps_public_Vs_TR4_Warmingtion_genome.blastn.out

Mimp 12 rm from query 

**Asia as subject�**�

All as query

blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_all_Chang_mimps.faa -db ../../R1_Asia_Raw_Data/SRMI01.1.fsa_nt -perc_identity 90.00 -out cubense_all_mimps_public_Vs_R1_Asia_genome.blastn.out

Mimp 12 rm from query 

blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_mimp12rm_Chang_mimps.faa -db ../../R1_Asia_Raw_Data/SRMI01.1.fsa_nt -perc_identity 90.00 -out cubense_mimp_12_rm_public_Vs_R1_Asia_genome.blastn.out

**Yun_58_TR4 as a subject�**�

All as query 

blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_all_Chang_mimps.faa -db ../../TR4_Yun_58_Raw_Data/GWHAASU00000000.genome.fasta -perc_identity 90.00 -out cubense_all_mimps_public_Vs_TR4_YUN_58_genome.blastn.out

Mimp 12 rm from query 

blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_mimp12rm_Chang_mimps.faa -db ../../TR4_Yun_58_Raw_Data/GWHAASU00000000.genome.fasta -perc_identity 90.00 -out cubense_mimp_12_rm_public_Vs_TR4_YUN_58_genome.blastn.out

**Yun_60_R1 as a subject�**�

All as query 

-query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_all_Chang_mimps.faa -db ../../R1_Yun_R1_60_Raw_Data/Yun_R1_60.faa -perc_identity 90.00 -out cubense_all_mimps_public_Vs_R1_YUN_60_genome.blastn.out

Mimp 12 rm from query 

blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_mimp12rm_Chang_mimps.faa -db ../../R1_Yun_R1_60_Raw_Data/Yun_R1_60.faa -perc_identity 90.00 -out cubense_mimp_12_rm_public_Vs_R1_YUN_60_genome.blastn.out

**FO_II5_V1 as a subject�**�

All as query 

 blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_all_Chang_mimps.faa -db ../../TR4_FO_II5_V1_Raw_Data/genome_assemblies_NUCL/ncbi-genomes-2020-02-25/GCF_000260195.1_FO_II5_V1_genomic.fna -perc_identity 90.00 -out cubense_all_mimps_public_Vs_TR4_FO_II5_V1_genome.blastn.out

Mimp 12 rm from query 

blastn -query ~/Fusarium_data/Chang_mimps/cubense_mimps_public/Focub_mimp12rm_Chang_mimps.faa -db ../../TR4_FO_II5_V1_Raw_Data/genome_assemblies_NUCL/ncbi-genomes-2020-02-25/GCF_000260195.1_FO_II5_V1_genomic.fna -perc_identity 90.00 -out cubense_mimp_12_rm_public_Vs_TR4_FO_II5_V1_genome.blastn.out

Command for calculating the total matches 

E.g. 

MACP80255:Blast_outputs u1983390$ grep -c ">" cubense_all_mimps_public_Vs_R1_Asia_genome.blastn.out 

27

Output from grep commands: 

Chang et al paper - 21 TR4 mimps, and 20 R1 mimps 

|

              |**Race**|**All mimps**|**With mimp 12 rm**|
|--------------|--------|-------------|-------------------|
|**Warmington**|TR4     |32           |29                 |
|**Asia**      |R1      |27           |23                 |
|**Yun_TR4_58**|TR4     |32           |29                 |
|**FO_II5_V1** |4       |53           |46                 |
|**Yun_R1_60** |R1      |23           |17                 |

BLASTp was performed using the **Candidate effector protein sequences** **provided in supplementary table 3 of the Chang paper**. Outputs were saved in the Blast_output folders corresponding to each genome. Not cut off was used. 

Commands: ![14bb42f4-43f8-408c-83e9-28766fd5c7a2.png](image/14bb42f4-43f8-408c-83e9-28766fd5c7a2.png)

NOTE: realised I had saved Run TR4 as 60, when its supposed to be 58, have since changed folder name but not file names. 

The outfmt 6 outputs have been saved in an excel file: ~/Documents/Lab_Book/Blast/Chang_Candidate_Effectors_vs_genomes.blastp.xlsx 
