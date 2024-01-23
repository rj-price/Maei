# Development of Mimp-associated Effector Identification Script (Maei.sh)

An mimp associated effector script was originally generated using BLASTN based on the method from Chang et al (2020). 

The original script/pipeline can be found:

My Machine:

/Users/u1983390/Scripts/MIMP_associated_effector_prediction.sh

Vettel:

/home/u1983390/Scripts/MIMP_associated_effector_prediction.sh

Shared:

/Volumes/shared-1/HCSS3/Shared258/Jamie/Scripts/MIMP_associated_effector_prediction.sh

The BLASTN and subsequent filtering step was removed and replaced with my [mimp_finditer.py](http://mimp_finditer.py/) script. Similarly, the BLASTCLUST step was removed and replaced with two CD-HIT steps to generate a non-redundant set of candidate effectors for each genome. then a clustered group of effectors to determine if they can be grouped by race. 

The new, adapted script for effector prediction script ([Maei.sh](http://maei.sh/)) can be found: 

My Machine:

/Users/u1983390/Scripts/Maei.sh 

Vettel:

/home/u1983390/Scripts/Maei.sh

Shared:

/Volumes/shared-1/HCSS3/Shared258/Jamie/Scripts/Maei.sh 

To run the script, [Maei.sh](http://maei.sh/), create a directory containing all of the genomes and a text file listing the genomes to run the analysis on. Run [Maei.sh](http://maei.sh/) from this directory. Ensure that are no spaces in FASTA headers. 

---

**Workflow for Maei.sh�**�

[MY_effector_pred_pipeline.pdf](./file/MY_effector_pred_pipeline.pdf)

---

**Maei.sh Run�**�

[Maei.sh](http://maei.sh/) was run on the 13th of August 2020 on genomes listed in Table 1. 

|Table 1: Genomes used in the identification of _mimps_ and _SIX_ genes[[1]](applewebdata://BC134C0E-EE64-4532-B542-43706AD757D6#_ftn1).|                                                                   |              |          |        |                   |                                  |                                                                                                            |
|---------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|--------------|----------|--------|-------------------|----------------------------------|------------------------------------------------------------------------------------------------------------|
|**Assembly Accession**                                                                                                                 |**Bioproject Accession**                                           |**Isolate**   |**Source**|**Race**|**Contig N50 (MB)**|**Predicted Protein Coding genes**|**BUSCO**** – intact single copy orthologs**[[2]](applewebdata://BC134C0E-EE64-4532-B542-43706AD757D6#_ftn2)|
|[GCA_000260195.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000260195.1)                                                               |[PRJNA73539](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA73539)   |54006 (II5)   |NCBI      |TR4     |0.350              |16.638                            |99.0%                                                                                                       |
|[GCA_000350345.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_000350345.1)                                                               |[PRJNA174274](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA174274) |Foc1_1.1 or N2|NCBI      |R1      |0.023              |18,065                            |89.6%                                                                                                       |
|GCA_000350365.1                                                                                                                        |[PRJNA174275](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA174275) |Foc4_1.0 or B2|NCBI      |R4      |0.016              |17,462                            |80.5%                                                                                                       |
|GWHAAST00000000                                                                                                                        |[PRJCA001282](https://bigd.big.ac.cn/bioproject/browse/PRJCA001282)|Foc1 60       |NGDC      |R1      |2.134              |15,865                            |97.7%                                                                                                       |
|GWHAASU00000000                                                                                                                        |[PRJCA001282](https://bigd.big.ac.cn/bioproject/browse/PRJCA001282)|FocTR4 58     |NGDC      |TR4     |4.378              |15,519                            |99.6%                                                                                                       |
|GCA_001696625.1                                                                                                                        |[PRJNA280907](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA280907/)|C1HIR_9889    |NCBI      |R4      |0.087              |No annotation available           |99.5%                                                                                                       |
|[GCA_005930515.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_005930515.1)                                                               |PRJNA529756                                                        |160527        |NCBI      |R1      |4.885              |16,536                            |99.1%                                                                                                       |
|                                                                                                                                       |                                                                   |              |          |        |                   |                                  |                                                                                                            |
|GCA_011316005.1                                                                                                                        |[PRJNA552447](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA552447/)|TC1-1         |NCBI      |Unknown |0.067              |No annotation available           |98.6%                                                                                                       |
|[GCA_007994515.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_007994515.1)                                                               |PRJNA556111                                                        |UK0001        |NCBI      |TR4     |4.494              |14,472                            |98.4%                                                                                                       |
|                                                                                                                                       |                                                                   |              |          |        |                   |                                  |                                                                                                            |

[[1]](applewebdata://BC134C0E-EE64-4532-B542-43706AD757D6#_ftnref1) NCBI refers to genomes from National Center for Biotechnology Information, sourced via GenBank website ([https://www.ncbi.nlm.nih.gov/genbank/](https://www.ncbi.nlm.nih.gov/genbank/)).   NDGC refers to genomes from National Genomics Data Centre website ([https://bigd.big.ac.cn/](https://bigd.big.ac.cn/)).

[[2]](applewebdata://BC134C0E-EE64-4532-B542-43706AD757D6#_ftnref2)BUSCO - based on genome - % Intact and single copy orthologs (_Sordariomycetes_ data). 

---

**Results**

The [mimp_finditer.py](http://mimp_finditer.py/) script found a similar number of mimps per genome to [FoEC.py](http://foec.py/) (Table 2). 

|Table 2: Number of mimps identified with mime_finditer.py as part of Maei.sh compared to the number of mimps found by FoEC.py |                                            |                                   |
|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------|-----------------------------------|
|**Genome/Isolate**                                                                                                            |**mimp_finditer.py number of mimps found�**�|**FoEC.py number of mimps found�**�|
|FocTR4 58                                                                                                                     |24                                          |23                                 |
|C1HIR_9889                                                                                                                    |19                                          |17                                 |
|Foc1 60                                                                                                                       |5                                           |5                                  |
|Foc1_1.1 or N2                                                                                                                |22                                          |20                                 |
|160527                                                                                                                        |33                                          |31                                 |
|UK0001                                                                                                                        |23                                          |22                                 |
|TC1-1                                                                                                                         |22                                          |21                                 |
|Foc4_1.0 or B2                                                                                                                |17                                          |17                                 |
|54006 (II5)                                                                                                                   |12                                          |11                                 |

Excluding the EffectorP(2.0.1) step, [Maei.sh](http://maei.sh/) found 62 candidate effectors. 

Results can be found: 

Vettel:

/home/u1983390/Fusarium_data/MIMPS/MIMP_Maei_Searches_13.08.2020

My Machine:

/Users/u1983390/Fusarium_data/MIMPS/MIMP_Maei_Searches_13.08.2020

Shared:

/Volumes/shared-1/HCSS3/Shared258/Jamie/MIMPS/MIMP_Maei_Searches_13.08.2020

The predicted effectors clustered into 62 groups (Table 3), some of which were homologous to the predicted effectors from FoEC run (04/08/2020). 

Table 3: MIMP ASSOCIATED EFFECTORS FROM Maei CD-HIT CLUSTERS FOR ALL GENOMES

[MIMP_ASSOCIATED_EFFECTORS_FROM_Maei_CD-HIT_CLUSTERS_FOR_ALL_GENOMES.xlsx](./file/MIMP_ASSOCIATED_EFFECTORS_FROM_Maei_CD-HIT_CLUSTERS_FOR_ALL_GENOMES.xlsx)

When the EffectorP step is included 14 candidate effectors are identified (Table 4). 

|Table 4: Effector clusters identified using the pipeline developed in house, their SIX homologs determined using BLASTP (1e-6), and the predicted effectors from the van Dam _et al._ (2016) method with which they cluster.|                      |            |        |                                                                  |               |
|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------|------------|--------|------------------------------------------------------------------|---------------|
|**Effector Cluster**                                                                                                                                                                                                        |**Genomes in cluster**|**Isolate** |**Race**|**Clusters with van Dam** **et al.** **(2016) predicted effector**|**SIX homolog**|
|Cluster_6                                                                                                                                                                                                                   |GCA_000350365.1       |Foc4_1 or B2|R4      |X0004..MAPYSMVLLGALSILGFGAYA                                      |SIX1           |
|GCA_005930515.1                                                                                                                                                                                                             |160527                |R1          |        |                                                                  |               |
|GCA_007994515.1                                                                                                                                                                                                             |UK0001                |TR4         |        |                                                                  |               |
|GWHAASU00000000                                                                                                                                                                                                             |FocTR4_58             |TR4         |        |                                                                  |               |
|Cluster_8                                                                                                                                                                                                                   |GCA_000350345.1       |Foc1_1 or N2|R1      |X0015..MAKSLKLLMLIACMFIACTRA                                      |               |
|GCA_005930515.1                                                                                                                                                                                                             |160527                |R1          |        |                                                                  |               |
|Cluster_9                                                                                                                                                                                                                   |GCA_000260195.2       |FO_II5_V1   |TR4     |                                                                  |               |
|GCA_000350345.1                                                                                                                                                                                                             |Foc1_1 or N2          |R1          |        |                                                                  |               |
|GCA_000350365.1                                                                                                                                                                                                             |Foc4_1 or B2          |R4          |        |                                                                  |               |
|GCA_005930515.1                                                                                                                                                                                                             |160527                |R1          |        |                                                                  |               |
|GCA_007994515.1                                                                                                                                                                                                             |UK0001                |TR4         |        |                                                                  |               |
|GWHAASU00000000                                                                                                                                                                                                             |FocTR4_58             |TR4         |        |                                                                  |               |
|Cluster 11                                                                                                                                                                                                                  |GCA_000260195.2       |FO_II5_V1   |TR4     |X0006..MHFTTAALSALLASAVSA                                         |               |
|GCA_000350365.2                                                                                                                                                                                                             |Foc4_1 or B2          |R4          |        |                                                                  |               |
|GCA_001696625.1                                                                                                                                                                                                             |C1HIR_9889            |R4          |        |                                                                  |               |
|GCA_007994515.1                                                                                                                                                                                                             |UK0001                |TR4         |        |                                                                  |               |
|Cluster 13                                                                                                                                                                                                                  |GCA_000350365.1       |Foc4_1 or B2|R4      |X0002..MTRFHLILLPLLFSWFSYCFG                                      |SIX13          |
|GCA_001696625.1                                                                                                                                                                                                             |C1HIR_9889            |R4          |        |                                                                  |               |
|GCA_007994515.1                                                                                                                                                                                                             |UK0001                |TR4         |        |                                                                  |               |
|GWHAASU00000000                                                                                                                                                                                                             |FocTR4_58             |TR4         |        |                                                                  |               |
|Cluster 15                                                                                                                                                                                                                  |GCA_005930515.1       |160527      |R1      |                                                                  |               |
|GCA_011316005.1                                                                                                                                                                                                             |TC1-1                 |Unknown     |        |                                                                  |               |
|Cluster 17                                                                                                                                                                                                                  |GCA_000260195.2       |FO_II5_V1   |TR4     |                                                                  |SIX9           |
|GCA_000350345.1                                                                                                                                                                                                             |Foc1_1 or N2          |R1          |        |                                                                  |               |
|GCA_001696625.1                                                                                                                                                                                                             |C1HIR_9889            |R4          |        |                                                                  |               |
|GCA_005930515.1                                                                                                                                                                                                             |160527                |R1          |        |                                                                  |               |
|GCA_007994515.1                                                                                                                                                                                                             |UK0001                |TR4         |        |                                                                  |               |
|GCA_011316005.1                                                                                                                                                                                                             |TC1-1                 |Unknown     |        |                                                                  |               |
|GWHAASU00000000                                                                                                                                                                                                             |FocTR4_58             |TR4         |        |                                                                  |               |
|Cluster 18                                                                                                                                                                                                                  |GCA_000350345.1       |Foc1_1 or N2|R1      |                                                                  |               |
|GCA_005930515.1                                                                                                                                                                                                             |160527                |R1          |        |                                                                  |               |
|Cluster 22                                                                                                                                                                                                                  |GWHAAST00000000       |Foc1_60     |R1      |                                                                  |               |
|Cluster 27                                                                                                                                                                                                                  |GCA_005930515.1       |160527      |R1      |                                                                  |               |
|Cluster 30                                                                                                                                                                                                                  |GCA_000350345.1       |Foc1_1 or N2|R1      |X0013..MVLSKLLSSAVMAAVQA                                          |               |
|GCA_005930515.1                                                                                                                                                                                                             |160527                |R1          |        |                                                                  |               |
|GCA_011316005.1                                                                                                                                                                                                             |TC1-1                 |Unknown     |        |                                                                  |               |
|Cluster 31                                                                                                                                                                                                                  |GCA_000260195.2       |FO_II5_V1   |TR4     |X0017..MRFSNIAVCLLTAISGADA                                        |SIX9           |
|GCA_000350365.1                                                                                                                                                                                                             |Foc4_1 or B2          |R4          |        |                                                                  |               |
|GCA_001696625.1                                                                                                                                                                                                             |C1HIR_9889            |R4          |        |                                                                  |               |
|GCA_007994515.1                                                                                                                                                                                                             |UK0001                |TR4         |        |                                                                  |               |
|GWHAASU00000000                                                                                                                                                                                                             |FocTR4_58             |TR4         |        |                                                                  |               |
|Cluster 53                                                                                                                                                                                                                  |GWHAAST00000000       |Foc1_60     |R1      |                                                                  |               |
|Cluster 56                                                                                                                                                                                                                  |GCA_005930515.1       |160527      |R1      |                                                                  |               |

---

**BREAKDOWN OF Maei.sh Script�**�

Jamie Pike - created 12/8/2020

Using a dircetory containing genomes and a .txt file containing a list of those genomes, find mimp sequences present and predicted mimp-associated effectors.

    

```
#PLEASE SET PATHS TO THE FOLLOWING DEPENDENCIES:
    Run_Bio_python=/home/u1983390/miniconda3/envs/biopythonEnv/bin/python3           #path to biopython
    Mimp_finditer=/home/u1983390/Scripts/Mimp_finditer.py                            #path to mimp_finditer.py script
    Run_Augustus=/home/u1983390/miniconda3/envs/augustusEnv/bin/augustus             #path to Augustus executable
    Run_getAnnoFasta=/home/u1983390/miniconda3/envs/augustusEnv/bin/getAnnoFasta.pl  #path to Augustus getAnnoFasta.pl
    Emboss_get_orf=/home/u1983390/miniconda3/envs/EmbossEnv/bin/getorf               #path to EMBOSS getorf executable
    Run_SignalP=/home/u1983390/apps/FoEC/signalp-4.1/signalp                         #path to SignalP executable
    Run_EffectorP=/home/u1983390/apps/EffectorP-2.0-2.0.1/Scripts/EffectorP.py       #path to EffectorP.py
```

Insert via commandline the name of file containing a list of genomes

```
    gL=${1?Error: no list of FASTAs given. \n Please provide a .txt file containing a list of FASTAs for analysis}

```

Run the mimp_finditer.py script

Loop through a list of genomes running the mimp_finditer.py script

```
    for i in $(cat ${gL}); do

      #Run the mimp_finder script.
      $Run_Bio_python $Mimp_finditer ${i} ;
      #Create an index file for bedtools
      samtools faidx ${i}; cut -f 1,2 ${i}.fai > ${i}.chrom.sizes ;
      #Sort the output from mimp_finditer.py
      bedtools sort -i ${i}_mimp_hits.bed > ${i}_mimp_hits.sorted.bed ;
      #Create a fasta containing the mimp sequence
      bedtools slop -i ${i}_mimp_hits.sorted.bed -g ${i}.chrom.sizes -b 0 > ${i}_mimp_hits.sorted.0.bed;
      bedtools getfasta -fi ${i} -bed ${i}_mimp_hits.sorted.0.bed -fo ${i}_mimp_hits.0.fasta;
      #create a fasta by expanding 2.5kb upstream and downstream of the mimp bases.
      bedtools slop -i ${i}_mimp_hits.sorted.bed -g ${i}.chrom.sizes -b 2500 > ${i}_mimp_hits.sorted.2500.bed;
      bedtools getfasta -fi ${i} -bed ${i}_mimp_hits.sorted.2500.bed -fo ${i}_mimp_hits.2500.fasta; done
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
    python -c "print('\n')"
    echo "Creating directories"
```

Taking all of the previous outputs, create a folder for each genome in gL placing the all of the ouputs for the genome into the folder.

```
      for i in $(cat ${gL});
        do mkdir ${i}_MIMPS; mv ./*${i}* ./${i}_MIMPS/; done
      echo "//Ignore the mv: rename error message"
```

This does produce an error message saying that the folder cannot be moved into itself. This can be ignored as the folder should not be moved into itself

```
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
```

Run AUGUSTUS on sequences extracted 2.5k upstream and downstream of the mimp

```
    for i in $(cat ${gL});
      do
      echo "//Running predictions for ${i}";
      ${Run_Augustus} --species=fusarium ${i}_MIMPS/${i}_mimp_hits.2500.fasta > ${i}_MIMPS/${i}_mimp_hits.2500.augustus.gff;
      echo "//Creating fasta for ${i}";
      ${Run_getAnnoFasta} ${i}_MIMPS/${i}_mimp_hits.2500.augustus.gff; done
    for i in $(cat ${gL});
      do
      sed "s/>/>${i}_Candidate_Sequence_/" ${i}_MIMPS/${i}_mimp_hits.2500.augustus.aa
      > ${i}_MIMPS/${i}_mimp_hits.2500.augustus.fsa_aa; done
    python -c "print('=' * 80)"
    echo "Running getorf on mimps"
    python -c "print('=' * 80)"
```

Run EMBOSS getorf on sequences extracted 2.5k upstream and downstream of the mimp

```
    for i in $(cat ${gL});
      do echo "//Running predictions for ${i}";
      ${Emboss_get_orf} -find 1 ${i}_MIMPS/${i}_mimp_hits.2500.fasta -outseq ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf ; done
    for i in $(cat ${gL});
      do
      sed "s/>/>${i}_Candidate_Sequence_/" ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf
      > ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.fsa_aa ; done
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
    echo "Running SignalP on predictions..."
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
    echo "AUGUSTUS predictions"
    python -c "print('=' * 80)"
```

AUGUSTUS SIGNALP ANALYSIS AND FILTERING

Loop through the genomes running SignalP on the predicted genes from Augustus

```
    for i in $(cat ${gL});
      do echo "//Running SignalP for ${i} Augustus predictions";
      ${Run_SignalP} -f summary  ./${i}_MIMPS/${i}_mimp_hits.2500.augustus.fsa_aa
      > ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.augustus.SignalP_raw.out;
      echo "//Filtering SignalP results for ${i} Augustus predictions";
      grep "^Name*" ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.augustus.SignalP_raw.out | awk '$2 ~ /YES/ {print $0}'
      > ./${i}_MIMPS/${i}_mimp_hits.2500.augustus.SignalP.out ; done
```

Filter the SignalP results generating a FASTA which contains the predicted genes with a signal peptide

```
    for i in $(cat ${gL});
      do echo "//Generating fasta of Augustus predictions for ${i} which contain a signal peptide";
      awk  '{gsub("=","\t",$0); print;}' ${i}_MIMPS/${i}_mimp_hits.2500.augustus.SignalP.out | awk '{print $2}'
      > ${i}_MIMPS/${i}_SignalPep_sequences_augustus_list.txt;
        for j in $(cat ${i}_MIMPS/${i}_SignalPep_sequences_augustus_list.txt;);
        do samtools faidx ${i}_MIMPS/${i}_mimp_hits.2500.augustus.fsa_aa ${j}; done
        > ${i}_MIMPS/${i}_mimp_hits.2500.augustus.SignalP.fsa_aa; done
    python -c "print('=' * 80)"
    echo "Getorf predictions"
    python -c "print('=' * 80)"
```

EMBOSS getorf SIGNALP ANALYSIS AND FILTERING

loop through genomes runnning SignalP on the getorf outputs

```
    for i in $(cat ${gL});
      do echo "//Running SignalP for ${i} getorf output";
      ${Run_SignalP} -f summary  ./${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.fsa_aa
      > ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.emboss.orf.SignalP_raw.out;
      echo "//Filtering SignalP results for ${i} getorf output";
      grep "^Name*" ./${i}_MIMPS/${i}_SignalP_Raw/${i}_mimp_hits.2500.emboss.orf.SignalP_raw.out
      | awk '$2 ~ /YES/ {print $0}' > ./${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.SignalP.out ; done
```

Filter the SignalP results generating a FASTA which contains the predicted genes with a signal peptide

```
    for i in $(cat ${gL});
      do echo "//Generating fasta of getorf output for ${i} which contain a signal peptide";
      awk  '{gsub("=","\t",$0); print;}' ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.SignalP.out
      | awk '{print $2}' > ${i}_MIMPS/${i}_SignalPep_sequences_emboss.orf_list.txt;
        for j in $(cat ${i}_MIMPS/${i}_SignalPep_sequences_emboss.orf_list.txt;);
        do samtools faidx ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.fsa_aa ${j}; done
        > ${i}_MIMPS/${i}_mimp_hits.2500.emboss.orf.SignalP.fsa_aa ; done
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
    echo "Tidying up"
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
```

Organise the files into directories

```
    for i in $(cat ${gL});
      do
      mv ./${i}_MIMPS/${i}_SignalPep_* ./${i}_MIMPS/${i}_SignalP_Raw/ ;
      mv ./${i}_MIMPS/${i}*.SignalP.out ./${i}_MIMPS/${i}_SignalP_Raw/ ;
      mv ./${i}_MIMPS/${i}*.2500.augustus.aa ./${i}_MIMPS/${i}_Predictions/ ;
      mv ./${i}_MIMPS/${i}*.2500.emboss.orf ./${i}_MIMPS/${i}_Predictions/ ;
      mv ./${i}_MIMPS/${i}*.augustus.fsa_aa.fai  ./${i}_MIMPS/${i}_Predictions/ ;
      mv ./${i}_MIMPS/${i}*.emboss.orf.fsa_aa.fai ./${i}_MIMPS/${i}_Predictions/ ;
      mv ./${i}_MIMPS/${i}*.2500.augustus.gff ./${i}_MIMPS/${i}_Predictions/ ; done
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
    echo "Creating non-redundant set of predicted protiens with signal peptide..."
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
```

Cluster genomes using CD-HIT

Create fasta of all the augustus and getorf results and place in one fasta file

```
    for i in $(cat ${gL});
     do
       echo "//Creating fasta from the augutsus and getorf predictions for ${i}" ;
       cat ${i}_MIMPS/${i}*.SignalP.fsa_aa > ${i}_MIMPS/${i}.SignalP.fsa_aa ; done
```

Create a non-redundant protein set from the grouped fasta created above.

```
     for i in $(cat ${gL});
      do
        echo "//Creating non-redundant fasta from the augutsus and getorf predictions for ${i}" ;
        cd-hit -i ${i}_MIMPS/${i}.SignalP.fsa_aa -d 0 -o ${i}_MIMPS/${i}.SignalP.non_redundant.fasta
        -c 1 -n 5  -G 1 -g 1 -b 20 -l 10 -s 0.0 -aL 1 -aS 0.0 ; done
```

Create a master fasta file containing all of the candidate effectors from all of the genomes.

```
    echo "//Creating master fasta file for candiate sequences "
     for i in $(cat ${gL});
      do
      cat ${i}_MIMPS/${i}.SignalP.non_redundant.fasta ; done  > ./all_mimp_associated_effectors_unclustered.fasta
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
    echo "Clustering the non-redundant protien set..."
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"

    cd-hit -i all_mimp_associated_effectors_unclustered.fasta -d 0 -o all_mimp_associated_effectors_CLUSTERED.fasta
    -c 0.9 -n 5  -G 1 -g 1 -b 20 -l 10 -s 0.0 -aL 0.0 -aS 0.0

    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
    echo 'Running EffectorP on all the  clustered candidate effectors...'
    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
```

Run EffectorP on outputs

  run EffectorP for effector prediction on all filtered hits with a signal peptide and produce a fasta with the output.

```
    python ${Run_EffectorP} -i all_mimp_associated_effectors_CLUSTERED.fasta
    -E all_mimp_associated_effectors_CLUSTERED.EffectorP.fasta

    python -c "print('=' * 130)"
    python -c "print('=' * 130)"
```

---

Script:

[Maei.sh](./file/Maei.sh)

Chang _et al.,_ (2020) paper: 

[Chang_et_al_Identification_of_mimp-associated_effector_genes_in_Fusarium_oxysporum_f-sp-cubense_race_1_and_race_4_and_virulence_confirmation_of_a_candidate_effector_gene.pdf](./file/Chang_et_al_Identification_of_mimp-associated_effector_genes_in_Fusarium_oxysporum_f-sp-cubense_race_1_and_race_4_and_virulence_confirmation_of_a_candidate_effector_gene.pdf)
