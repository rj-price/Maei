# NCBI BLASTP for Predicted Effectors

Predicted effectors from Maei.sh and FoEC.py were searched using the NCBI database via BLASTP. 

NCBI Search link: [https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=ENH74795.1&LINK_LOC=protein&PAGE_TYPE=BlastSearch](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=ENH74795.1&LINK_LOC=protein&PAGE_TYPE=BlastSearch) 

Database searched:

Title:All non-redundant GenBank CDS translations+PDB+SwissProt+PIR+PRF excluding environmental samples from WGS projects

Molecule Type:Protein

All query files for NBCI BLASTP are stored in:

/Users/u1983390/Fusarium_data/MIMPS/MIMP_Maei_Searches_13.08.2020/Compare_Maei_to_FoEC

or

/home/u1983390/Fusarium_data/MIMPS/FoEC_and_Maei_downstream/Maei_and_FoEC_FASTAs

    drwxr-xr-x  10 u1983390  734225131   320B  7 Sep 13:29 Original_Clusters

    drwxr-xr-x   8 u1983390  734225131   256B  7 Sep 13:07 EffectorP_Clusters

    drwxr-xr-x   9 u1983390  734225131   288B  7 Sep 13:06 non_effectorP_results

---

**Output�**�

The 90 sequences were divided/searched in the following groups to generate an .xml, .txt, and .hit_table.txt  for each: 

    all_Maei_and_FoEC_comb.NON-REDUNDANT.CLUSTERED_90.Alignment.txt

    All of the sequences produced from Maei and FoEC in a non-redundant set.

    In total there are 90 sequences.

    all_Maei_and_FoEC_comb.REDUNDANT.CLUSTERED_75

    All of the sequences produced from Maei and FoEC clustered using CD-HIT (-c 0.9 -n 5  -G 1 -g 1 -b 20 -l 10 -s 0.0 -aL 0.0 -aS 0.0). 

    There are 75 sequences in total.

    Maei_and_FoEC_EffectorPs_comb.CLUSTERED.fasta_17

    All of the sequences from Maei and FoEC predicted to be effectors by EffectorP 

    17 sequences in total. 

    all_Maei_and_FoEC_comb.NON-EFFECTORS.CLUSTERED_59

    All of the sequences from Maei and FoEC not predicted to be effectors by EffectorP

    59 sequences in total 

The files below can be found in the directories:

Vettel:

/home/u1983390/Fusarium_data/MIMPS/FoEC_and_Maei_downstream/NCBI_BLASTP

My Machine:

/Users/u1983390/Fusarium_data/MIMPS/MIMP_Maei_Searches_13.08.2020/Compare_Maei_to_FoEC/NCBI_BLASTP

**Output Files:**

-rw-r--r-- 1 u1983390 u1983390 299K Sep  7 15:22 all_Maei_and_FoEC_comb.NON-REDUNDANT.CLUSTERED_90.Alignment.hit_table.txt

-rw-r--r-- 1 u1983390 u1983390 3.2M Sep  7 15:22 all_Maei_and_FoEC_comb.NON-REDUNDANT.CLUSTERED_90.Alignment.txt

-rw-r--r-- 1 u1983390 u1983390 4.0M Sep  7 15:22 all_Maei_and_FoEC_comb.NON-REDUNDANT.CLUSTERED_90.Alignment.xml

-rw-r--r-- 1 u1983390 u1983390 222K Sep  7 15:22 all_Maei_and_FoEC_comb.REDUNDANT.CLUSTERED_75.Alignment.hit_table.txt

-rw-r--r-- 1 u1983390 u1983390 2.6M Sep  7 15:22 all_Maei_and_FoEC_comb.REDUNDANT.CLUSTERED_75.Alignment.txt

-rw-r--r-- 1 u1983390 u1983390 3.2M Sep  7 15:22 all_Maei_and_FoEC_comb.REDUNDANT.CLUSTERED_75.Alignment.xml

-rw-r--r-- 1 u1983390 u1983390  45K Sep  7 15:22 Maei_and_FoEC_EffectorPs_comb.CLUSTERED.fasta_17.Alignment.hit_table.txt

-rw-r--r-- 1 u1983390 u1983390 466K Sep  7 15:22 Maei_and_FoEC_EffectorPs_comb.CLUSTERED.fasta_17.Alignment.txt

-rw-r--r-- 1 u1983390 u1983390 602K Sep  7 15:22 Maei_and_FoEC_EffectorPs_comb.CLUSTERED.fasta_17.Alignment.xml

-rw-r--r-- 1 u1983390 u1983390 132K Sep  7 15:22 all_Maei_and_FoEC_comb.NON-EFFECTORS.CLUSTERED_59.Alignment.hit_table.txt

-rw-r--r-- 1 u1983390 u1983390 1.8M Sep  7 15:22 all_Maei_and_FoEC_comb.NON-EFFECTORS.CLUSTERED_59.Alignmnent.txt

-rw-r--r-- 1 u1983390 u1983390 2.1M Sep  7 15:22 all_Maei_and_FoEC_comb.NON-EFFECTORS.CLUSTERED_59.Alignment.xml
