# EffectorP Random Sample Test

Following on from the meeting with David, Dan, and Murray (see Meetings: Project  Update - 25/09/2020) I tested a random sample of predicated proteins from the genomes listed below. The predicted proteins were randomly sampled using reformat.sh from bbmap and filtered to include a signal; peptide using SignalP. They were then combined into one FASTA. (See script below). EffectorP 2.0.1 was run on the combine FASTA and the Summary format and a predicted effector FASTA were generated.  

Predicted Proteomes used: 

GCA_000350365.1_Foc4_1.0_B2

GWHAASU00000000_FocTR4_58

GCA_007994515.1_UK0001

GCA_000350345.1_Foc1_1.0_N2

GCA_000260195.2_FO_II5_V1

GWHAAST00000000_Foc1_60

GCA_005930515.1_160527

Script used: Can be found on Vettel /home/u1983390/Fusarium_data/EFFECTORP_RANDOM_SAMPLE_TEST/[Extract_random_signal_pep_sequences.sh](http://Extract_random_signal_pep_sequences.sh) 

```
#!/bin/bash 

#Created 29/09/2020 by Jamie Pike 

#Extract a random sample, identify if it has a signal peptide and generate a FASTA containing those sequences. 

#SignalP path

Run_SignalP=/home/u1983390/apps/FoEC/signalp-4.1/signalp                         

#Create for loop to loop though genomes with annotations/protein sequences 

for i in $(cat File_list.txt); do 

  #Extract random sample using bbmap's [reformat.sh](http://reformat.sh) (idea from [https://www.biostars.org/p/241972/](https://www.biostars.org/p/241972/)). 

  [reformat.sh](http://reformat.sh) in=${i}.protein.faa out=stdout.fa samplereadstarget=200 ignorejunk | sed 's/ /_/g' > ${i}_random_prot_sample.faa ;

  #Identify sequences with a signal peptide, and generate a list of the FASTA headers for these sequences. 

  $Run_SignalP -f summary ${i}_random_prot_sample.faa | grep "^Name*" | awk '$2 ~ /YES/ {print $0}' | awk  '{gsub("=","\t",$0); print;}' | awk '{print $2}' > ${i}_samtools.list ;

  #Extract these sequences using samtools faidx.

  for j in $(cat ${i}_samtools.list); do samtools faidx ${i}_random_prot_sample.faa ${j} ; done > ${i}_random_prot_sample_with_signal_peptide.faa; done 

#Create a masta FASTA file with all of these sequences. 

for i in $(cat File_list.txt); do cat ${i}_random_prot_sample_with_signal_peptide.faa; done > all_random_prot_sample_with_signal_peptide.faa
```

EffectorP outputs (Vettel):

/home/u1983390/Fusarium_data/EFFECTORP_RANDOM_SAMPLE_TEST/all_random_prot_sample_with_signal_peptide.EffectorP_summary

/home/u1983390/Fusarium_data/EFFECTORP_RANDOM_SAMPLE_TEST/all_random_prot_sample_with_signal_peptide.EffectorP.fsa_aa 

Table of results comparing the different EffectorP runs:

(Vettel) /home/u1983390/Fusarium_data/EFFECTORP_RANDOM_SAMPLE_TEST/Table_of_EffectorP_pred.txt

Run          Date        Numbers  Percentage_pred

Maei         13/08/2020  14/62    23%

Maei         21/09/2020  15/64    23%

[FoEC.py](http://FoEC.py)      04/08/2020  9/28     32%

Rand_sample  30/09/2020  29/150   19%

I also checked the paper with which EffectorP is associated (Link:[https://bsppjournals.onlinelibrary.wiley.com/doi/full/10.1111/mpp.12682](https://bsppjournals.onlinelibrary.wiley.com/doi/full/10.1111/mpp.12682)). They don’t give a clear/strict/decisive false positive or sensitivity indication. They do present a table in which they have tried different collections of proteins against different classifiers (including EffectorP 2.0) which is copied below:

**Table 3.**Independent validation of EffectorP's prediction accuracy.

|                                                                                            |             |                       |                |                     |                     |                               |
|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|-----------------------|----------------|---------------------|---------------------|-------------------------------|
|                                                                                                                                                                        |             |**Predicted effectors**|                |                     |                     |                               |
|Dataset                                                                                                                                                                 |# of proteins|EffectorP 2.0          |EffectorP 1.0   |EffectorP 1.0 and 2.0|Small size classifier|Small, cysteine‐rich classifier|
|Fungal saprophyte secreted proteins                                                                                                                                     |24 432       |2865 (11.7%)           |4774 (19.5%)    |2444 (10%)           |10 529 (43.1%)       |4961 (20.3%)                   |
|Fungal, plant and mammalian proteins with signal peptide and localization to endoplasmic reticulum, Golgi, membranes or with glycosylphosphatidylinositol (GPI) anchors |2631         |220 (8.4%)             |294 (11.2%)     |164 (6.2%)           |654 (24.9%)          |307 (11.7%)                    |
|Fungal proteins with unaffected pathogenicity phenotype                                                                                                                 |938          |45 (4.8%)              |59 (6.3%)       |36 (3.8%)            |128 (13.6%)          |60 (6.4%)                      |
|                                                                                                                                                                        |**28 001**   |**3130 (11.2%)**       |**5127 (18.3%)**|**2644 (9.4%)**      |**11 311 (40.4%)**   |**5328 (19%)**                 |
|Fungal effector positive training set                                                                                                                                   |94           |89 (94.7%)             |80 (85.1%)      |79 (84%)             |88 (93.6%)           |53 (56.4%)                     |
|Fungal effector independent test set                                                                                                                                    |21           |16 (76.2%)             |16 (76.2%)      |16 (76.2%)           |19 (90.5%)           |10 (47.6%)                     |
|**Accuracy**                                                                                                                                                            |             |88.8%                  |81.7%           |**90.5%**            |59.8%                |80.9%                          |
