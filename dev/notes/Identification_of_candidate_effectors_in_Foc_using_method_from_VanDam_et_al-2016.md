# Identification of candidate effectors in Foc using method from VanDam et al  2016

A tool for effector identification used in the paper: Effector profiles distinguish formae speciales of fusarium oxysporum (Bottom of this page) was downloaded from GitHub: [https://github.com/pvdam3/FoEC](https://github.com/pvdam3/FoEC). Installed in the directory /home/u1983390/apps/FoEC/. 

A default python environment was established as [FoEC.py](http://foec.py/) is written in an older version of python (2.7). Similarly, [FoEC.py](http://foec.py/) must be run using a different version of R from the default version on Vettel. Both the python virtual environment and correct version of R must be enabled by following the text in the file Notes_for_running.txt in the FoEC directory. Due to [FoEC.py](http://foec.py/) being written in ~2016, earlier versions of BLAST and SignalP were installed using miniconda. Paths to these earlier versions can be found in the head (-n 20) of the [FoEC.py](http://foec.py/) file. 

Installation was carried out following the instructions for installation found on the FoEC github page. 

---

_Early outputs 4 genomes included in my work and that from van Dam et al., (2016) differed (Table 1). �_�

|

                                                       |Genomes/Strain/Isolate|Number of Complete mimps |Number of Incomplete mimps |Number of Scaffolds/contigs |
|-------------------------------------------------------|----------------------|-------------------------|---------------------------|----------------------------|
|My Runs of FoEC.py                                     |Fol4287               |0                        |19                         |88                          |
|Focub II5                                              |3                     |11                       |418                        |                            |
|Focub B2                                               |3                     |12                       |840                        |                            |
|Focub N2                                               |0                     |3                        |1341                       |                            |
|Published Van Dam _et al.,_ (2016) supplementary info. |Fol4287               |40                       |78                         |144                         |
|Focub II5                                              |11                    |19                       |418                        |                            |
|Focub B2                                               |2                     |7                        |840                        |                            |
|Focub N2                                               |0                     |3                        |1341                       |                            |

It was determined that the result for Fol4287 changed due to the genomic being updated (note scaffold/contig number).

The [FoEC.py](http://foec.py/) script uses a mimp finder script: [01a.mimpfinder_combine_to_putefflist_MetStop.py](http://01a.mimpfinder_combine_to_putefflist_metstop.py/). The mimp finder within this script is case-sensitive. It only search for capitalised mimp-TIRs. The genomes for Foc had all been downloaded with soft-masking therefore many mimps will have been missed by the script.

The B2 and N2 to isolates used in the van Dam (2016) paper had also been soft-masked, hence the lower number mimps reported for these genomes in the paper. After capitalising the genomes to remove soft-masking, the number of mimps found was more representative of those found in the other genomes used in the van Dam et al (2016) paper. For more information on bug fixing the FoEC mimp finder, see the email chain between Jamie pike and Dr Like Fokkens (subj. Mimp sequences in F. o. cubense) and the directory: 

/home/u1983390/Fusarium_data/FoEC_MIMPS/Like_Fokkens_Bug_Fixing. 

Email Chain: 

[F72AEEF7-6268-4242-B3CB-3CA4C702C4E1.olk15Message](./file/F72AEEF7-6268-4242-B3CB-3CA4C702C4E1.olk15Message)

For [FoEC.py](http://foec.py/) runs on my genomes (Table 2), all genomes were capitalised to remove soft masking in order to allow for the bug discovered in the script from van Dam _et al.,_ (2016). 

|Table 2: Genomes used in the identification of _mimps_ and _SIX_ genes[[1]](applewebdata://90485BA2-28F5-4990-9ABA-F5CCE786D6EA#_ftn1).|                                                                   |              |          |        |                   |                                  |                                                                                                            |
|---------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|--------------|----------|--------|-------------------|----------------------------------|------------------------------------------------------------------------------------------------------------|
|**Assembly Accession**                                                                                                                 |**Bioproject Accession**                                           |**Isolate**   |**Source**|**Race**|**Contig N50 (MB)**|**Predicted Protein Coding genes**|**BUSCO**** – intact single copy orthologs**[[2]](applewebdata://90485BA2-28F5-4990-9ABA-F5CCE786D6EA#_ftn2)|
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

---

[[1]](applewebdata://90485BA2-28F5-4990-9ABA-F5CCE786D6EA#_ftnref1) NCBI refers to genomes from National Center for Biotechnology Information, sourced via GenBank website ([https://www.ncbi.nlm.nih.gov/genbank/](https://www.ncbi.nlm.nih.gov/genbank/)).   NDGC refers to genomes from National Genomics Data Centre website ([https://bigd.big.ac.cn/](https://bigd.big.ac.cn/)).

[[2]](applewebdata://90485BA2-28F5-4990-9ABA-F5CCE786D6EA#_ftnref2)BUSCO - based on genome - % Intact and single copy orthologs (_Sordariomycetes_ data). 

The outputs, log and error files, and command from the final [FoEC.py](http://foec.py/) run can be found on Vettel in: 

Output:

/home/u1983390/Fusarium_data/FoEC_MIMPS/output_20.08.04_21h25m15

.log and .err files:

/home/u1983390/Fusarium_data/FoEC_MIMPS

Command:

/home/u1983390/apps/FoEC/Notes_for_running.txt

[FoEC.py](http://foec.py/) produced a FASTA containing a collection of predicted effectors and a file of those effectors clustered (details on clustering in FoEC README). 

/home/u1983390/Fusarium_data/FoEC_MIMPS/output_20.08.04_21h25m15/02.cluster_putative_effectors/

A heirclust plot was produced by [FoEC.py](http://foec.py/) 

/home/u1983390/Fusarium_data/FoEC_MIMPS/output_20.08.04_21h25m15/04.cluster_and_plot

[hierclust_plot.pdf](./file/hierclust_plot.pdf)

The predicted/candidate effectors from [FoEC.py](http://foec.py/) may allow for race prediction as some effectors cluster within specific races: 

[hierclust_plot_20-08-04_21h51m15_edited.pdf](./file/hierclust_plot_20-08-04_21h51m15_edited.pdf)

Van Dam _ et al.,_ (2016) paper. 

[Van_Dam_Effector_Profiles_distinguish_f-sp-of_Fo.pdf](./file/Van_Dam_Effector_Profiles_distinguish_f-sp-of_Fo.pdf)
