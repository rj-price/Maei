# *Mimp*-associated Effector Identification (Maei)

## Identifying effectors associated with transposable elements in *Fusarium oxysporum* genomes. 

### Pipeline Structure:

Initilally, two methods of *mimp* searching are employed. The first method, uses searching by regular expression, whereby the *mimp* TIR sequences, "CAGTGGG..GCAA[TA]AA" and "TT[TA]TTGC..CCCACTG", are used as a search pattern. Where sequences matching this pattern occur within 400 nucleotides of each other a mimp is recorded. The script for this can be found in the "Supplementary_Scripts" directory as mimp_finditer.py.

The second method, employs a Hidden Markov Model (HMM) which was developed using the HMM tool HMMER (3.3.1) (Eddy, N.D.). Briefly, publicly available *mimp* sequences and *mimp* sequences identified using the regular expression method were used to build a *mimp* profile-HMM. This profile-HMM was is as the input for an NHMMER search of each genome. The profile HMM (Mimp.hmm) should be maintained in the directory in which Maei.sh is to be run, and can be found in the GitHub Example_Directory, alongside two test genomes and a fasta list. 

Using *mimps* identified by both *mimp* searching methods, sequences 2.5kb upstream and downstream of each *mimp* are extracted, mereged and stored in {GenomeName}_mimp_region.gff. This region is also used to generate an Augustus region GFF (*mimp* region + 20kb either side, to prevent truncation of gene models). The GFFSa are used to generate a *mimp* region only fasta and an Augustus region only fasta, whereby all nucleotides outside of the *mimp* or augustus region are hard masked. 

The Augustus region fasta is then subjected to gene prediction using AUGUSTUS (3.3.3) (Stanke, et al., 2006) with the “Fusarium” option enabled, and open reading frames (ORFs) withing the *mimp* regions are identified using the EMBOSS (6.6.0.0) tool, getorf (https://www.bioinformatics.nl/cgi-bin/emboss/getorf) and the *mimp* region fasta. 

The Augustus gene models and *mimp* region gff files are intersected and merged, to ensure that all gene models are associated with a *mimp* region, but that no gene models are truncated. The intersected Augustus and *mimp*-region gff file ({genome file name}_Results/Augustus/{genome file name}.augustus-mimp_intersections.gff) is used as input for AGAT, which extracts Augustus *mimp*-associated gene models in to a fasta file. The *mimp*-associated gene model fasta and getORF fasta are merged to generate a *mimp*-associated gene model and ORF fasta. The output from getORF is also converted to GFF format using the getORF2gff.py script in the Supplementary_Scripts directory. These gffs are then mergered to generate a *mimp*-associated gene model and ORF gff.

The *mimp*-associated gene model and ORF fasta is then submitted to SignalP (4.1), and sequences which are not predicted to conatin a signal peptide are removed. Sequences are then filtered based on size, with sequences >30aa and <450aa parsed to EffectorP (2.0.1) (Sperschneider, et al., 2018) for effector prediction. This output is used to generate a genome sepcific candidate effector fasta and gff file. 

The candidate effector fasta from each genome is then combined and clusterd  using CD-HIT (4.8.1) (Fu, et al., 2012) (percentage identity determined using command line), to generate a candidate effector clusters, whereby differences in effector profile can be determined, and specific effector clusters assessed. Data from CD-HIT is also parsed to a custom python script, Processingcdhit.py, which generates an overview table and a datamatrix which can be used to build a heatmap.

[Maie_v5_Figure.pdf](https://github.com/JamiePike/Maei/files/10708866/Maie_v5_Figure.pdf)

### Dependencies

Maei makes use of multiple bioinformatic tools and python modules, **all of which will have to be configured at the top of the Maei.sh script**.

All dependcies can be installed using Bioconda, and I recommend you create a conda enviroment specifically for Maei. 

The dependencies:
* Bedtools (v2.25.0)
* Samtools (0.1.19-96b5f2294a)
* python3
    -> including biopython and the pandas modules.
* AUGUSTUS (3.3.3) 
* EMBOSS (6.6.6)
* SignalP (4.1)
* EffectorP (2.1)
* HMMER (3.3.1) 
* AGAT (v1.0.0)
* CD-HIT (version 4.6)

To install dependencies (except SignalP and EffectorP) into a new Conda called `maei`:
```bash
# Get pipeline
git clone https://github.com/rj-price/Maei.git 
# Change to directory
cd Maei
# Create conda environment 
conda env create -f environment.yml
# Activate enviroment before running
conda activate maei
```

### Usage 

./Maei_v5.sh [List_of_Genome_FASTAs]  [Sequence identity threshold for cd-hit, e.g. 0.8]

### Notes

Duplicate your orginal running directory before you run Maei.sh. This will save time should you wish to run again. 
