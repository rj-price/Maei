# Running Maei to identify candidate effectors in *F. oxysporum*

All scripts referenced and used in this doc can be found in the [`bin`](https://github.com/JamiePike/Maei/tree/main/bin) directory for this repository.

## *F. oxysporum* genome database

| Species                        | Race | Isolate ID    | Accession           | Size (Mb) | Contig No. | Contig N50 (Mb) | GC (%) | Completeness (% complete BUSCOs) |
|--------------------------------|------|---------------|---------------------|-----------|------------|-----------------|--------|-----------------------------------|
| *F. graminearum*               |      | PH-1          | GCA_000240135.3     | 38        | 5          | 9.3             | 48.2   | --                                |
| *Fo.* fsp. *apii*               | 2    | AJ720[^a]      | --                  | 64.6      | 29         | 4.1             | 47.9   | 97.3                             |
| *Fo.* fsp. *apii*               | 2    | 207.A         | GCA_014843455.1     | 64.7      | 49         | 3.5             | 47.5   | 98.7                             |
| *Fo.* fsp. *apii*               | 3    | NRRL38295     | GCA_014843565.1     | 65.3      | 75         | 4               | 47.5   | 98.8                             |
| *Fo.* fsp. *apii*               | 4    | AJ498[^a]      | --                  | 64.6      | 58         | 2.4             | 47.8   | 96.2                             |
| *Fo.* fsp. *apii*               | 4    | 274.AC        | GCA_014843555.1     | 67.3      | 114        | 4.4             | 47.5   | 98.8                             |
| *Fo.* fsp. *capsici*            |      | Focpep1       | GCA_016801315.1     | 54.5      | 34         | 5               | 47.5   | 93.5                             |
| *Fo.* fsp. *cepae*              | 2    | FoC_Fus2      | GCA_003615085.1     | 53.4      | 34         | 4.1             | 47.5   | 99                               |
| *Fo.* fsp. *conglutinans*       |      | Fo5176        | GCA_014154955.1     | 68        | 25         | 3.4             | 48     | 99.1                             |
| *Fo.* fsp. *coriandrii*         |      | AJ615[^a]      | --                  | 69.3      | 45         | 3               | 48     | 97.5                             |
| *Fo.* fsp. *coriandrii*         |      | 3-2           | GCA_014843415.1     | 65.4      | 49         | 5               | 47.5   | 98.7                             |
| *Fo.* fsp. *coriandrii*         |      | GL306         | GCA_014843445.1     | 65        | 50         | 4.9             | 47.5   | 98.8                             |
| *Fo.* fsp. *cubense*            | 1    | 160527        | GCA_005930515.1     | 51.1      | 12         | 4.9             | 47     | 99.1                             |
| *Fo.* fsp. *cubense*            | 1    | 60            | GWHAAST00000000     | 48.6      | 35         | 2.1             | 47.6   | 95.2                             |
| *Fo.* fsp. *cubense*            | 1    | N2            | GCA_000350345.1     | 47.7      | 2,185      | 0.1             | 48     | --                               |
| *Fo.* fsp. *cubense*            | 4    | C1HIR_9889    | GCA_001696625.1     | 46.7      | 1,318      | 0.09            | 48.5   | --                               |
| *Fo.* fsp. *cubense*            | 4    | B2            | GCA_000350365.1     | 52.9      | 3,834      | 0.02            | 48     | --                               |
| *Fo.* fsp. *cubense*            | TR4  | UK0001        | GCA_007994515.1     | 48.6      | 15         | 4.5             | 47.5   | 98.4                             |
| *Fo.* fsp. *cubense*            | TR4  | 58            | GWHAASU00000000     | 48.2      | 29         | 4.4             | 47.5   | 96.9                             |
| *Fo.* fsp. *cubense*            | TR4  | Pers4         | GCA_021237285.1     | 46.4      | 115        | 1.6             | 47.5   | 97.7                             |
| *Fo.* fsp. *cubense*            | TR4  | II5_54006    | GCA_000260195.2     | 46.6      | 716        | 0.3             | 47.5   | --                               |
| *Fo.* fsp. *cubense*            |      | VPRI44082     | GCA_025216935.1     | 46.3      | 666        | 0.3             | 47     | --                               |
| *Fo.* fsp. *cubense*            |      | VPRI44083     | GCA_025216865.1     | 46.3      | 666        | 0.3             | 47     | --                               |
| *Fo.* fsp. *cubense*            |      | VPRI44081     | GCA_025216985.1     | 47.2      | 902        | 0.4             | 47     | --                               |
| *Fo.* fsp. *cubense*            |      | VPRI44079     | GCA_025216905.1     | 49.5      | 1,801      | 0.3             | 47.5   | --                               |
| *Fo.* fsp. *cubense*            |      | VPRI44084     | GCA_025216845.1     | 50.2      | 2,807      | 0.3             | 47.5   | --                               |
| *Fo.* endophyte                 | np   | Fo47[^d]      | GCA_013085055.1     | 50.4      | 12         | 4.5             | 47.5   | 99                               |
| *Fo.* fsp. *lactucae*           | 1    | AJ865[^a]     | --                  | 62.7      | 38         | 2.7             | 47.7   | 95.3                             |
| *Fo.* fsp. *lactucae*           | 1    | AJ718[^a]     | --                  | 62.1      | 39         | 2.5             | 47.6   | 95.6                             |
| *Fo.* fsp. *lactucae*           | 1    | AJ520[^a][^d]   | --                  | 62.2      | 40         | 2.9             | 47.6   | 95.1                             |
| *Fo.* fsp. *lactucae*           | 4    | AJ705[^a][^d]   | --                  | 66.2      | 32         | 3               | 47.7   | 97.7                             |
| *Fo.* fsp. *lactucae*           | 4    | AJ592[^a]     | --                  | 66        | 36         | 2.6             | 47.7   | 97.5                             |
| *Fo.* fsp. *lactucae*           | 4    | AJ516[^a][^d]   | --                  | 68.8      | 37         | 3               | 47.6   | 97.6                             |
| *Fo.* fsp. *lini*               |      | 39            | GCA_012026625.1     | 59.2      | 34         | 3.4             | 47.5   | 99.5                             |
| *Fo.* fsp. *lycopersici*        | 2    | 4287          | GCA_001703175.2     | 56.2      | 47         | 4.1             | 47.5   | 99.5                             |
| *Fo.* fsp. *matthiolae*         |      | AJ260[^d]     | GCA_020796175.1     | 60.3      | 40         | 4.5             | 47.5   | 97.8                             |
| *Fo.* fsp. *matthiolae*         |      | PHW726_1      | GCA_009755825.1     | 57.2      | 585        | 0.7             | 47     | --                               |
| *Fo.* fsp. *narcissi*           |      | FON63[^a]     | --                  | 60        | 34         | 4               | 47.9   | --                               |
| *Fo.* fsp. *niveum*             | np   | 110407-3-1-1  | GCA_019593455.1     | 49.7      | 33         | 2.8             | 47     | 99.8                             |
| *Fo.* fsp. *rapae*              |      | Tf1208        | GCA_019157295.1     | 59.8      | 25         | 4.2             | 47.5   | 99                               |
| *Fo.* from rocket               |      | AJ174[^a]     | --                  | 62.6      | 30         | 2.7             | 47.9   | 97.8                             |
| *Fo.* fsp. *vasinfectum*        | 1    | TF1           | GCA_009602505.1     | 50        | 17         | 4.2             | 47     | 98.8                             |
| *F. sacchari*                  | 1    | FS66          | GCA_017165645.1     | 47.5      | 47         | 2               | 48     | 99.5                             |
| *F. sacchari*                  |      | NRRL 66326    | GCA_013759005.1     | 42.8      | 515        | 0.2             | 49     | --                               |
| *Fusarium*[^c]                 |      | S6[^b][^e]      | --                  | 47.2      | 6,048      | 0.07            | 47.9   | 97.5                             |
| *Fusarium*[^c]                 |      | S16[^b][^e]     | --                  | 44.9      | 768        | 0.2             | 47.6   | 97.4                             |
| *Fusarium*[^c]                 |      | S32[^b][^e]     | --                  | 40.9      | 2,443      | 0.09            | 48.9   | 97.4                             |
| *Fusarium*[^c]                 |      | SY-2[^b]      | --                  | 44.2      | 408        | 0.2             | 48     | 99.6                             |

*Unfortunately, md format puts footnotes at the bottom of the doc. They are hyperlinked, however, so you can jump between here and the bottom of the doc by clicking them.*

## Maei Set up

Once the assembly database (db) had been generated, genome assemblies were acquired from NCBI (manual download from the assembly pages: https://www.ncbi.nlm.nih.gov/assembly/) or NDCG. Additonal assemblies were provided by collaborators at NAIB.

The assemblies were place in a folder ```/Maei/exp/AllFusAnalysis-Chap3/Maei```

As most of my scripts are designed to work with the `.fna` extension, and as it makes it easier for one-liners downstream, I moved any files ending in `.fasta` to `.fna` using `mv`.

The files were then listed for input into the Maei pipeline

```bash
# create reference fasta list.
ls *.fna > FastaList.txt
```

I need to remove any spaces in the FASTA headers to prevent truncated headers which can cause issues when using certain packages/software. I used the custom bash script `Space_for_underscore.sh` for this.

```bash
# loop through the genome assemblies and replace spaces for underscores. 
for i in $(cat FastaList.txt); do 
Space_for_underscore.sh -i ${i} ; done 

# this script creates a backup ".bak" file for each genome, where spaces have not been replace. I just deleted these files using rm.
rm *.bak
```

I then copied the [reference *mimp* HMM](https://github.com/JamiePike/Maei/blob/main/dev/notes/Mimp_Searching_Methods.md) to the working dir, and ran the `Maei_v5.sh` pipeline.

```bash
# copy hmm built earlier 
cp /[LOCAL PATH TO]/Maei/data/Example_Working_Dir/realn_All.aln.trimmed.hmm ./Mimp.hmm

# run maei 
nohup Maei_V5.sh FastaList.txt 0.65 1>Maei.log &
```

After 8 day, Maei was finished running.

## Filtering the candidates

Now I had a set of candidate effectors for each genome assembly and could cluster the sequences to try and identify shared candidates, but there may be sequences which have been missed using the pipeline in some assemblies, and we need to check to see if they have been missed before developing diagnostics or conducting further analysis. Therefore, I generated a non-redundant set of all candidate effectors to search back against all the genome assemblies using BLAST.

The pipeline generates an output file `./AllCandidateEffectorSets.fasta`, which combines all of the individual candidate effector FASTA files from each genome assembly into one FASTA file for all genomes.

```bash
# extract of lines 192-198 in Maei_v5.sh
touch ./AllCandidateEffectorSets.fasta #Create an empty file to hold all of the candidate effectors.
echo "Clustering final effector sets..."
for i in $(cat ${gL}); 
  do
  sed -i "s/^>/>${i};;_/" "$i" ${i}_Results/${i}_CandidateEffector.fasta ; #Add genome filename to start of Fasta Headers so we know which isolate this came from. 
  cat ${i}_Results/${i}_CandidateEffector.fasta >> ./AllCandidateEffectorSets.fasta ; #Combine all of the individual candidate effector sets by adding them to the empty candidate effector FASTA.
done
```

I used this FASTA and cd-hit to generate the non-redundant set of candidate effectors, where the `-i` flag (for % identity) is set to `1` (100% identity).

```bash
# use cd-hit to generate a non-redundant set of candiate effectors.
cd-hit -i ./AllCandidateEffectorSets.fasta -d 0 -o ./AllCandidateEffectorSets-nonredundant -c 1 -n 5  -G 1 -g 0 -b 20 -l 10 -s 0.0 -aL 0.0 -aS 0.0 1> ./cd-hit_nr.log
```

Now I have a non-redundant set of candidate effectors from all of the genome assemblies I used, I can BLAST them against the input genomes and find candidates which may have been initially missed. For this, I use my `BigBlast.sh` script. I used tBLASTn and a 65% identity and coverage threshold. E-value is 1e-6 as default.

```bash
# blast search the non-redundant set of candidate effectors against all of the genomes used with tblastn. 
BigBlast.sh -q AllCandidateEffectorSets-nonredundant.fasta -g FastaList.txt -i 65 -c 65 -t tblastn -b

# move the file out of the current dir 
mv BigBlast_tblastn_29-01-23 ../BlastOfNRsetAgaintsGenomes
```

Now I have my BLAST results, but I do not necessarily have candidate effectors. Using this BLAST searching method, I often find small hits for sequences that are unlikely to be candidate effectors. To refine the set and increase confidence in my candidates, I filter the BLAST results using SignalP (v5.0b) and EffectorP (v2.0.1).

I have another shell script for this, `FilterBlastHits2.sh`. This extracts the sequences from the BLAST hit locations, translates them, then passes them to signalP (v5.0b) followed by EffectorP (v2.0.1).

```bash
# filter blast hits 
FilterBlastResults.sh FastaList.txt
```

The output of this shell script is an FASTA file containing all the filtered BLAST results (candidate effectors). I use this to quantify they total number of candidate effectors identified per genome assembly.

```bash
#command to count candidate effectors
grep -c ">" */*-EffectorP.filtered.fasta
```

| Species | fsp  | race | isolate code | No. mimps | No. candidate effectors  | Genome size (Mb) |
|---|---|---|---|---|---|---|
| *graminearum* |  |  | PH-1 |  | 6 | 38 |
| *oxysporum* |*apii*| 2 | 207A | 420 | 388 | 64.6 |
| *oxysporum* |*apii*| 2 | AJ720 | 442 | 399 | 64.7 |
| *oxysporum* |*apii*| 3 | NRRL38295 |  | 328 | 65.3 |
| *oxysporum* |*apii*| 4 | AJ498 | 539 | 332 | 64.6 |
| *oxysporum* |*apii*| 4 | 274AC |  | 357 | 67.3 |
| *oxysporum* |*cepae*| 2 | FoC_Fuc2 | 325 | 359 | 53.4 |
| *oxysporum* |*conglutinans*|  | Fo5176 | 442 | 385 | 68 |
| *oxysporum* |*coriandrii*|  | 3-2 | 478 | 315 | 65.4 |
| *oxysporum* |*coriandrii*|  | AJ615 | 675 | 603 | 69.3 |
| *oxysporum* |*cubense*| 1 | 60 | 39 | 45 | 48.6 |
| *oxysporum* |*cubense*| 1 | N2 | 141 | 63 | 47.7 |
| *oxysporum* |*cubense*| 1 | 160527 | 183 | 95 | 51.1 |
| *oxysporum* |*cubense*| 4 | B2 | 142 | 60 | 52.9 |
| *oxysporum* |*cubense*| 4 | C1HIR_9889 | 145 | 70 | 46.7 |
| *oxysporum* |*cubense*| TR4 | Pers4 | 87 | 62 | 46.4 |
| *oxysporum* |*cubense*| TR4 | NRRL_54006 | 105 | 70 | 46.6 |
| *oxysporum* |*cubense*| TR4 | 58 | 165 | 108 | 48.2 |
| *oxysporum* |*cubense*| TR4 | UK0001 | 160 | 127 | 48.6 |
| *oxysporum* |*cubense*|  | VPRI44081 | 137 | 45 | 47.2 |
| *oxysporum* |*cubense*|  | VPRI44082 | 146 | 50 | 46.3 |
| *oxysporum* |*cubense*|  | VPRI44083 | 146 | 50 | 46.3 |
| *oxysporum* |*cubense*|  | VPRI44084 | 179 | 104 | 49.5 |
| *oxysporum* |*cubense*|  | VPRI44079 | 149 | 107 | 49.5 |
| *oxysporum* | endophyte | np | Fo47 | 70 | 81 | 50.4 |
| *oxysporum* | from rocket |  | AJ174 | 419 | 169 | 62.6 |
| *oxysporum* |*lactucae*| 1 | AJ718 | 536 | 260 | 62.1 |
| *oxysporum* |*lactucae*| 1 | AJ865 | 569 | 296 | 62.7 |
| *oxysporum* |*lactucae*| 1 | AJ520 | 533 | 343 | 62.2 |
| *oxysporum* |*lactucae*| 4 | AJ592 | 615 | 482 | 66 |
| *oxysporum* |*lactucae*| 4 | AJ705 | 614 | 490 | 66.2 |
| *oxysporum* |*lactucae*| 4 | AJ516 | 522 | 548 | 68.8 |
| *oxysporum* | *lini* |  | 39 | 263 | 237 | 59.2 |
| *oxysporum* | *lycopersici* | 2 | 4287 | 332 | 337 | 56.2 |
| *oxysporum* | *matthiolae* |  | AJ260 | 301 | 209 | 60.3 |
| *oxysporum* | *narcissus* |  | FON63 | 555 | 251 | 60 |
| *oxysporum* | *niveum* | np | 110407-3-1-1 | 52 | 39 | 49.7 |
| *oxysporum* | *rapae* |  | Tf1208 | 377 | 278 | 59.8 |
| *oxysporum* | *vasinfectum* | 1 | TF1 | 199 | 91 | 50 |
| *sacchari* |  |  | FS66 | 9 | 12 | 47.5 |
| *sacchari* |  |  | NRRL_66326 | 10 | 15 | 42.8 |
| Unconfirmed |  |  | SY2 | 3 | 12 | 44.2 |

## Clustering of candidate effectors

Now I have my candidate effectors for each isolate, I need to determine the candidate effector profile of all isolates; that is cluster the sequences at a given %ID and see what is shared or unique. First, I must extract all of the candidate effectors into a single FASTA for cd-hit.

```bash
#copy the candidate effector fasta files into the current directory so that the full dir name isn't added to the start of each sequence.
cp */*-EffectorP.filtered.fasta ./

#use awk to add the fasta header to the start of each sequence so we know which isolate it came from and generate a final file. 
for i in *-EffectorP.filtered.fasta ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' ${i} ; done > AllCandidateEffectors-AfterBLASTandFilter.fasta

#I can now remove the copied effector fasta files. This will save space. 
rm *-EffectorP.filtered.fasta
```

I can now use my `Maei_Clustering.sh` shell script to cluster the candidate effectors at various %ID and see how this effects the candidate effector profiles.

```bash

```

### Table footnotes

[^a]: Assemblies generated by collaborators at NIAB. GenBank accession numbers are not currently available.
[^b]: Assemblies generated in association with TNAU. GenBank accession numbers are not currently available.
[^c]: Species not confirmed.
[^d]: Included in RNA-seq analysis.
[^e]: Not included in MAEI pipeline analysis.
