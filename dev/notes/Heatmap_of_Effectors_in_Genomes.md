# Heatmap of Effectors in Genomes

Using R, a clustered heatmap was built do identify if candidate effectors are specific to race. 

[Transposed_heatmap_for_candid_effectors-nb.html](./file/Transposed_heatmap_for_candid_effectors-nb.html)

[Transposed_heatmap_for_candid_effectors.Rmd](./file/Transposed_heatmap_for_candid_effectors.Rmd)

The following script was used: 

```
/Users/u1983390/Documents/Lab\ Book/Candidate\ Effectors/Transposed_heatmap_for_candid_effectors.Rmd 
```

This csv was used: 

```
/Users/u1983390/Documents/Lab\ Book/Candidate\ Effectors/transposed_list_of_effectors_in_genomes.csv
```

The CSV was generated using the excel doc (

```
/Users/u1983390/Documents/Lab\ Book/Candidate\ Effectors/List\ of\ effectors\ and\ the\ genomes\ they\ belong\ to\ .xlsx
```

) prepared in Effector Set Analysis. 
![image.png](image/image.png)

A new fasta file was generated which contained the corresponding Card_Eff_X naming: 

```
Cand_Eff.fasta
```

It was done using the python script below in the following directory on Vettel: `/home/u1983390/Fusarium_data/EFFECTOR_SET_ANALYSIS`

```
#!/usr/bin/env python

fasta= open('All_effector_sets_combined.non_redundant.fsa_aa')
newnames= open('Cand_Eff_header_list.txt')
newfasta= open('Cand_Eff.fasta', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
```
