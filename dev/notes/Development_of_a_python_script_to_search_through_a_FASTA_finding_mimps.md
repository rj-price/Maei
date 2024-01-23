# Development of a python script to search through a FASTA finding mimps

`mimp_finder`

A script was developed to search through a FASTA file and generate a bed file containing the location of mimps.

This script (mimp_finditer.py) can be found:

Vettel:

/home/u1983390/Scripts/Mimp_finditer.py 

Shared:

/Volumes/shared/HCSS3/Shared258/Jamie/Scripts/[Mimp_finditer.py](http://mimp_finditer.py/) 

On Github (Private):

[https://github.com/JamiePike/Effector_Pred_Maei/blob/master/mimp_finditer.py](https://github.com/JamiePike/Effector_Pred_Maei/blob/master/mimp_finditer.py)

Example command:

(biopythonEnv) u1983390: python ~/Scripts/[Mimp_finditer.py](http://mimp_finditer.py/) [FASTA INFILE] 

Created in Atom. Biopython required. 

Default python version becomes 3.8.3 when this condo environment is activated. 

---

MIMP FINDER - SEARCHES FOR ALL MIMP IN SEQUENCE

Set up and import correct modules.

```
    from datetime import datetime
    startTime = datetime.now()
    import re, sys

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
```

Establish input file (FASTA) and start time.

    

```
infile = sys.argv[1]
    print(f'\n\nRunning Mimp_finder.py on {infile} at {startTime}.')
```

Open a .bed file using the name of the infile. The loctation of the mimp finder will be desposited within.

```
    with open(f'{infile}_mimp_hits.bed'.format(), "a") as file_object:
```

Parse and loop through the FASTA input file. Search through each scaffold/contig of the FASTA looking for sequences which match the regex/mimp TIR "CAGTGGG..GCAA[TA]AA" (re.IGNORECASE remove case sensitivity issues).

mimp_length = Create a variable for mimp length - Mimps cannot exceed 400 nucleotides.

mimp_count = Create a variable to count the number of mimps found on each contid/scaffold.

```
        for sequence in SeqIO.parse(infile, "fasta"):
            print("="*80 + "\n")
            print(f'Searching {sequence.description} for mimps.\n')
            hit = re.finditer(r"CAGTGGG..GCAA[TA]AA", str(sequence.seq), re.IGNORECASE)
            mimp_length = 400
            mimp_count = 0
```

For every instance where the mimp TIR "CAGTGGG..GCAA[TA]AA" is found, print the location of the hit and search that contig for the other mimp TIR.

h_start = record the start location of the mimp

```
            for h in hit:
                if h:
                    h_start = h.start()
                    print(f'--Mimp TIR found: "{h.group()}" at position {h.start()} to {h.end()}--')
                    hit_rc = re.finditer(r"TT[TA]TTGC..CCCACTG", str(sequence.seq), re.IGNORECASE)
```

If the an other mimp TIR is found, ensure that the second mimp TIR does not come before the first mimp TIR.

Then, if TRUE, ensure that the distance between the start of the first mimp TIR and end of the second mimp TIR is not greater than mimp_length (400 nucleotides).

If TRUE, add 1 to the mimp_count variable and print the contig/scaffold description along with the start location of the mimp and the end location of the mimp in the .bed file opened (ensure fields are \t seperated to maintain .bed file convention).

If the distance between the start of the first mimp TIR and the end of the second mimp TIR is greater than 400 nucleotides, it is not a mimp.

If the second mimp TIR comes before the first mimp TIR found (second mimp TIR end - first mimp TIR start < 0), it is not a mimp.

```
                    for h_rc in hit_rc:
                        print('Looking for reverse complement sequence(s)...')
                        h_rc_end = h_rc.end()
                        if h_rc:
                            print(f'Mimp reverse complement TIR found: "{h_rc.group()}" at position {h_rc.start()} to {h_rc.end()}.  Is this complete mimp?')
                            length = h_rc_end - h_start
                            if length > 0:
                                if length < mimp_length:
                                    print(f'Yes! Full mimp found at position {h.start()} to {h_rc.end()}.')
                                    print("The mimp length is:")
                                    print(f'{length}bp\n')
                                    mimp_count = mimp_count+1
                                    print(sequence.description, h.start(), h_rc.end(), file=file_object, sep="\t")
                                else:
                                    print("No, TIRs are at a greater distance than 400bp. \nOnly a partial mimp.\n")
                            else:
                                print("No, the TIR's reverse complement is from the previous mimp as it comes before the start location of this TIR\n")
```

Print the number of mimps found on that contig/scaffold.

When whole FASTA has been looped through, print the file searched and time finished.

```
                    print(f'\n{sequence.description} contains {mimp_count} mimp(s)')

        print("="*80 + "\n")

        print(f'\nMimp_finder.py complete. \n\nFile searched: {infile}\n\nFinished at: {startTime}')
```

---

**Example Output:**

Printed to the the screen: 

(biopythonEnv) MACP80255:cubense_mimps_public u1983390-> python ~/Scripts/Mimp_finditer.py Focub_all_Chang_mimps.faa 

Running [Mimp_finder.py](http://mimp_finder.py/) on Focub_all_Chang_mimps.faa at 2020-09-08 15:33:34.071948.

================================================================================

Searching Focub_II5_mimp_1 _contig_1.16(656599:656809) for mimps.

--Mimp TIR found: "CAGTGGGATGCAAAAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTTTTGCATCCCACTG" at position 192 to 208.  Is this complete mimp?

Yes! Full mimp found at position 2 to 208.

The mimp length is:

206bp

Focub_II5_mimp_1 _contig_1.16(656599:656809) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_2 _contig_1.47(41315:41540) for mimps.

--Mimp TIR found: "CAGTGGGAGGCAATAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTATTGCCCCCCACTG" at position 207 to 223.  Is this complete mimp?

Yes! Full mimp found at position 2 to 223.

The mimp length is:

221bp

Focub_II5_mimp_2 _contig_1.47(41315:41540) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_3 _contig_1.65(13656:13882) for mimps.

--Mimp TIR found: "CAGTGGGAGGCAATAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTATTGCCCCCCACTG" at position 208 to 224.  Is this complete mimp?

Yes! Full mimp found at position 2 to 224.

The mimp length is:

222bp

Focub_II5_mimp_3 _contig_1.65(13656:13882) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_4 _contig_1.70(61591:61809) for mimps.

--Mimp TIR found: "CAGTGGGATGCAATAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTATTGCATCCCACTG" at position 200 to 216.  Is this complete mimp?

Yes! Full mimp found at position 2 to 216.

The mimp length is:

214bp

Focub_II5_mimp_4 _contig_1.70(61591:61809) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_5 _contig_1.70(63407:63595) for mimps.

--Mimp TIR found: "CAGTGGGATGCAAAAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTTTTGCACCCCACTG" at position 170 to 186.  Is this complete mimp?

Yes! Full mimp found at position 2 to 186.

The mimp length is:

184bp

Focub_II5_mimp_5 _contig_1.70(63407:63595) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_6 _contig_1.86(922:1140) for mimps.

--Mimp TIR found: "CAGTGGGATGCAATAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTATTGCATCCCACTG" at position 200 to 216.  Is this complete mimp?

Yes! Full mimp found at position 2 to 216.

The mimp length is:

214bp

Focub_II5_mimp_6 _contig_1.86(922:1140) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_7 _contig_1.89(690:899) for mimps.

--Mimp TIR found: "CAGTGGGATGCAAAAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTTTTGCATCCCACTG" at position 191 to 207.  Is this complete mimp?

Yes! Full mimp found at position 2 to 207.

The mimp length is:

205bp

Focub_II5_mimp_7 _contig_1.89(690:899) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_8 _contig_1.94(19371:19577) for mimps.

--Mimp TIR found: "CAGTGGGATGCAAAAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTTTTGCATCCCACTG" at position 188 to 204.  Is this complete mimp?

Yes! Full mimp found at position 2 to 204.

The mimp length is:

202bp

Focub_II5_mimp_8 _contig_1.94(19371:19577) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_9 _contig_1.107(11011:11230) for mimps.

================================================================================

Searching Focub_II5_mimp_10 _contig_1.107(11270:11516) for mimps.

--Mimp TIR found: "CAGTGGGGGGCAAAAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTTTTGCATCCCACTG" at position 228 to 244.  Is this complete mimp?

Yes! Full mimp found at position 2 to 244.

The mimp length is:

242bp

Focub_II5_mimp_10 _contig_1.107(11270:11516) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_11 _contig_1.128(3654:3884) for mimps.

--Mimp TIR found: "CAGTGGGGGGCAATAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTATTGCCCCCCACTG" at position 212 to 228.  Is this complete mimp?

Yes! Full mimp found at position 2 to 228.

The mimp length is:

226bp

Focub_II5_mimp_11 _contig_1.128(3654:3884) contains 1 mimp(s)

================================================================================

Searching Focub_II5_mimp_12 _contig_1.194(3335:3541) for mimps.

--Mimp TIR found: "CAGTGGGATGCAAAAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTTTTGCATCCCACTG" at position 188 to 204.  Is this complete mimp?

Yes! Full mimp found at position 2 to 204.

The mimp length is:

202bp

Focub_II5_mimp_12 _contig_1.194(3335:3541) contains 1 mimp(s)

================================================================================

Searching Focub_B2_mimp_1 _contig_44(75544:75770) for mimps.

--Mimp TIR found: "CAGTGGGAGGCAATAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTATTGCCCCCCACTG" at position 208 to 224.  Is this complete mimp?

Yes! Full mimp found at position 2 to 224.

The mimp length is:

222bp

Focub_B2_mimp_1 _contig_44(75544:75770) contains 1 mimp(s)

================================================================================

Searching Focub_B2_mimp_2 _contig_135(1675161:1675386) for mimps.

--Mimp TIR found: "CAGTGGGGGGCAATAA" at position 2 to 18--

Looking for reverse complement sequence(s)...

Mimp reverse complement TIR found: "TTATTGCCTCCCACTG" at position 207 to 223.  Is this complete mimp?

Yes! Full mimp found at position 2 to 223.

The mimp length is:

221bp

Focub_B2_mimp_2 _contig_135(1675161:1675386) contains 1 mimp(s)

================================================================================

[Mimp_finder.py](http://mimp_finder.py/) complete. 

File searched: Focub_all_Chang_mimps.faa

Finished at: 2020-09-08 15:33:34.071948

.bed file created:

(biopythonEnv) MACP80255:cubense_mimps_public u1983390-> cat Focub_all_Chang_mimps.faa_mimp_hits.bed 

Focub_II5_mimp_1 _contig_1.16(656599:656809)2208

Focub_II5_mimp_2 _contig_1.47(41315:41540)2223

Focub_II5_mimp_3 _contig_1.65(13656:13882)2224

Focub_II5_mimp_4 _contig_1.70(61591:61809)2216

Focub_II5_mimp_5 _contig_1.70(63407:63595)2186

Focub_II5_mimp_6 _contig_1.86(922:1140)2216

Focub_II5_mimp_7 _contig_1.89(690:899)2207

Focub_II5_mimp_8 _contig_1.94(19371:19577)2204

Focub_II5_mimp_10 _contig_1.107(11270:11516)2244

Focub_II5_mimp_11 _contig_1.128(3654:3884)2228

Focub_II5_mimp_12 _contig_1.194(3335:3541)2204

Focub_B2_mimp_1 _contig_44(75544:75770)2224

Focub_B2_mimp_2 _contig_135(1675161:1675386)2223

Note it does not find a mimp on mimp_9 as the TIR deviates from the TIR/regex in the script.

Mimp_9:

```
>Focub_II5_mimp_9 _contig_1.107(11011:11230)
TTCAGTGAGGTGCAATAAGTTTGAATTTGAGATGAAGTATCTGCCTTGTCGCTCTAGTTCTGTCAGAATGGTGACTAGCTATATCGGGTAGTTATCTAGCCTAGGTAGGCAACTAGCTAACTTGATGAGGCACAAGCAGAGCAAAGTGATAGCTACCGATCGGCTTCGCAAGGAGGCGGGTACTTGCGCCTGGATTCAAACTTCTTGCATCTCACTGTA
```

---

Script:

[Mimp_finditer.py](./file/Mimp_finditer.py)

 
