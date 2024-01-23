# Changes to my mimp finder script

`mimp_finder`

Updates can be seen on Github: [https://github.com/JamiePike/Effector_Pred_Maei/commits?author=JamiePike](https://github.com/JamiePike/Effector_Pred_Maei/commits?author=JamiePike)

Following a meeting with Dr Laura Baxter (See: meetings, Update_on_effectors with Laura). Modifications were made to the mimp_finditer.py script. 

Instead of only searching for mimps in one direction:

```
                    
#If the an other mimp TIR is found, ensure that the second mimp TIR does not come before the first mimp TIR.
#Then, if TRUE, ensure that the distance between the start of the first mimp TIR and end of the second mimp TIR is not greater than mimp_length (400 nucleotides).
#If TRUE, add 1 to the mimp_count variable and print the contig/scaffold description along with the start location of the mimp and the end location of the mimp in the .bed file opened (ensure fields are \t seperated to maintain .bed file convention).
#If the distance between the start of the first mimp TIR and the end of the second mimp TIR is greater than 400 nucleotides, it is not a mimp.
#If the second mimp TIR comes before the first mimp TIR found (second mimp TIR end - first mimp TIR start < 0), it is not a mimp.

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

Mimps are now searched for in each direction of the first mimp hit. 

```
#Ensure that the distance between the start of the first mimp TIR and end of the second mimp TIR is not greater than mimp_length (400 nucleotides).
#If TRUE, add 1 to the mimp_count variable and print the contig/scaffold description along with the start location of the mimp and the end location of the mimp in the .bed file opened (ensure fields are \t seperated to maintain .bed file convention).
#If the distance between the start of the first mimp TIR and the end of the second mimp TIR is greater than 400 nucleotides, it is not a mimp.
#If the second mimp TIR comes before the first mimp TIR found (second mimp TIR end - first mimp TIR start < 0), it is not a mimp.

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

**OUTPUTS:**

After running on genome [GCA_007994515.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_007994515.1) (isolate: UK0001, Race: TR4). 

(biopythonEnv) MACP80255:VMNF01.1.fsa_nt_MIMPS u1983390-> wc -l *output.bed

      23 VMNF01.1.fsa_nt_mimp_hits_original_output.bed

      27 VMNF01.1.fsa_nt_mimp_hits_updated_output.bed

(biopythonEnv) MACP80255:VMNF01.1.fsa_nt_MIMPS u1983390-> head -n 27 *output.bed

==> VMNF01.1.fsa_nt_mimp_hits_original_output.bed <==

VMNF01000005.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_5, whole genome shotgun sequence484756484978

VMNF01000007.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_7, whole genome shotgun sequence62944276294648

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence837480837682

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11918651192079

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11918651192223

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11930231193247

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11934631193670

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence27986272798841

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28004432800627

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28723022872509

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28905002890721

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28910332891239

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence29002252900435

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence31922083192434

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence32258773226079

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence35091253509331

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence35854853585690

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36798453680059

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36798453680203

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36810033681227

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36814433681650

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36924703692712

VMNF01000015.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_15, whole genome shotgun sequence17694671769673

==> VMNF01.1.fsa_nt_mimp_hits_updated_output.bed <==

VMNF01000005.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_5, whole genome shotgun sequence484756484978

VMNF01000007.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_7, whole genome shotgun sequence62944276294648

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence837480837682

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11918651192079

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11918651192223

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11930231193247

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11934631193247

VMNF01000013.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_13, whole genome shotgun sequence11934631193670

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence27986272798841

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28004432800627

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28723022872509

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28905002890721

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28910332890721

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence28910332891239

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence29002252900435

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence31922083192434

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence32258773226079

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence35091253509331

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence35854853585690

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence35896963589452

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36798453680059

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36798453680203

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36810033681227

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36814433681227

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36814433681650

VMNF01000014.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_14, whole genome shotgun sequence36924703692712

VMNF01000015.1 Fusarium oxysporum f. sp. cubense strain TR4 isolate UK0001 scf_28419_15, whole genome shotgun sequence17694671769673
