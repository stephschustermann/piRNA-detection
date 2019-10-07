# piRNA-detection

piRNA detection via multidimensional small RNA clustering
Piwi-interacting RNAs (piRNA) are a novel class of small RNA molecules thought to mediate retrotransposon silecing. Since transposable elements (TEs) can serve as mutagenic factors that contribute to genomic instability, piRNA function is extremely important for genome defense.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites - File structures

You will need the following files to use all the features in this project:

1. Transposable Elements file in order to get all the mappings to TE

```
trans_id        trans_coord     target_seq      probe_id        probe_seq       num_mismatch    strand
TAR1    Repeat sequence Homo sapiens    170     GGGAGGAGTAGGCTG HWI-ST1391:140326:C441WACXX:7:2114:14096:41705 1:N:0: CAACTAA   GGGAGGAGTAGGATG 1       -
TAR1    Repeat sequence Homo sapiens    789     GTCCCCGGCGGCGCAGAGACGAGTGGAACCTG        HWI-ST1391:140326:C441WACXX:7:1110:14668:100150 1:N:0: CACCTAA  GTACCCGGCGGCGCAGAGACGAGTGGAACCTG  1       -
MamGypLTR2c_LTR Gypsy   Mammalia        527     TCTTGGGGGGCTCAGG        HWI-ST1391:140326:C441WACXX:7:1114:17283:59379 1:N:0: CAACAAT   TCTTGGGGGGCTCAGG        0-
MLT1H   MaLR family     Homo sapiens    188     CTTGAGGCCTAGGCT HWI-ST1391:140326:C441WACXX:7:1113:15313:71495 1:N:0: CAACTAA   CTTGAGGACTAGGCT 1       -
LTR29   LTR     Homo sapiens    43      ATGCTGAACTGAAGAAGCCTCAAGGT      HWI-ST1391:140326:C441WACXX:7:2310:17684:30008 1:N:0: CAACTAA   ATGCTGAACTGAAGAAGCCTCAAGGT0       +
MLT2F   LTR     Homo sapiens    582     TCTGATTGAATCCTGACTGATACA        HWI-ST1391:140326:C441WACXX:7:1213:21078:77562 1:N:0: CAACTAA   TCTGGTTGAATCCTGACTGATACA 1+
MER68A  LTR Retrotransposon     Eutheria        267     CCTGGGCACTGAGTC HWI-ST1391:140326:C441WACXX:7:2216:14800:78111 1:N:0: CAACTAA   CCGGGGCACTGAGTC 1       +
Kanga1  Mariner/Tc1     Eutheria        266     TTTCATTTCAATGCC HWI-ST1391:140326:C441WACXX:7:2304:3686:56287 1:N:0: CAACTAA    TATCATTTCAATGCC 1       -
HERVK3I endogenous retrovirus   Homo sapiens    332     GTTAATTCTCAGACA HWI-ST1391:140326:C441WACXX:7:2201:19850:22283 1:N:0: CAACTAA   GTTAATTCTCAGACC 1       +

```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the program

Explain how to run the automated tests for this system

### Files

Explain what these tests test and why

```
Give an example
``` 

## Authors

* **Stephanie Schustermann** - *Initial work* - [StephSchustermann](https://github.com/steohschustermann)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.


## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
