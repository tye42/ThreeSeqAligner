# ThreeSeqAligner

[Project site](https://tye42.github.io/2016/05/31/Reduce-Space-Complexity-in-Three-Sequence-Alignment.html)

Three DNA sequence alignment using divide-and-conquer and dynamic programming.

Usage:

```
./threeSeqAligner seq1.txt seq2.txt seq3.txt > output.txt
```

The output format is BLAST-like

```
# Sequences: seq1.txt: 556  seq2.txt: 569  seq3.txt: 627
# Running time: 9.55284
#===================================================================================================================
#
# Length: 631
# Score: 4752
# Identity: 427/631 (67.67%)
#
#===================================================================================================================
#
# 1     -----ACA-T-------T--C--------TC-CTTCTG------ATAGACTC--AG-GAAGCAATCATGGTGCTCTCTGCAGATGACAAAACCAACATCA    67
# 1     ----GACACT-------T--C--------TG-ATTCTG------ACAGACTC--AG-GAAGAAACCATGGTGCTCTCTGGGGAAGACAAAAGCAACATCA    69
# 1     CATAAACCCTGGCGCGCTCGCGGCCCGGCACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCACCATGGTGCTGTCTCCTGCCGACAAGACCAACGTCA   100
#            **  *       *  *            *****      * ******  ** ***   * ********* ***   *  ***** * **** ***
#
# 68    AGAACTG-CTGGGGGAAGATTGGTGGCCATGGTGGTGAATATGGCGAGGAGGCCCTACAGAGGATGTTCGCTG-CCTTCCCCACCACCAAGACCTACTTC   165
# 70    AGG-CTGCCTGGGGGAAGATTGGTGGCCATGGTGCTGAATATGGAGCTGAAGCCCTGGAAAGGATGTTTGCTA-GCTTCCCCACCACCAAGACCTACTTC   167
# 101   AGG-CCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTC-CTGTCCTTCCCCACCACCAAGACCTACTTC   198
#       **  * * ****** *** * ** *  ** * **  ** ***** *  ** *****  * ********  **   *************************
```
