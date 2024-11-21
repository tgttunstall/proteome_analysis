#!/usr/bin/perl

perl -ne 'print if $. == 1 || rand() < .01' Specie_protein_cluster.tsv > sample_Specie_protein_cluster.tsv; wc -l sample_Specie_protein_cluster.tsv

#perl -ne 'print if $. == 1 || rand() < .01' FILE1 > FILE2; wc -l FILE2.csv
#perl -ne 'print if $. == 1 || rand() < .01' af2-models-split-ifresid_.tsv > sample_af2.tsv; wc -l sample_af2.tsv

