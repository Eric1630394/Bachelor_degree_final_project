# Content stored in this repository
Hey!! 👋

Let me introduce the content inside this repository. I am Eric Garcia, an undergraduate student of the Genetic's degree at the Universitat Autònoma de Barcelona (UAB) 🏛️. This year, as I face the completion of my studies, I am elaborating my bachelor's degree final project. 

My reserach has been focused in studying association between alleles of Single Nucleotide Polymorphisms (SNPs) and Structurals Variants (SVs), specifically inversions. With this objective in mind, I created a collection of shell 🐚 and R 📊 scripts to run my analyses and automatize the treatment of big-data datasets. 

This repository has been divided in different folders 📁 (introduced here) depending on the material used in every phase. In each folder one can find the code itself and other files passed to the scripts or generated from them. 

## 📁 Input files

This directory contains all the files that are passed to the different scripts (stored in the other two folders), necessary to generate the results. A general description and the usage of every file in this folder is attached here, ordered according to the different phases of the analysis:
| Input file | Description |
|------------|-------------|
|invappend2_v2.in|A list of all the inversions (n=254, after manual filtering) analysed in this study. This input file is passed to 🐚`Monomorphic_class.sh`.  
|Inversions_coordinates_porubsky_v3.csv|For every inversion, this file contains information on the chromosome (2nd column), genomic coordinates for the breakpoints (BP1-1, BP1-2, BP2-1, BP2-2, columns 3-6) and the flanking class (7th column). This file is passed to 📊 `Inversions_filtering.r` and 📊 `United3_plot_statitsics.r`. 
|Genotyping_info.txt|Output generated from 🐚 `Monomorphic_class.sh`. It contains genotyping information for every inversion depending on the genotypes of the 43 samples analysed in this study. Inversions are classified as: only 1 sample genotyped, chrY inversions, monomorphic and polymorphic. This file is the passed to 📊 `Inversions_filtering.r`. 
|United3.ld|An unified file where, for every inversion, it is stored the chromosome (2nd column), the genomic position for the SNP (3rd column), the biallelic variants for that site (4th column) and the LD value (5th column). This file is the passed to 📊 `United3_plot_statistics.r`.

## 📁 Inversions_filtering

## 📁 LD_measures


