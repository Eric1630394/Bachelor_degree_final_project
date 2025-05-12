# Content stored in this repository
Hey!! 👋

Let me introduce the content inside this repository. I am Eric Garcia, an undergraduate student of the Genetic's degree at the Universitat Autònoma de Barcelona (UAB) 🏛️. This year, as I face the completion of my studies, I am elaborating my bachelor's degree final project. 

My reserach has been focused in studying association between alleles of Single Nucleotide Polymorphisms (SNPs) and Structurals Variants (SVs), specifically inversions. With this objective in mind, I created a collection of shell 🐚 and R 📊 scripts to run my analyses and automatize the treatment of big-data datasets. 

This repository has been divided in different folders 📁 (introduced here) depending on the material used in every phase. In each folder one can find the code itself and other files passed to the scripts or generated from them. 

## 📁 Input files

This directory contains all the files that are passed to the different scripts (stored in the other two folders), necessary to generate the results. A general description and the usage of every file in this folder is attached here, ordered according to the different phases of the analysis:
| Input file 📖 | Description ✍️|
|------------|-------------|
|invappend2_v2.in|A list of all the inversions (n=254, after manual filtering) analysed in this study 
|Porubsky_invs_callset_GT_v4.csv|Genotypes for all the samples in this study
|Inversions_coordinates_porubsky_v3.csv|For every inversion, this file contains information on the chromosome (2nd column), genomic coordinates for the breakpoints (BP1-1, BP1-2, BP2-1, BP2-2, columns 3-6) and the flanking class (7th column)
|Genotyping_info.txt|Output generated from 🐚 `Monomorphic_class.sh`. It contains genotyping information for every inversion depending on the genotypes of the 43 samples analysed in this study. Inversions are classified as: only 1 sample genotyped, chrY inversions, monomorphic and polymorphic
|United3.ld|An unified file where, for every inversion, it is stored the chromosome (2nd column), the genomic position for the SNP (3rd column), the biallelic variants for that site (4th column) and the LD value (5th column)

## 📁 Inversions_filtering

This directory stores those scripts and results generated for the inversions filtering depenidng on the genotypes information for all the samples studied in this research. 

|Scripts|Description|Necessary input files|
|-------|-----------|---------------------|
|`Monomorphic_class.sh`|reviews information for the genotypes of the samples and generates an output, `Genotyping_info.txt`, with information for genotyping category.|`invappend2_v2.in`, `Porubsky_invs_callset_GT_v4.csv` and `Inversions_coordinates_porubsky_v3.csv`|
|`Inversions_filtering.r`|makes the proper visualization of the results generated, shown in `Barplot_filtering.png`.|`Genotyping_info.txt` and `Inversions_coordinates_porubsky_v3.csv`|

Results coming from this filtering can be visualized in `Barplot_filtering.png`, where every category of inversions has been dividied in a different colour. The y-axis represents the proportion for the inversions in every category and the x-axis refers to every category depending on flanking information. 

## 📁 LD_measures


