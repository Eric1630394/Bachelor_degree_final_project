# Content stored in this repository
Hey!! 👋

Let me introduce the content inside this repository. I am Eric Garcia, an undergraduate student of the Genetic's degree at the Universitat Autònoma de Barcelona (UAB) 🏛️. This year, as I face the completion of my studies, I am elaborating my bachelor's degree final project. 

My research consisted in the evaluation of association between alleles of Single Nucleotide Polymorphisms (SNPs) and both orientations of polymorphic inversions (standard and inverted) found in the human genome. This was done by using the correlation coefficient ($r^2) as a measure of linkage disequilibrium (LD). So, with this objective in mind, I created a collection of shell 🐚 and R 📊 scripts to run my analyses and automatize the treatment of big-data datasets. 

These scripts can be found in this repository by their classification on different folders 📁 depending on the material used in every phase. In each folder one can find the code itself and other files passed to the scripts or generated from them. 

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
|`Inversions_filtering.r`|makes the proper visualization of the results generated, shown in `Barplot_filtering.png` and `Barplot_filtering_no_low.png`|`Genotyping_info.txt` and `Inversions_coordinates_porubsky_v3.csv`|

Results coming from this filtering can be visualized in `Barplot_filtering.png` and `Barplot_filtering_no_low.png` (depending on if low confidence genotypes have been considered or not), where every category of inversions has been dividied in a different colour. The y-axis represents the proportion for the inversions in every category and the x-axis refers to every category depending on flanking information. 

## 📁 LD_measures

This directory stores information on the results generated after LD has already been calculated (information contained in `United3.ld`). A single script is used to make the proper visualization of the results and perform the appropriate statistic test:

|Scripts|Description|Necessary input files|
|-------|-----------|---------------------|
|`United3_plot_statistics.r`|reviews information for the LD values between polymorphic inversions (n=203) and SNPs (`United3.ld`)| `United3.ld` and `Inversions_coordinates_porubsky_v3.csv`|

As for the results, four images are attached here:
  - `United3_ld_plot.png`: graph representing all the LD values (y-axis) for every SNP-inversion pair. A total of 180 inversions is shown in the x-axis (those with at least 1 SNP in LD > 0.35). These results have been divided by categories on flanking information.
  - `Histogram.png`: distribution of the different LD values, divided by categories on flanking information. 
  - `r2_values_perfect_ld.png`: inversions with SNPs in perefct LD. The y -axis represents the absolute count of SNPs for every inversion. These results have been divided by categories on flanking information.  
  - `Barplot_comparison.png`: comparison between the original proportion of inversions for each flanking category and the proportion among inversions with at least 1 tagSNP. Asterisks indicate significance for the associated adjpval coming from the Fisher's Exact Test. 


