# Potpourri
#### An Epistasis Test Prioritization Algorithm via Diverse SNP Selection

## Getting Started:
Potpourri provides a MATLAB interface for ease of use. These instructions will guide you to build and run Potpourri on MATLAB.

## Requirements:
Building Potpourri requires the Boost C++ library.

## Installation:
In order to build Potpourri for MATLAB, just type make on the terminal!:
```
make
```
or directly run the MATLAB script for building mex files:
```
build_mex.m
```

## Input Format:
```
@ Feature Matrix: 
This should consist of a grid {0, 1, 2} characters, representing homozygous major, heterozygous and homozygous
minor genotypes respectively for all samples. Each row corresponds to a sample.
@ Labels:
This should consist of {0, 1} binary labels representing control and case respectively.
@ SNP Information:
This should consist of three columns: unique SNP_id, chromosome and position
@ Regulatory/Coding Information:
This should consist of {0, 1} binary labels representing whether the corresponding SNP (with columns of feature matrix) is in the regulatory region (1) or not (0).
@ Network Matrix:
An adjacency matrix for SNP-SNP interaction. 
```
## Parameters:
```
@ Maximum marginal significance:
Takes an integer value from 1-6, representing the maximum marginal significance of loci for consideration in pairwise testing as a -log10(p-value).
@ outputFileName:
Prefix for output files.
@ Omega: 
A float parameter of Potpourri to reward regulatory region.
@ b:
Number of neighbors that should be included in the epistasis test for each selected SNP.
@ k:
Number of features to be selected 
```
## Examples:
How to run SPADIS on MATLAB. 
Simply run the demo file:
```
demo_potpourri.m
```

The example data is adapted from Atwell et. al. (2010). The genotype and phenotype data of Arabidopsis Thaliana (AT) obtained from [Atwell et. al. (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20336072) and adapted to the algorithm accordingly. For descriptions and format of the data, check the [readme file for data](data/readme_data.txt).
## Output
After running Potpourri, files above will be either created, or appended to, in the output directory.
###.Summary
The summary file contains several statistics such as statistical tests performed, number of reciprocal pairs found, 
Each row corresponds to a separate run.
###.Cutoff Pairs
Contains all detected pairs above the dynamic significance threshold at the conclusion of a run.
###.Reciprocal Pairs
Contains a subset of the pairs in the cutoff pairs set such that both loci in a pairing had the other locus selected
as its most significant interaction.
###.Reciprocal Pairs Formatted
Each row represents a single reciprocal locus pairing with chi-squared significance, locus 1 and 2, chromosome 1 and 2, base pair 1 and 2. 


## License
This project is licensed under GNU GPL v3 - see the [LICENSE.md](LICENSE.md) file for details.

## References
Yilmaz, Serhan, Tastan, Oznur & Cicek, A. Ercument (2018). [SPADIS: An Algorithm for Selecting Predictive and Diverse SNPs in GWAS](https://www.biorxiv.org/content/early/2018/01/30/256677). bioRxiv

Atwell, S. et al. (2010) [Genome-wide association study of 107 phenotypes
in Arabidopsis thaliana inbred lines](https://www.ncbi.nlm.nih.gov/pubmed/20336072). Nature, 465(7298), 627–631.

Wu, M. C. et al. (2011) [Rare-variant association testing for sequencing
data with the sequence kernel association test](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135811/). The American Journal of Human Genetics, 89(1), 82–93.

Cowman T, Koyutürk M. (2017) [Prioritizing tests of epistasis through hierarchical representation of genomic redundancies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5737499/). Nucleic acids research, 45(14), e131.

