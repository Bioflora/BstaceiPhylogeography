# Admixture analysis and LD filtering

In this section you can find the scripts used and input file for admixture analysis.

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium stacei* data included in the paper:
>
> "Repeated migration, interbreeding, and bottlenecking  shaped the phylogeography and niche adaptation of the ancestral selfing Mediterranean grass *Brachypodium stacei*" by Miguel Campos Cáceres.
>
> And co-authored by Ernesto Pérez-Collazos, Antonio Díaz-Pérez, Diana López-Alvarez, Luis Mur, John Vogel and Pilar Catalán. 

## Table of Contents
* [LD Filtering](#ld_filtering)
* [Input generation](#input_generation)
* [Run Admixture](#run_admixture)

## LD Filtering
In order to avoid the linkage disequilibrium (LD) we filtered the data using the next commands in R:
```
## Package installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")

# Load libraries
library(SNPRelate)
library(MASS)

# Set directory and file
setwd("~/LD")
vcf.fn <- "Stacei.vcf" ## archivo vcf

# Filter the SNPs based on the LD statistics
snpgdsVCF2GDS(vcf.fn, "Stacei.gds", method="biallelic.only")
genofile <- snpgdsOpen("Stacei.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, remove.monosnp=TRUE, slide.max.bp=10000)
snpset.id <- unlist(snpset)
snpgdsGDS2BED(genofile,"Stacei_Filtered", snp.id=snpset.id, verbose=TRUE)
# With plink outside R we can transform it in a vcf again:
plink2 --bfile Stacei_Filtered --recode vcf --out Stacei_Filtered
```

## Input generation
Admixture requires a bed file as input. We started from a Variant Call File (.vcf) and we follow the next steps: 
- `Stacei_Filtered.vcf` file that contains the genomic data of B. stacei, filtered from LD.

Once we have the input data, we need to transform it to a bed file in order to use it in Admixture:
We used conda to install `admixture` and `plink`.
```
conda activate
conda install -c bioconda plink
conda install -c bioconda admixture
```
## Run Admixture 
Now, we used a loop to perform the Admixture analysis using K=1 to K=10. First we create the folder Run1 because we have to repeat this loop another 10 times, then we run the loop with a cv value of 10 and use 20 threads (-j20).
```
mkdir Run1
cd Run1
for i in {1..10}
do
 admixture -j20 --cv=10 Stacei_Filtered.bed $i > log${i}.out
done
```

Then, the last step is to get the cross-validation errors for each run, for that we use this command:
```
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > Stacei_Filtered.cv.error
```
