# Admixture analysis

In this section you can find the scripts used and input file for admixture analysis.

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium stacei* data included in the paper:
>
> "Repeated migration, interbreeding, and bottlenecking  shaped the phylogeography and niche adaptation of the ancestral selfing Mediterranean grass *Brachypodium stacei*" by Miguel Campos Cáceres.
>
> And co-authored by Ernesto Pérez-Collazos, Antonio Díaz-Pérez, Diana López-Alvarez, Luis Mur, John Vogel and Pilar Catalán. 

## Table of Contents
* [Input generation](#input_generation)
* [Run Admixture](#run_admixture)

## Input generation
Admixture requires a bed file as input. We started from a Variant Call File (.vcf) and we follow the next steps: 
- `Stacei.vcf` file that contains the genomic data of B. stacei, obtained from Ipyrad.

Once we have the input data, we need to transform it to a bed file in order to use it in Admixture:
We used conda to install `admixture` and `plink`.
```
conda activate
conda install -c bioconda plink
conda install -c bioconda admixture
```
Once we've installed the softwares, we need to make the bed file:
```
plink --vcf Stacei.vcf --make-bed --out Stacei --allow-extra-chr
```
The last command is to allow plink to use different chromosome names (if we dont use this only accepts human like names)

## Run Admixture 
The software does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
```
awk '{$1="0";print $0}' Stacei.bim > Stacei.bim.tmp
mv Stacei.bim.tmp Stacei.bim
```

Now, we used a loop to perform the Admixture analysis using K=1 to K=10. First we create the folder Run1 because we have to repeat this loop another 10 times, then we run the loop with a cv value of 10 and use 20 threads (-j20).
```
mkdir Run1
cd Run1
for i in {1..10}
do
 admixture -j20 --cv=10 Stacei.bed $i > log${i}.out
done
```

Then, the last step is to get the cross-validation errors for each run, for that we use this command:
```
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > Stacei.cv.error
```
