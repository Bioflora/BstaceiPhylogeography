# ENM, niche breadth, overlap, and Bioclimatic analysis

In this section, you can find the scripts and input data for the environmental niche analysis.

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium stacei* data included in the paper: 
>
> "Repeated migration, interbreeding, and bottlenecking  shaped the phylogeography and niche adaptation of the ancestral selfing Mediterranean grass *Brachypodium stacei*" by Miguel Campos Cáceres.
>
> And co-authored by Ernesto Pérez-Collazos, Antonio Díaz-Pérez, Diana López-Alvarez, Luis Mur, John Vogel and Pilar Catalán. 

## Table of Contents
* [Environmental Niche Modeling](#environmental_niche_modeling)
* [Bioclimatic Table](#bioclimatic_table)
* [PCA](#pca)

## Environmental_Niche_Modeling
For the environmental niche modeling, we need 3 files: 
- `Coordinates.csv` that contains the name of populations, and the coordinates of each one.
- `ENM_MC.R` the main script with all the commands used in R, includes the niche breadth and overlap.
- `functionsENM.r`contains the functions used for bootstrap analysis.

Additionally to, you may need to download the bioclimatic rasters, the current layers can be downloaded using the command line (included in the R script) but the past projections had to be downloaded directly from WorldClim:
https://www.worldclim.org/data/worldclim21.html

## Bioclimatic_Table
In order to generate a bioclimatic table with the values of all variables for our populations we need the following:
- `Coordinates.csv` that contains the coordinates and the name of our populations
- `Bioclimatic_Table.R` is the main script to generate the table
After obtaining the table, the script has the statistical comparison of the variables.

## PCA
To generate the bioclimatic PCA and detect potential climatic groupings we need 2 files:
- `Table_Sta.csv` the table generated previously with the values for all our populations.
- `PCA.R` is the main script to generate the PCA and save the plots.
