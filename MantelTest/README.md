# Isolation by Distance (IBD) and Isolation by Environment (IBE)

In this section, you can find the scripts and input data for the mantel test performed to test IBD and IBE.

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium stacei* data included in the paper:
>
> "Repeated migration, interbreeding, and bottlenecking  shaped the phylogeography and niche adaptation of the ancestral selfing Mediterranean grass *Brachypodium stacei*" by Miguel Campos Cáceres.
>
> And co-authored by Ernesto Pérez-Collazos, Antonio Díaz-Pérez, Diana López-Alvarez, Luis Mur, John Vogel and Pilar Catalán. 

## Isolation by Distance (IBD) and Isolation by Environment (IBE)
For the mantel test performed, we need 4 files: 
- `Coordinates.txt` that contains the name of populations, and the coordinates of each one.
- `MantelTest.R` is the main script with all the commands used in R.
- `Fst.txt` genetic differentiation values obtained from Arlequin.
- `Bioclimaticas.csv` mean values of bioclimatic variables for each population, using only non-correlated variables.
