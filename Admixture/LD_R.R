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
