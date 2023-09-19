#######################################
##### Mantel test for IBD and IBE #####
##### Miguel Campos Cáceres 2023 #####
######################################

library(vegan)
library(geosphere)

# Define the coordinates for each population
coordinates <- read.table("Coordinates.txt", sep="\t", header=TRUE)

# Calculate geographic distances
distances <- distm(coordinates[, c("Longitude", "Latitude")], fun = distHaversine)

# Create a genetic distance matrix
genetic_distances <- read.table("Fst.txt", sep="\t", header=TRUE, row.names=1)

# Perform Mantel test with linearized FST distances and geographic distances
mantel_result_GD <- mantel(genetic_distances, distances, method="pearson", permutations=10000)
print(mantel_result_GD)

# Calculate bioclimatic distances using the mean of populations and non correlated variables
data <- read.csv("Bioclimatic.csv", row.names = 1, sep="\t")
bioclim_vars <- data[, -1]
bioclim_dist <- vegdist(bioclim_vars)

# Perform Mantel test with linearized FST distances and climatic distances
mantel_result_BC <- mantel(genetic_distances, bioclim_dist, method="pearson", permutations=10000)
print(mantel_result_BC)

