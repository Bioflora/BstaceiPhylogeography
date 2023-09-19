#### Script to generate a bioclimatic table with the values per population
#### Modified from https://damariszurell.github.io/EEC-MGC/a1_SpatialData.html#34_Extract_raster_data


##Install packages and load library

#install.packages('geodata')
library(geodata)

# Download bioclimatic variables
(clim <- geodata::worldclim_global(var = 'bio', res = 2.5, download = F, path = 'data'))
plot(clim[[1]])

# Load our population data and select only coordinates
data<-read.csv("Coordinates.csv",sep=";",dec=".")
presences <- dplyr::select(data, -Pop)
head(presences)

# Extract the values of our coordinates from bioclimatic data
table <- terra::extract(clim, presences)
na.exclude(table)

tablelonglat <- cbind(presences, terra::extract(clim, presences))
head(tablaconlonglat)

# Export table
write.csv(table,"./Table_Sta.csv") 


