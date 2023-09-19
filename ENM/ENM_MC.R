#######################################################################################
#######################################################################################
################### ENM general Script Miguel Campos Enero 2022 #######################
#######################################################################################
#######################################################################################




#######################
## Install packages ##
######################
# install.packages("plotmo", dep = T)
# install.packages("dismo", dep = T)
# install.packages("raster", dep = T)
# install.packages("rgeos", dep=TRUE)
# install.packages("HH", dep=TRUE)
# install.packages("vegan", dep = T)
# install.packages("rgdal", dep=TRUE)
# install.packages("ecospat", dep=TRUE)
# install.packages("devtools", dep=TRUE)
# install.packages("maxnet", dep=TRUE)
# install_github("danlwarren/ENMTools")



#####################
## Load libraries ##
####################
library(plotmo)
library(dismo) 
library(raster) 
library(rgeos)
library(HH) 
library(vegan)
library(rgdal)
library(ecospat) 
library(devtools)
library(ENMTools)
library(maxnet)
library(dplyr)
library(dismo)
library(maptools)


source("functionsENM.R") #Para el análisis de bootstrap


#Create the folder for the results
setwd("/YOUR/PATH/ENM")
dir.create("Results")

# Load coordinates
data<-read.csv("Coordinates.csv",sep=";",dec=".")


#remove duplicated data
duplicated<-duplicated(data[ , c("Latitude", "Longitude")])
length(duplicated[duplicated==TRUE])
data<-data[!duplicated, ]

#Check for more duplicated data
duplicated<-duplicated(data[ , c("Latitude", "Longitude")])
length(duplicated[duplicated==TRUE])

# VGenerate a variable with only the coordinates
nrow(data)
presences <- dplyr::select(data, -Pop)

########################################### 
##      Download bioclimatic data        ##
## If you already have it, only load it  ##
###########################################

#bioclimatic <- stack(getData('worldclim', var='bio', res=2.5))

#Load the data and visualice them to check
bioclimatics <- stack(list.files(path="bioclimatics",pattern='*.bil', full.names=TRUE))

# Crop the map to the desired extension and check
e<-extent(-20,60,20,60)
bioclimatics<-crop(bioclimatics,e)
x11(width = 15, height = 10)
plot(bioclimatics[[1]])


# Generate a map with the presence points and save it in PDF
pdf("./Results/Presence.pdf", width=15, height=10, pointsize=20)
plot(bioclimatics[[1]], main="Pops of B. stacei in bio1")
points(presences)
dev.off()


########################
## Variable selection ##
########################

# Transform maps to a table and discarg null values
bioclimatics.table<-as.data.frame(bioclimatics)
bioclimatics.table<-na.omit(bioclimatics.table)

# Create the correlation matrix
bioclimatics.correlation<-cor(bioclimatics.table)
bioclimatics.correlation

# Create the distance matrix in absolute values
bioclimatics.dist<-as.dist(abs(bioclimatics.correlation))
bioclimatics.dist

# Build the cluster and save
bioclimatics.cluster<-hclust(1-bioclimatics.dist)
plot(bioclimatics.cluster)
pdf("./Results/Bioclimatic correlation.pdf", width=15, height=10, pointsize=20)
plot(bioclimatics.cluster)
dev.off()

# Select non-correlated variables
bioclimatics.selected<-c("bio7", "bio8", "bio6", "bio14", "bio19")
bioclimatics.table2<-bioclimatics.table[ , bioclimatics.selected]

# Check the variance inflation index for doing the final selection (should be less than 5)
bioclimatics.df <- na.omit(raster::as.data.frame(bioclimatics.table2))
vif(bioclimatics.df)   
bioclimatics<-brick(bioclimatics[[bioclimatics.selected]])



#####################
## Model Creation ## 
####################

names(presences)<-c("x","y")
presences.sp<-presences
coordinates(presences.sp)<-c("x","y")
backg<-randomPoints(bioclimatics[[1]],n=20000,p=presences.sp, excludep=T)

plot(bioclimatics[[1]])
points(backg)


## Create a new list object
sp <- list()

#Prepare presence data
sp$presence <- data.frame(
  presences,                                              
  raster::extract(
    x = bioclimatics,
    y = presences,             
    df = TRUE,
    cellnumbers = FALSE
  )
)
sp$presence$ID <- NULL
sp$presence$presence <- 1
names(sp$presence)[names(sp$presence) == "Longitude"] <- "x"
names(sp$presence)[names(sp$presence) == "Latitude"] <- "y"


# Add background data
sp$background <- data.frame(
  backg,
  raster::extract(
    x = bioclimatics,
    y = backg,
    df = TRUE,
    cellnumbers = FALSE
  )
)
sp$background$ID <- NULL
sp$background$presence <- 0

pb <- rbind(
  sp$presence, 
  sp$background
)

pb<-pb[c(3:8)] #### Select only variables and presence/absence columns  

write.csv(pb,"./Results/matrix_for_maxent.csv",row.names = F) 


############
## MAXENT ##
############

pb<-read.csv("./Results/matrix_for_maxent.csv")

temp.maxent <- maxnet::maxnet(
  p = pb$presence, 
  data = pb[,1:5]  ### select only variables
)

#Predictor coefficients
data.frame(
  importance = sort(
    abs(temp.maxent$betas)[1:9], 
    decreasing = TRUE
  )
)

######################
## Response curves ##
#####################
x11(width = 15, height = 10)
pdf("./Results/ResponseCurves.pdf", width=15, height=10, pointsize=20)
plot(temp.maxent, type="logistic")
dev.off()

################
## Prediction ##
################

temp.maxent.map <- raster::predict(
  object = bioclimatics, 
  model = temp.maxent, 
  type="logistic"
)

#Plot la rediction
x11(width = 15, height = 10)
pdf("./Results/Model.pdf", width=15, height=10, pointsize=20)
raster::plot(temp.maxent.map, col = viridis::turbo(100, direction = 1))
dev.off()

writeRaster(temp.maxent.map, filename='./Results/ModelMap', format='EHdr') #Save raster

#AUC and Threshold sensitivity
pres<-extract(temp.maxent.map,presences)
pres<-pres[!is.na(pres)]
aus<-extract(temp.maxent.map,backg)
evaluate=dismo::evaluate(p=pres,a=aus)
dismo::threshold(evaluate, sensitivity=0.9)

bootstrap.output <-bootstrapMaxent(data = pb,
                                   columna.presencia = "presence",
                                   columnas.variables = names(bioclimatics),
                                   iteraciones = 100, # at least 30 recommended
                                   porcentaje.evaluacion = 25, # between 20 and 40%
                                   variables=bioclimatics)
names(bootstrap.output)
bootstrap.output$auc.mean
bootstrap.output$auc.sd
plot(bootstrap.output$models.mean)
plot(bootstrap.output$models.sd)
writeRaster(bootstrap.output$models.mean, filename='./Results/ModelMean', format='EHdr')
writeRaster(bootstrap.output$models.sd, filename='./Results/ModelSD', format='EHdr')

#Threshold and save raster
Model <- raster("./Results/ModelMap")
Model <- bootstrap.output$models.mean
Threshold <- Model
Threshold[Threshold<.2180155] <- 0 #we've to set the number that we obtained in dismo::threshold(evaluate, sensitivity=0.9)
plot(Threshold, main="Model with Threshold 10%", col = viridis::turbo(100, direction = 1))
writeRaster(Threshold, filename='./Results/ModelThreshold', format='EHdr')


######################
## Past Projections ##
######################
# Load variables
list.variables<-names(bioclimatics)
variables.current<-bioclimatics[[list.variables]]


list.variables.midholocene <- list.files(path="mhbioclimatics",pattern='*.tif', full.names=TRUE)
list.variables.lgm <- list.files(path="lgmbioblimatics",pattern='*.tif', full.names=TRUE)
list.variables.lig <- list.files(path="ligbioclimatics",pattern='*.bil', full.names=TRUE)


variables.midholocene <- stack(list.variables.midholocene)
variables.lgm <- stack(list.variables.lgm)
variables.lig <- stack(list.variables.lig)

#Crop and select variables
e<-extent(-20,60,20,60)
variables.midholocene<-crop(variables.midholocene,e)
variables.lgm<-crop(variables.lgm,e)
variables.lig<-crop(variables.lig,e)

variables.midholocene<-resample(x=variables.midholocene, y=variables.current, method="bilinear")
variables.lgm<-resample(x=variables.lgm, y=variables.current, method="bilinear")
variables.lig<-resample(x=variables.lig, y=variables.current, method="bilinear")

# Select only same variables than current
variables.midholocene<-variables.midholocene[[list.variables]]
variables.lgm<-variables.lgm[[list.variables]]
variables.lig<-variables.lig[[list.variables]]

# Calibrate the model
# We've to use the current model and then project
mh.maxent.map<-predict(variables.midholocene, temp.maxent, type="logistic")
lgm.maxent.map<-predict(variables.lgm, temp.maxent, type="logistic")
lig.maxent.map<-predict(variables.lig, temp.maxent, type="logistic")

#Plot, AUC and threshold 
bootstrap.output.mh <-bootstrapMaxent(data = pb, 
                                      columna.presencia = "presence",
                                      columnas.variables = names(variables.midholocene), 
                                      iteraciones = 100, 
                                      porcentaje.evaluacion = 25,
                                      variables=variables.midholocene)
names(bootstrap.output)
bootstrap.output$auc.mean
bootstrap.output$auc.sd

pres<-extract(bootstrap.output.mh$models.mean,presences)
pres<-pres[!is.na(pres)]
aus<-extract(bootstrap.output.mh$models.mean,backg)
evaluate1=dismo::evaluate(p=pres,a=aus)
dismo::threshold(evaluate1, sensitivity=0.9)
Threshold1 <- bootstrap.output.mh$models.mean
Threshold1[Threshold1<.1190381] <- 0
pdf("./Results/Model MH Threshold.pdf", width=15, height=10, pointsize=20)
plot(Threshold1, main="MH Model with Threshold 10%", col = viridis::turbo(100, direction = 1))
dev.off()
writeRaster(Threshold1, filename='./Results/Model MH Threshold', format='EHdr')



bootstrap.output.lgm <-bootstrapMaxent(data = pb,
                                       columna.presencia = "presence",
                                       columnas.variables = names(variables.lgm), 
                                       iteraciones = 100,
                                       porcentaje.evaluacion = 25,
                                       variables=variables.lgm)
names(bootstrap.output)
bootstrap.output$auc.mean
bootstrap.output$auc.sd

pres<-extract(bootstrap.output.lgm$models.mean,presences)
pres<-pres[!is.na(pres)]
aus<-extract(bootstrap.output.lgm$models.mean,backg)
evaluate2=dismo::evaluate(p=pres,a=aus)
dismo::threshold(evaluate2, sensitivity=0.9)
Threshold2 <- bootstrap.output.lgm$models.mean
Threshold2[Threshold2<.2495656] <- 0
pdf("./Results/Model LGM Threshold.pdf", width=15, height=10, pointsize=20)
plot(Threshold2, main="LGM Model with Threshold 10%", col = viridis::turbo(100, direction = 1))
dev.off()
writeRaster(Threshold2, filename='./Results/Model LGM Threshold', format='EHdr')




bootstrap.output.lig <-bootstrapMaxent(data = pb, 
                                       columna.presencia = "presence",
                                       columnas.variables = names(variables.lig), 
                                       iteraciones = 100, 
                                       porcentaje.evaluacion = 25,
                                       variables=variables.lig)

names(bootstrap.output)
bootstrap.output$auc.mean
bootstrap.output$auc.sd

pres<-extract(bootstrap.output.lig$models.mean,presences)
pres<-pres[!is.na(pres)]
aus<-extract(bootstrap.output.lig$models.mean,backg)
evaluate3=dismo::evaluate(p=pres,a=aus)
dismo::threshold(evaluate3, sensitivity=0.9)
Threshold3<- bootstrap.output.lig$models.mean
Threshold3[Threshold3<.06307618] <- 0
pdf("./Results/Model LIG Threshold.pdf", width=15, height=10, pointsize=20)
plot(Threshold3, main="LIG Model with Threshold 10%", col = viridis::turbo(100, direction = 1))
dev.off()
writeRaster(Threshold3, filename='./Results/Model LIG Threshold', format='EHdr', overwrite=TRUE)

# Niche Breadth
Stacei <- raster('./Results/ModelThreshold')
StaceiMH <- raster('./Results/Model MH Threshold')
StaceiLGM <- raster('./Results/Model LGM Threshold')
StaceiLIG <- raster('./Results/Model LIG Threshold')


nicheBreadth(model.map=Stacei)
nicheBreadth(model.map=StaceiMH)
nicheBreadth(model.map=StaceiLGM)
nicheBreadth(model.map=StaceiLIG)

# Niche Overlap Present and Past
dismo::nicheOverlap(Stacei, StaceiMH, stat='D')
dismo::nicheOverlap(Stacei, StaceiMH, stat='I')

dismo::nicheOverlap(StaceiMH, StaceiLGM, stat='D')
dismo::nicheOverlap(StaceiMH, StaceiLGM, stat='I')

dismo::nicheOverlap(StaceiLGM, StaceiLIG, stat='D')
dismo::nicheOverlap(StaceiLGM, StaceiLIG, stat='I')

