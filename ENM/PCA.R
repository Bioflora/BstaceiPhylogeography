################################
##      Bioclimatic PCA       ##
## Miguel Campos Cáceres 2023 ##
################################


########## Packages
library(dplyr)
library(ade4)
library(factoextra)
library(magrittr)
library(HH) 
source("vif_fun.R")

########## Load the data
data<-read.csv ("Table_Sta.csv" , header=T , sep=";" , dec=".")
names(data)
summary(data)
group<-data[,1] # Grouping factor
group=as.factor(group)
variables<-data[,5:23] # Variables the PCA is based on


########### Select non-correlated variables
variables.tabla<-as.data.frame(variables)
variables.tabla<-na.omit(variables.tabla)
variables.correlacion<-cor(variables.tabla)
variables.dist<-as.dist(abs(variables.correlacion))
variables.cluster<-hclust(1-variables.dist)
plot(variables.cluster)

pdf("./Bioclimatic correlation.pdf", width=15, height=10, pointsize=20)
plot(variables.cluster)
dev.off()

variables <- dplyr::select(variables, bio13, bio14, bio3, bio9, bio8) 

########## PCA
pca<-dudi.pca(variables) #Select 3 axes
fviz_eig(pca) #Bars with the explained variance

s.label(pca$li,clab=0) 
s.class(pca$li, fac = group, col = c('darkgreen', 'darkorange', 'purple'), 
        add.plot = TRUE, cstar = 0, cellipse=0, cpoint=2)          
text(pca$li,labels=data$Pop,adj=c(.1,.8), cex=0.8)


### Contribution values
inertia.dudi(pca) #Axis contributions
inertia.dudi(pca,col.inertia=T) #Variables contributions
inertia.dudi(pca,row.inertia=T) #Points contributions


# Results for individuals
pca.ind <- get_pca_ind(pca)
pca.ind$coord          # Coordinates


#Saving Plots
pdf("./Percentage of explained variances.pdf", width=15, height=10, pointsize=20)
fviz_eig(pca)
dev.off()

pdf("./Variables.pdf", width=15, height=10, pointsize=20)
s.corcircle(pca$co)
dev.off()

pdf("./PCA K3.pdf", width=15, height=15, pointsize=20)
s.label(pca$li,clab=0)
s.class(pca$li, fac = group, col = c('darkgreen', 'darkorange', 'purple'), 
        add.plot = TRUE, cstar = 0, cellipse=0, cpoint=2)
text(pca$li,labels=data$Pop,adj=c(.1,.8), cex=0.8)
dev.off()
