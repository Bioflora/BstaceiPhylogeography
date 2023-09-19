###########################################
#### DAPC Analysis with plastomic data ####
####### Miguel Campos Caceres 2023 #######
###########################################

library(adegenet)
library(dartR)

Plastomes <- gl.read.fasta("Plastomes.fasta")
Pops <- read.table("Pops_plastomes.txt", header = TRUE, sep = "\t")
Pop <- Pops$Pop # Replace 'Population' with the column name in your file containing population labels
pop(Plastomes) <- Pop

## PCoA
pcoa <- gl.pcoa(Plastomes)

## Group clustering usign a PCoA, we've to select the number of clusters that are above the AIC line
grp<-find.clusters(Plastomes, max.n.clust=10, method='ward', stat='AIC', glPca = pcoa)


# If we've previously stablished populations we can check if they follow that grouping for example geographic
table<-table(pop(Plastomes),grp$grp)
table.value(table(pop(Plastomes),grp$grp),col.lab=paste("inferred",1:5), row.lab=paste("original",1:10))

# Now the DAPC analysis
dapc1<-dapc(Plastomes,grp$grp)
dapc1

# Nex step is to check the optimal number of PC's for the DAPC, if its less we can repeat the previous step
temp<-optim.a.score(dapc1)

# Plot the DAPC
scatter(dapc1)

myCol<-c("blue","purple","darkgreen","orange", "darkred")
scatter(dapc1,scree.da=FALSE,bg="white",pch=20,cell=0,cstar=0,col=myCol,solid=.4, cex=3,clab=0,leg=TRUE,txt.leg=paste("Cluster",1:4))

#Plot with a minimum spanning network
scatter(dapc1,ratio.pca=0.3,bg="white",pch=20,cell=0, cstar=0,col=myCol,solid=.4,cex=4,clab=0, mstree=TRUE,scree.da=FALSE,posi.pca="bottomright", leg=TRUE,txt.leg=paste("Cluster",1:5)) 
par(xpd=TRUE) 
points(dapc1$grp.coord[,1],dapc1$grp.coord[,2],pch=4, cex=2,lwd=5,col="black") 
points(dapc1$grp.coord[,1],dapc1$grp.coord[,2],pch=4, cex=2,lwd=2,col=myCol)
myInset<-function(){ 
  temp<-dapc1$pca.eig 
  temp<-100*cumsum(temp)/sum(temp) 
  plot(temp,col=rep(c("black","lightgrey"), c(dapc1$n.pca,1000)),ylim=c(0,100), xlab="PCAaxis",ylab="Cumulatedvariance(%)", cex=1,pch=20,type="h",lwd=2) 
} 
add.scatter(myInset(),posi="bottomleft", inset=c(-0.03,-0.01),ratio=.28, bg=transp("white"))

# Structure-like plot
compoplot(dapc1,posi="bottomright", txt.leg=paste("Cluster",1:5),lab="", ncol=1,xlab="individuals",col=funky(5))

# To see member assignation for other analysis 
round(dapc1$posterior,1)



