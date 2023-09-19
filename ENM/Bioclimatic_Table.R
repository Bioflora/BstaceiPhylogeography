###########################################################
#### Statistical differences in bioclimatic variables ####
############## Miguel Campos June '23 ###################
##########################################################


##Install packages and load library
library(geodata)
library(nortest)
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(dplyr)
library(car)
library(ggplot2)
library(patchwork)
library(MASS)
library(multcompView)

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
write.csv(tablaconlonglat,"./Table_Sta.csv") 

# Once we have the data we can start the test:
datos <- tablaconlonglat
################################## 
## Normality test of Lilliefors ##
##################################
#Less than 0.05 no normal
lillie.test(datos$Bio1) 
lillie.test(datos$Bio2) 
lillie.test(datos$Bio3) 
lillie.test(datos$Bio4)
lillie.test(datos$Bio5) 
lillie.test(datos$Bio6)
lillie.test(datos$Bio7)
lillie.test(datos$Bio8) 
lillie.test(datos$Bio9) 
lillie.test(datos$Bio10)
lillie.test(datos$Bio11)
lillie.test(datos$Bio12)
lillie.test(datos$Bio13)
lillie.test(datos$Bio14)
lillie.test(datos$Bio15)
lillie.test(datos$Bio16)
lillie.test(datos$Bio17)
lillie.test(datos$Bio18)
lillie.test(datos$Bio19)


######################################
## Levene's test of homocedasticity ##
######################################
#Superior to 0.05 is homocedastic
datos$Group <- as.factor(datos$Group)
leveneTest(datos$Bio1, group=datos$Group, na.action = na.omit) 
leveneTest(datos$Bio2, group=datos$Group, na.action = na.omit) 
leveneTest(datos$Bio3, group=datos$Group, na.action = na.omit) 
leveneTest(datos$Bio5, group=datos$Group, na.action = na.omit) 
leveneTest(datos$Bio8, group=datos$Group, na.action = na.omit) 
leveneTest(datos$Bio9, group=datos$Group, na.action = na.omit) 
leveneTest(datos$Bio10, group=datos$Group, na.action = na.omit)
leveneTest(datos$Bio12, group=datos$Group, na.action = na.omit)
leveneTest(datos$Bio13, group=datos$Group, na.action = na.omit)
leveneTest(datos$Bio15, group=datos$Group, na.action = na.omit)
leveneTest(datos$Bio16, group=datos$Group, na.action = na.omit)
leveneTest(datos$Bio18, group=datos$Group, na.action = na.omit)
leveneTest(datos$Bio19, group=datos$Group, na.action = na.omit)

################
## Test Anova ##
################
#Superior to 0.05 no differences
Bio1 = aov(lm(datos$Bio1 ~ datos$Group))
summary(Bio1) 

Bio2 = aov(lm(datos$Bio2 ~ datos$Group))
summary(Bio2) 

Bio3 = aov(lm(datos$Bio3 ~ datos$Group))
summary(Bio3)

Bio5 = aov(lm(datos$Bio5 ~ datos$Group))
summary(Bio5) 

Bio8 = aov(lm(datos$Bio8 ~ datos$Group))
summary(Bio8) 

Bio9 = aov(lm(datos$Bio9 ~ datos$Group))
summary(Bio9) 

Bio10 = aov(lm(datos$Bio10 ~ datos$Group))
summary(Bio10)

Bio12 = aov(lm(datos$Bio12 ~ datos$Group))
summary(Bio12)

Bio13 = aov(lm(datos$Bio13 ~ datos$Group))
summary(Bio13) 

Bio15 = aov(lm(datos$Bio15 ~ datos$Group))
summary(Bio15)

Bio16 = aov(lm(datos$Bio16 ~ datos$Group))
summary(Bio16) 

Bio19 = aov(lm(datos$Bio19 ~ datos$Group))
summary(Bio19) 

##########################
## Test Kruskall-Wallis ##
##########################
#More than a 0.05 no differences

kruskal.test(Bio4 ~ Group, data = datos) 
kruskal.test(Bio6 ~ Group, data = datos) 
kruskal.test(Bio7 ~ Group, data = datos) 
kruskal.test(Bio11 ~ Group, data = datos)
kruskal.test(Bio14 ~ Group, data = datos)
kruskal.test(Bio17 ~ Group, data = datos)
kruskal.test(Bio18 ~ Group, data = datos)

#########################
## Post-Hoc Bonferroni ##
#########################
#More than a 0.05 are equal
pairwise.t.test(datos$Bio1, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio3, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio5, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio8, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio9, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio10, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio12, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio13, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio15, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio16, datos$Group, p.adjust.method = "bonferroni")
pairwise.t.test(datos$Bio19, datos$Group, p.adjust.method = "bonferroni")

pairwise.wilcox.test(x = datos$Bio4, g = datos$Group, p.adjust.method = "bonf" )
pairwise.wilcox.test(x = datos$Bio6, g = datos$Group, p.adjust.method = "bonf" )
pairwise.wilcox.test(x = datos$Bio7, g = datos$Group, p.adjust.method = "bonf" )
pairwise.wilcox.test(x = datos$Bio11, g = datos$Group, p.adjust.method = "bonf" )
pairwise.wilcox.test(x = datos$Bio14, g = datos$Group, p.adjust.method = "bonf" )
pairwise.wilcox.test(x = datos$Bio17, g = datos$Group, p.adjust.method = "bonf" )
pairwise.wilcox.test(x = datos$Bio18, g = datos$Group, p.adjust.method = "bonf" )



