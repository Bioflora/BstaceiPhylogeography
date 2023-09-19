#############################################
##### PCoA 3D Generation with SNPs data #####
######## Miguel Campos Cáceres 2023 ########
#############################################

library(dartR)
library(plotly)
library(stats)
library("SNPRelate")
library(ape)

#Load data
Stacei <- gl.read.vcf('Stacei_Final.vcf', verbose = 2)

## PCoA
pcoa <- gl.pcoa(Stacei)

# Export the values for the PCoA 3D you've to add a new column with colors manually
pcoa$scores
write.csv(pcoa$scores,"./PCoAscores.csv",row.names = F)

# To obtain the percentaje of variance of each axis
var_frac <- pcoa$eig/sum(pcoa$eig)
signif(sum(var_frac[1:3]) * 100, 4)
signif(sum(var_frac[1:2]) * 100, 4)
signif(sum(var_frac[1:1]) * 100, 4)

###  PCoA 3D
library(plotly)
PCoA <- read.csv(file="PCoAscores_mod.csv", header=TRUE, sep=";")
trace1 <- list(
  mode = "markers",
  name = "Brachypodium stacei",
  type = "scatter3d",
  x = PCoA$PC1,
  y = PCoA$PC2,
  z = PCoA$PC3,
  marker = list(
    size = 6,
    color = PCoA$Colour,
    opacity = 1
  ),
  inherit = FALSE,
  text = PCoA$Pop
)

data <- list(trace1)
layout <- list(
  scene = list(
    xaxis = list(
      title = "PA1 (10.04%)",
      mirror = TRUE,
      showline = TRUE,
      gridcolor = "rgba(211,211,211,1)",
      gridwidth = 2,
      linecolor = "rgba(169,169,169,1)",
      linewidth = 4,
      showticksuffix = "none"
    ),
    yaxis = list(
      title = "PA2 (9.45%)",
      mirror = TRUE,
      showline = TRUE,
      gridcolor = "rgba(211,211,211,1)",
      gridwidth = 2,
      linecolor = "rgba(169,169,169,1)",
      linewidth = 4
    ),
    zaxis = list(
      title = "PA3 (8.06%)",
      mirror = TRUE,
      showline = TRUE,
      gridcolor = "rgba(211,211,211,1)",
      gridwidth = 2,
      linecolor = "rgba(169,169,169,1)",
      linewidth = 4
    ),
    camera = list(eye = list(
      x = 0,
      y = 0,
      z = 2.5
    )),
    bgcolor = "#E1DCE9HH"
  ),
  title = "PCoA based on SNP markers"
)
p <- plot_ly()
p <- add_trace(p, mode=trace1$mode, name=trace1$name, type=trace1$type, x=trace1$x, y=trace1$y, z=trace1$z, marker=trace1$marker, inherit=trace1$inherit, text=trace1$text)
p <- layout(p, scene=layout$scene, title=layout$title)
p 
