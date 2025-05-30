library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics CRAN v3.5.0
library(GGally) # Extension to 'ggplot2' CRAN v2.2.1
library(CCA) # Canonical Correlation Analysis CRAN v1.2.2
library(MASS) # Support Functions and Datasets for Venables and Ripley's MASS CRAN v7.3-60
library(vegan) # Community Ecology Package CRAN v2.6-4
library(ggforce) # Accelerating 'ggplot2' CRAN v0.4.2
library(png) # Read and write PNG images CRAN v0.1-8 
library(patchwork) # The Composer of Plots CRAN v1.2.0 
library(cowplot) # Streamlined Plot Theme and Plot Annotations for 'ggplot2' CRAN v1.1.3
library(grid)
library(ggrepel)

## Análisis bajo parámetros de abundancia

# Read data
DB <- read.csv("NMDS_DungBeetles.csv", header = T, sep = ";")
BE <- read.csv("NMDS_Bees.csv", header = T, sep = ";")
BU <- read.csv("NMDS_Butterflies.csv", header = T, sep = ";")

imgDB <- readPNG("DungBeetle.png", native = TRUE) 
imgBE <- readPNG("Bee.png", native = TRUE) 
imgBU <- readPNG("Butterflies.png", native = TRUE) 

DB <- BU
# Dung Beetles ----
dimDB <- dim(DB)
dimDB
DB3 <- DB[,2:length(DB)]

# ANOSIM ----
# Fire regime
analisis <- anosim(DB3, DB$VegCover, distance = "bray", permutations = 9999)
analisis

# Manipulación de datos
DB2 <- DB[1:dimDB[1],]
yDB <- DB2[,2:dimDB[2]]




# NMDS
nmdsDB <- metaMDS(yDB, k = 2, distance = "bray")
nmdsDB

# Cálculo de la bondad de ajuste para NMDS
bdanDB <- goodness(nmdsDB)
print(bdanDB)
sum(bdanDB)
mean(bdanDB)  

nmdsDB$points
nmds.scoresDB <-scores(nmdsDB)
nmds.scoresDB

SitioDB <- nmdsDB$points
sppDB <- nmdsDB$species
sppDB

# Determinar sitios
nmdsdataDB <- data.frame(nmds.scoresDB$sites)
nmdsdataDB$sites <- c(1:12)
nmdsdataDB$sites[1:4] <- c("Secondary forest")
nmdsdataDB$sites[5:8] <- c("Pasture")
nmdsdataDB$sites[9:12] <- c("Cacao agroforestry systems")

nmdsdataDB$sites <- factor(nmdsdataDB$sites, 
                           levels = c("Secondary forest", 
                                      "Cacao agroforestry systems", 
                                      "Pasture"))

# Plot
ggplot(nmdsdataDB, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(shape = sites, color = sites), size = 2.7) +
  geom_mark_ellipse(expand = 0, aes(color = sites, fill = sites))+
  scale_color_manual(values= c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  scale_fill_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  annotate("text", x = 1.4, y = 0.8, label = "Strees: 0.0774", size = 4)+
  theme(axis.title.x = element_text(family = "serif", size = 16, face = "bold"))+
  theme(axis.text.x = element_text(family = "serif", size = 12, face = "bold"))+
  theme(axis.title.y = element_text(family = "serif", size = 16, face = "bold"))+
  theme(axis.text.y = element_text(family = "serif", size = 12, face = "bold"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  labs(x = "NMDS1", y = "NMDS2",
       family = "Arial", size = 14)+
  theme(panel.background = element_blank(), 
        legend.text = element_text(family = "sans", size = 12), 
        legend.title = element_text(family = "sans", size = 12), 
        legend.position = "top")+
  labs(color = "Vegetation cover", shape = "Vegetation cover", 
       fill = "Vegetation cover")

ggsave("plotNMDS.png", width = 12, height = 7)

stressplot(nmdsDB)

View(nmdsdataDB)

# Bees ----
dimBE <- dim(BE)

# Manipulación de datos
BE2 <- BE[1:dimBE[1],]
yBE <- BE2[,2:dimBE[2]]

# NMDS
nmdsBE <- metaMDS(yBE, k = 2, distance = "bray")
nmdsBE

# Cálculo de la bondad de ajuste para NMDS
bdanBE <- goodness(nmdsBE)
print(bdanBE)
sum(bdanBE)
mean(bdanBE)  

nmdsBE$points
nmds.scoresBE <-scores(nmdsBE)
nmds.scoresBE

SitioBE <- nmdsBE$points
sppBE <- nmdsBE$species
sppBE

# Determinar sitios
nmdsdataBE <- data.frame(nmds.scoresBE$sites)
nmdsdataBE$sites <- c(1:12)
nmdsdataBE$sites[1:4] <- c("Secondary forest")
nmdsdataBE$sites[5:8] <- c("Pasture")
nmdsdataBE$sites[9:12] <- c("Cacao agroforestry")


# Plot
PlotBE <-ggplot(nmdsdataBE, aes(x = NMDS1, y = NMDS2, color = sites))+
  geom_point(aes(shape = sites), size=4) +
  geom_mark_ellipse(expand = 0, aes(color = sites))+
  scale_color_manual(values=c("cyan3", "cyan4", "black"))+
  annotate("text", x = -0.9, y = 0.7, label = "Strees: 0.1895")+
  theme(axis.text.x = element_text(size = 14, colour = "black"))+
  theme(axis.text.y = element_text(size = 14, colour = "black"))+
  theme(axis.title.x = element_text(size = 14, colour = "black"))+
  theme(axis.title.y = element_text(size = 14, colour = "black"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "none")+
  labs(x = "NMDS1", y = "NMDS2",
       family = "Arial", size = 14)+
  theme(legend.title=element_blank(), axis.title.x=element_blank(),
        text=element_text(size=20), legend.key=element_blank())+
  inset_element(p = imgBE, left = 1, 
                bottom = 0.85, right = 0.9, top = 1.05)

stressplot(nmdsBE)

View(nmdsdataBE)

# Butterflies ----
dimBU <- dim(BU)

# Manipulación de datos
BU2 <- BU[1:dimBU[1],]
yBU <- BU2[,2:dimBU[2]]

# NMDS
nmdsBU <- metaMDS(yBU, k = 2, distance = "bray")
nmdsBU

# Cálculo de la bondad de ajuste para NMDS
bdanBU <- goodness(nmdsBU)
print(bdanBU)
sum(bdanBU)
mean(bdanBU)  

nmdsBU$points
nmds.scoresBU <-scores(nmdsBU)
nmds.scoresBU

SitioBU <- nmdsBU$points
sppBU <- nmdsBU$species
sppBU

# Determinar sitios
nmdsdataBU <- data.frame(nmds.scoresBU$sites)
nmdsdataBU$sites <- c(1:12)
nmdsdataBU$sites[1:4] <- c("Secondary forest")
nmdsdataBU$sites[5:8] <- c("Pasture")
nmdsdataBU$sites[9:12] <- c("Cacao agroforestry")


# Plot
PlotBU <- ggplot(nmdsdataBU, aes(x = NMDS1, y = NMDS2, color = sites))+
  geom_point(aes(shape = sites), size=4) +
  geom_mark_ellipse(expand = 0, aes(color = sites))+
  annotate("text", x = -1, y = 0.7, label = "Strees: 0.1394")+
  scale_color_manual(values=c("cyan3", "cyan4", "black"))+
  theme(axis.text.x = element_text(size = 14, colour = "black"))+
  
  theme(axis.text.y = element_text(size = 14, colour = "black"))+
  theme(axis.title.x = element_text(size = 14, colour = "black"))+
  theme(axis.title.y = element_text(size = 14, colour = "black"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "none")+
  labs(x = "NMDS1", y = "NMDS2",
       family = "Arial", size = 14)+
  theme(legend.position="bottom", legend.title=element_blank(),
        text=element_text(size=20), legend.key=element_blank())+
  inset_element(p = imgBU, left = 1, 
                bottom = 0.80, right = 0.9, top = 0.90)

stressplot(nmdsBU)

View(nmdsdataBU)

## Plots -------------------------------------------------------------------

##generate multiple plot
final_plotNMDS <- plot_grid(PlotDB, PlotBE, PlotBU,
                         labels = c('A', 'B', 'C'),
                         align = 'vh', ncol = 1, nrow = 3)
png("PlotNMDS.png", width = 1500, height = 1858, res = 130)
print(final_plotNMDS)
dev.off()
