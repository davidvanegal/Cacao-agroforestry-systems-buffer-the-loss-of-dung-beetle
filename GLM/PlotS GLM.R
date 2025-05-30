library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics CRAN v3.5.0
library(cowplot) # Streamlined Plot Theme and Plot Annotations for 'ggplot2' CRAN v1.1.3
library(ggrepel) # Automatically Position Non-Overlapping Text Labels with 'ggplot2' CRAN v0.9.5
library(ggimage) # Use Image in 'ggplot2' CRAN v0.3.3
library(png) # Read and write PNG images CRAN v0.1-8 
library(patchwork) # The Composer of Plots CRAN v1.2.0 

DB <- read.csv("GLM_DungBeetles.csv", header = T, sep = ";")
BE <- read.csv("GLM_Bees.csv", header = T, sep = ";")
BU <- read.csv("GLM_Butterflies.csv", header = T, sep = ";")
To <- read.csv("GLM_Total.csv", header = T, sep = ";")

imgDB <- readPNG("DungBeetle.png", native = TRUE) 
imgBE <- readPNG("Bee.png", native = TRUE) 
imgBU <- readPNG("Butterflies.png", native = TRUE) 

# Dung beetles -----
DB$VegCover <- ordered(DB$VegCover, levels = c("Secondary forest", 
                                               "Cacao agroforestry", 
                                               "Pastures"))

## Plot (A) Number of individuals----
IndDB <- ggplot(DB, aes(x = VegCover, y = Abundance))+
        geom_boxplot(width = 0.6, outlier.shape = NA)+
        geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
        theme(axis.title.x = element_blank())+
        theme(axis.text.x = element_blank())+
        scale_y_continuous("Number of individuals")+
        xlab("Vegetation cover")+
        theme(axis.text.x = element_text(size = 14))+
        theme(axis.title.x = element_text(size = 12))+
        theme(axis.text.y = element_text(size = 12))+
        theme(axis.title.y = element_text(size = 14))+
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"))+
        theme(panel.background = element_blank())+
        theme(legend.position = "none")+
        annotate("text", label = "A", size = 6, x = 1, y = 30)+
        annotate("text", label = "B", size = 6, x = 2, y = 45)+
        annotate("text", label = "C", size = 6, x = 3, y = 15)+
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank())+
  inset_element(p = imgDB, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.95)
IndDB

## Plot (B) Species -----------------------------------------------------------
SpDB <- ggplot(DB, aes(x = VegCover, y = Species))+
  geom_boxplot(width = 0.6, outlier.shape = NA)+
  geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous("Number of species")+
  xlab("Vegetation cover")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  theme(legend.position = "none")+
        annotate("text", label = "A", size = 6, x = 1, y = 8)+
        annotate("text", label = "A", size = 6, x = 2, y = 8)+
        annotate("text", label = "B", size = 6, x = 3, y = 3.5)+
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank())+
  inset_element(p = imgDB, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.95)
SpDB

## Plot (C) Longitude --------------------------------------------------------
LonDB <- ggplot(DB, aes(x = VegCover, y = Longitude))+
  geom_boxplot(width = 0.6, outlier.shape = NA)+
  geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous("Longitude (mm)")+
  xlab("Vegetation cover")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  theme(legend.position = "none")+
  annotate("text", label = "A", size = 6, x = 1, y = 20)+
  annotate("text", label = "B", size = 6, x = 2, y = 12)+
  annotate("text", label = "A", size = 6, x = 3, y = 20)+
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank())+
  inset_element(p = imgDB, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.95)
LonDB

## Plot (D) Biomass -------------------------------------------------------
BioDB <- ggplot(DB, aes(x = VegCover, y = Biomass))+
  geom_boxplot(width = 0.6, outlier.shape = NA)+
  geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous("Biomass (g)")+
  xlab("Vegetation cover")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  theme(legend.position = "none")+
  annotate("text", label = "A", size = 6, x = 1, y = 16)+
  annotate("text", label = "B", size = 6, x = 2, y = 5)+
  annotate("text", label = "B", size = 6, x = 3, y = 5)+
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank())+
  inset_element(p = imgDB, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.95)
BioDB

# Bees -----
BE$VegCover <- ordered(BE$VegCover, levels = c("Secondary forest", 
                                               "Cacao agroforestry", 
                                               "Pasture"))

## Plot (A) Number of individuals ----
IndBE <- ggplot(BE, aes(x = VegCover, y = Abundance))+
  geom_boxplot(width = 0.6, outlier.shape = NA)+
  geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous("Number of individuals")+
  xlab("Vegetation cover")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  theme(legend.position = "none")+
  annotate("text", label = "A", size = 6, x = 1, y = 3.2)+
  annotate("text", label = "A", size = 6, x = 2, y = 4)+
  annotate("text", label = "A", size = 6, x = 3, y = 4.5)+
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank())+
  inset_element(p = imgBE, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.95)
IndBE

## Plot (B) Species -----------------------------------------------------------
SpBE <- ggplot(BE, aes(x = VegCover, y = Species))+
  geom_boxplot(width = 0.6, outlier.shape = NA)+
  geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous("Number of species")+
  xlab("Vegetation cover")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  theme(legend.position = "none")+
  annotate("text", label = "B", size = 6, x = 1, y = 1.7)+
  annotate("text", label = "A", size = 6, x = 2, y = 2.5)+
  annotate("text", label = "B", size = 6, x = 3, y = 1.9)+
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank())+
  inset_element(p = imgBE, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.95)
SpBE

# Butterflies -----
BU$VegCover <- ordered(BU$VegCover, levels = c("Secondary forest", 
                                               "Cacao agroforestry", 
                                               "Pasture"))

## Plot (A) Number of individuals ----
IndBU <- ggplot(BU, aes(x = VegCover, y = Abundance))+
  geom_boxplot(width = 0.6, outlier.shape = NA)+
  geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous("Number of individuals")+
  xlab("Vegetation cover")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  theme(legend.position = "none")+
  annotate("text", label = "A", size = 6, x = 1, y = 3.2)+
  annotate("text", label = "A", size = 6, x = 2, y = 3)+
  annotate("text", label = "A", size = 6, x = 3, y = 4.5)+
  inset_element(p = imgBU, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.85)
IndBU

## Plot (B) Species -----------------------------------------------------------
SpBU <- ggplot(BU, aes(x = VegCover, y = Species))+
  geom_boxplot(width = 0.6, outlier.shape = NA)+
  geom_point(alpha = 0.3, 
             position=position_jitter(width = .2), color = "cyan4")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  scale_y_continuous("Number of species")+
  xlab("Vegetation cover")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank())+
  theme(legend.position = "none")+
  annotate("text", label = "A", size = 6, x = 1, y = 2.5)+
  annotate("text", label = "A", size = 6, x = 2, y = 3.2)+
  annotate("text", label = "A", size = 6, x = 3, y = 2.5)+
  inset_element(p = imgBU, left = 1, 
                bottom = 0.75, right = 0.9, top = 0.85)
SpBU

## Plots -------------------------------------------------------------------

##generate multiple plot
final_plot2 <- plot_grid(IndDB, SpDB, LonDB, BioDB, 
                         IndBE, SpBE,
                         IndBU, SpBU,
                        labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'),
                        align = 'vh', ncol = 2, nrow = 4)
png("Plot.png", width = 1800, height = 1626, res = 130)
print(final_plot2)
dev.off()





