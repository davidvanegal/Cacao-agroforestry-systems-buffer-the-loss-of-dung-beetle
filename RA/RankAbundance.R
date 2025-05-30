source("RA/DetAbu.R")
source("RA/UndAbu.R")
source("RA/SpecDist.R")
library(iNEXT)
library(Rcpp)
library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(readxl)

## Read data ----
ra <- read.csv('RA/RangeAbundance.csv', sep=';')
View(ra)
ra2 <- data.frame(ra)
names(ra2)

# Group data per sites
ra2 <- ra2%>%
  group_by(Species)%>%
  summarise(Secondary.forest = sum(Secondary.forest), 
            Cacao.agroforestry.systems = sum(Cacao.agroforestry.systems), 
            Pasture = sum(Pasture))

# Convert to data frame
ra2 <- data.frame(ra2)

## PA ----
out1 <- SpecDist(ra2$Secondary.forest, "abundance")
out1 <- subset(out1, probability>0)
totalSecondary.forest <- dim(out1)[1]
i=1
for (i in 1:totalSecondary.forest) {
  out1$number[i] <- i
  out1$sp[i] <- as.numeric(rownames(out1))[i]
}
out1$VegetationCover  <- "Secondary forest"

## SPP ----
out2 <- SpecDist(ra2$Cacao.agroforestry.systems, "abundance")
out2 <- subset(out2, probability>0)
totalCacao.agroforestry.systems <- dim(out2)[1]
i=1
for (i in 1:totalCacao.agroforestry.systems) {
  out2$number[i] <- i
  out2$sp[i] <- as.numeric(rownames(out2))[i]
}
out2$VegetationCover  <- "Cacao agroforestry systems"

## VS ----
out3 <- SpecDist(ra2$Pasture, "abundance")
out3 <- subset(out3, probability>0)
totalPasture <- dim(out3)[1]
i=1
for (i in 1:totalPasture) {
  out3$number[i] <- i
  out3$sp[i] <- as.numeric(rownames(out3))[i]
}
out3$VegetationCover  <- "Pasture"

## Create a single base
outtotal <- rbind(out1, out2, out3)
outtotal <-subset(outtotal, method %in% c("detected"))

outtotal$VegetationCover <- factor(outtotal$VegetationCover, 
                               levels = c("Secondary forest", 
                                          "Cacao agroforestry systems", 
                                          "Pasture"))

## Plot ----
ggplot(outtotal, aes(x = number, y = probability, 
                     color = VegetationCover, shape = VegetationCover))+ 
  geom_point(size = 4)+ 
  geom_line(lwd = 1, linetype = 3)+
  scale_color_manual(values = c("lightgoldenrod2", "orange","orangered"))+
  scale_x_continuous("Species rank")+
  scale_y_continuous("Proportional abundance")+
  theme(axis.title.x = element_text(family = "serif", size = 16, face="bold"))+
  theme(axis.text.x = element_text(family = "serif", size = 12, face="bold"))+
  theme(axis.title.y = element_text(family = "serif", size = 16, face="bold"))+
  theme(axis.text.y = element_text(family = "serif", size = 12, face="bold"))+
  facet_grid(cols = vars(VegetationCover))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), text = element_text(size = 14))+
  theme(strip.background = element_blank(), strip.text.x = element_blank())+
  theme(legend.position="top")+
  labs(col = "Vegetation cover", shape = "Vegetation cover")+
  theme(panel.grid.major.y = element_line(color = "gray95", linetype = "dashed"))

ggsave("RA/RangeAbundance.png", width = 12, height = 7)
