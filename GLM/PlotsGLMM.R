library(multcomp) # Simultaneous Inference in General Parametric Models CRAN v1.4-25
library(lme4) # Linear Mixed-Effects Models using 'Eigen' and S4 CRAN v1.1-35.2
library(lmtest) # Testing Linear Regression Models CRAN v0.9-40
library(lattice) # Trellis Graphics for R CRAN v0.21-9
library(tidyverse) # Easily Install and Load the 'Tidyverse' CRAN v2.0.0 
library(emmeans) # Estimated Marginal Means, aka Least-Squares Means CRAN v1.10.1
library(ggeffects) # Create Tidy Data Frames of Marginal Effects for 'ggplot' from Model Outputs CRAN v1.6.0
library(psych) # Procedures for Psychological, Psychometric, and Personality Research CRAN v2.4.3
library(ggridges) # Ridgeline Plots in 'ggplot2' CRAN v0.5.6
library(lazyWeave) # LaTeX Wrappers for R Users CRAN v3.0.2 
library(cowplot) # Streamlined Plot Theme and Plot Annotations for 'ggplot2' CRAN v1.1.3
library(data.table) # Extension of `data.frame` CRAN v1.15.4
library(ggpubr) # 'ggplot2' Based Publication Ready Plots CRAN v0.6.0
library(MASS) # Support Functions and Datasets for Venables and Ripley's MASS CRAN v7.3-60
library(car) # Companion to Applied Regression CRAN v3.1-2
library(stats) 
library(nlme) # Linear and Nonlinear Mixed Effects Models CRAN v3.1-166
library(vcd) # Visualizing Categorical Data CRAN v1.4-13
library(tidyr) # Tidy Messy Data CRAN v1.3.1
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics CRAN v3.5.1
library(dplyr) # A Grammar of Data Manipulation CRAN v1.1.4
library(gg.layers) # ggplot layers [github::rpkgs/gg.layers] v0.1.3

# Read data
GLMMDB <- read.csv("GLM/GLM_DungBeetles.csv", header = T, sep = ";")

# Convert factor -----
GLMMDB$VegCover <- factor(GLMMDB$VegCover, levels = c("Secondary forest", 
                                                      "Cacao agroforestry systems",
                                                      "Pasture"))

#Plot1 -----

## Abundance ----
plotDB <- ggplot(GLMMDB, aes(x = VegCover, y = Abundance, 
                            fill = VegCover))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_colour_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  scale_fill_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  theme(axis.text.y = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.y = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  scale_y_continuous("Number of individuals", breaks = seq(0, 50, by = 10))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "none")+
  theme(legend.title = element_blank(), legend.key = element_blank())+
  annotate("text", label = "A", size = 4, x = 1, y = 30, family = "serif")+
  annotate("text", label = "B", size = 4, x = 2, y = 45, family = "serif")+
  annotate("text", label = "C", size = 4, x = 3, y = 10, family = "serif")+
  annotate("text", label = "A", size = 6, x = 3.5, y = 45, family = "serif")+
  xlab("Vegetal cover")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Species ----
plotSpe <- ggplot(GLMMDB, aes(x = VegCover, y = Species, 
                             fill = VegCover))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_colour_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  scale_fill_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  theme(axis.text.y = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.y = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  scale_y_continuous("Number of species", breaks = seq(0, 4, by = 1))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "none")+
  theme(legend.title = element_blank(), legend.key = element_blank())+
  annotate("text", label = "A", size = 4, x = 1, y = 4, family = "serif")+
  annotate("text", label = "A", size = 4, x = 2, y = 4, family = "serif")+
  annotate("text", label = "B", size = 4, x = 3, y = 2, family = "serif")+
  annotate("text", label = "B", size = 6, x = 3.5, y = 4, family = "serif")+
  xlab("Vegetal cover")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Longitude ----
plotLon <- ggplot(GLMMDB, aes(x = VegCover, y = Longitude, 
                              fill = VegCover))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_colour_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  scale_fill_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  theme(axis.text.y = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.y = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.x = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.x = element_text(size = 16, family = "serif", face="bold"))+
  scale_y_continuous("Longitude (mm)", breaks = seq(0, 20, by = 4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "none")+
  theme(legend.title = element_blank(), legend.key = element_blank())+
  annotate("text", label = "A", size = 4, x = 1, y = 19, family = "serif")+
  annotate("text", label = "B", size = 4, x = 2, y = 12, family = "serif")+
  annotate("text", label = "A", size = 4, x = 3, y = 20, family = "serif")+
  annotate("text", label = "C", size = 6, x = 3.5, y = 20, family = "serif")+
  xlab("Vegetation cover")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Biomass ----
plotBiom <- ggplot(GLMMDB, aes(x = VegCover, y = Biomass, 
                              fill = VegCover))+
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_colour_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  scale_fill_manual(values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  theme(axis.text.y = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.y = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.x = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.x = element_text(size = 16, family = "serif", face="bold"))+
  scale_y_continuous("Biomass (g)", breaks = seq(0, 20, by = 4))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "none")+
  theme(legend.title = element_blank(), legend.key = element_blank())+
  annotate("text", label = "A", size = 4, x = 1, y = 15.5, family = "serif")+
  annotate("text", label = "B", size = 4, x = 2, y = 5, family = "serif")+
  annotate("text", label = "B", size = 4, x = 3, y = 3.5, family = "serif")+
  annotate("text", label = "D", size = 6, x = 3.5, y = 16, family = "serif")+
  xlab("Vegetation cover")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Join plots ----
ggarrange(plotDB, plotSpe, plotLon, plotBiom, 
          ncol = 2, nrow = 2, align = "v")

## Save Plot ----
ggsave("plot1.png", width = 14, height = 8)

# Boxplots ----

## Size ----
Size <- GLMMDB |> 
  pivot_longer(cols = c(Large, Medium, Small), 
               names_to = "Types", 
               values_to = "Total") |>
  mutate(Types = factor(Types, 
                        levels = c("Large", "Medium", "Small"))) |> 
  ggplot(aes(x = VegCover, y = Total, fill = Types)) + 
  geom_boxplot2(width = 0.6, width.errorbar = 0.1, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_fill_manual(labels = c("Large", "Medium", "Small"),
                    values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.y = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.y = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.x = element_blank())+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "top", 
        legend.text = element_text(family = "sans", size = 12), 
        legend.title = element_text(family = "sans", size = 12))+
  annotate("text", label = "A", size = 4, x = 0.8, y = 12, family = "serif")+
  annotate("text", label = "a", size = 4, x = 1, y = 10, family = "serif")+
  annotate("text", label = "A", size = 4, x = 1.2, y = 9, family = "serif", 
           fontface = "italic")+
  annotate("text", label = "B", size = 4, x = 1.8, y = 4, family = "serif")+
  annotate("text", label = "a", size = 4, x = 2, y = 8, family = "serif")+
  annotate("text", label = "B", size = 4, x = 2.2, y = 38, family = "serif",
           fontface = "italic")+
  annotate("text", label = "B", size = 4, x = 2.8, y = 3, family = "serif")+
  annotate("text", label = "b", size = 4, x = 3, y = 6, family = "serif")+
  annotate("text", label = "A", size = 4, x = 3.2, y = 3, family = "serif", 
           fontface = "italic")+
  annotate("text", label = "A", size = 6, x = 3.5, y = 40, family = "serif")+
  ylab("Number of individuals")+
  labs(fill = "Size")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Nesting ----
NestDB <- GLMMDB |> 
  pivot_longer(cols = c(Dweller, Roller, Tunneler), 
               names_to = "Types", 
               values_to = "Total") |>
  mutate(Types = factor(Types, 
                        levels = c("Dweller", "Roller", "Tunneler"))) |> 
  ggplot(aes(x = VegCover, y = Total, fill = Types)) + 
  geom_boxplot2(width = 0.6, width.errorbar = 0.1, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_fill_manual(labels = c("Dweller", "Roller", "Tunneler"),
                    values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.y = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "top", 
        legend.text = element_text(family = "sans", size = 12), 
        legend.title = element_text(family = "sans", size = 12))+
  annotate("text", label = "A", size = 4, x = 0.8, y = 12, family = "serif")+
  annotate("text", label = "a", size = 4, x = 1, y = 2, family = "serif")+
  annotate("text", label = "A", size = 4, x = 1.2, y = 27, family = "serif", 
           fontface = "italic")+
  annotate("text", label = "A", size = 4, x = 1.8, y = 7, family = "serif")+
  annotate("text", label = "a", size = 4, x = 2, y = 2, family = "serif")+
  annotate("text", label = "B", size = 4, x = 2.2, y = 40, family = "serif",
           fontface = "italic")+
  annotate("text", label = "B", size = 4, x = 2.8, y = 3, family = "serif")+
  annotate("text", label = "a", size = 4, x = 3, y = 2, family = "serif")+
  annotate("text", label = "C", size = 4, x = 3.2, y = 8, family = "serif", 
           fontface = "italic")+
  annotate("text", label = "B", size = 6, x = 3.5, y = 40, family = "serif")+
  ylab("Number of individuals")+
  labs(fill = "Food relocation")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Dwellers ----
Dwell <- GLMMDB |> 
  pivot_longer(cols = c(MediumDweller, SmallDweller), 
               names_to = "Types", 
               values_to = "Total") |>
  mutate(Types = factor(Types, 
                        levels = c("MediumDweller", "SmallDweller"))) |> 
  ggplot(aes(x = VegCover, y = Total, fill = Types)) + 
  geom_boxplot2(width = 0.6, width.errorbar = 0.1, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_fill_manual(labels = c("Medium", "Small"),
                    values = c("darkgoldenrod3","darkolivegreen2"))+
  theme(axis.title.x = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.y = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.y = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.x = element_text(size = 12, family = "serif", face="bold"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "top", 
        legend.text = element_text(family = "sans", size = 12), 
        legend.title = element_text(family = "sans", size = 12))+
  scale_y_continuous(breaks = seq(0, 8, by = 2))+
  annotate("text", label = "A", size = 4, x = 0.85, y = 9, family = "serif")+
  annotate("text", label = "a", size = 4, x = 1.15, y = 1.2, family = "serif")+
  annotate("text", label = "A", size = 4, x = 1.85, y = 5.3, family = "serif")+
  annotate("text", label = "a", size = 4, x = 2.15, y = 1, family = "serif")+
  annotate("text", label = "B", size = 4, x = 2.85, y = 2.3, family = "serif")+
  annotate("text", label = "a", size = 4, x = 3.15, y = .5, family = "serif")+
  annotate("text", label = "C", size = 6, x = 3.5, y = 9, family = "serif")+
  ylab("Number of individuals")+
  labs(fill = "Dwellers size")+
  xlab("Vegetation cover")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Tunnellers ----
Tun <- GLMMDB |> 
  pivot_longer(cols = c(LargeTunn, MediumTunn, SmallTunn), 
               names_to = "Types", 
               values_to = "Total") |>
  mutate(Types = factor(Types, 
                        levels = c("LargeTunn", "MediumTunn", "SmallTunn"))) |> 
  ggplot(aes(x = VegCover, y = Total, fill = Types)) + 
  geom_boxplot2(width = 0.6, width.errorbar = 0.1, alpha = 0.6)+
  geom_point(alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.6), color = "cyan4")+
  scale_fill_manual(labels = c("Large", "Medium", "Small"),
                    values = c("darkolivegreen", "darkgoldenrod3","darkolivegreen2"))+
  theme(axis.title.x = element_text(size = 16, family = "serif", face="bold"))+
  theme(axis.text.y = element_text(size = 12, family = "serif", face="bold"))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.x = element_text(size = 12, family = "serif", face="bold"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(panel.background = element_blank(), legend.position = "top", 
        legend.text = element_text(family = "sans", size = 12), 
        legend.title = element_text(family = "sans", size = 12))+
  scale_y_continuous(breaks = seq(0, 40, by = 10))+
  annotate("text", label = "A", size = 4, x = 0.8, y = 12, family = "serif")+
  annotate("text", label = "a", size = 4, x = 1, y = 4, family = "serif")+
  annotate("text", label = "A", size = 4, x = 1.2, y = 9, family = "serif", 
           fontface = "italic")+
  annotate("text", label = "A", size = 4, x = 1.8, y = 5, family = "serif")+
  annotate("text", label = "a", size = 4, x = 2, y = 4, family = "serif")+
  annotate("text", label = "B", size = 4, x = 2.2, y = 38, family = "serif",
           fontface = "italic")+
  annotate("text", label = "B", size = 4, x = 2.8, y = 4, family = "serif")+
  annotate("text", label = "a", size = 4, x = 3, y = 4, family = "serif")+
  annotate("text", label = "A", size = 4, x = 3.2, y = 4, family = "serif", 
           fontface = "italic")+
  annotate("text", label = "D", size = 6, x = 3.5, y = 38, family = "serif")+
  ylab("Number of individuals")+
  labs(fill = "Tunnelers size")+
  xlab("Vegetation cover")+
  theme(panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))


## Plot Total ----
ggarrange(Size, NestDB, Dwell, Tun,
          ncol = 2, nrow = 2, align = "v")

## Save plot ----
ggsave("plot2.png", width = 12, height = 8)

