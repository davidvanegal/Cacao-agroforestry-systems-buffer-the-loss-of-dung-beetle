library(multcomp) # Simultaneous Inference in General Parametric Models CRAN v1.4-25
library(lme4) # Linear Mixed-Effects Models using 'Eigen' and S4 CRAN v1.1-35.2
library(lmtest) # Testing Linear Regression Models CRAN v0.9-40
library(car)

# Read data
GLMDB <- read.csv("GLM_DungBeetles.csv", header = T, sep = ";")
GLMBE <- read.csv("GLM_Bees.csv", header = T, sep = ";")
GLMBU <- read.csv("GLM_Butterflies.csv", header = T, sep = ";")

# Dung beetles -----
GLMDB$VegCover <- factor(GLMDB$VegCover)

## Analysis## Analfactorize()ysis
par(mfrow=c(2,4))
print(summary(GLMDB))
plot(GLMDB)

## Abundance 
IndDB <- glm(Abundance ~ VegCover, family = gaussian(), data = GLMDB)
print(summary(IndDB))
par(mfcol = c(2, 2))
plot(IndDB)

shapiro.test(IndDB$residuals)
leveneTest(log(Abundance+1) ~ VegCover, data = GLMDB)

anvIndDB <- anova(IndDB, test = "Chisq")
print(anvIndDB)
summary(anvIndDB)

tukIndDB <- glht(IndDB, linfct = mcp(VegCover = "Tukey"))
summary(tukIndDB)

## Species 
EspDB <- glm(log(Species+1) ~ VegCover, family = gaussian(), data = GLMDB)
print(summary(EspDB))
par(mfcol = c(2, 2))
plot(EspDB)

shapiro.test(EspDB$residuals)
leveneTest(log(Species+1) ~ VegCover, data = GLMDB)

anvEspDB <- anova(EspDB, test = "Chisq")
print(anvEspDB)
summary(anvEspDB)

tukEspDB <- glht(EspDB, linfct = mcp(VegCover = "Tukey"))
summary(tukEspDB)

# Longitude 
LonDB <- glm(log(Longitude+1) ~ VegCover, family = gaussian(), data = GLMDB)
print(summary(LonDB))
par(mfcol = c(2, 2))
plot(LonDB)

shapiro.test(LonDB$residuals)
leveneTest(log(Longitude+1) ~ VegCover, data = GLMDB)

anvLonDB <- anova(LonDB, test = "Chisq")
print(anvLonDB)
summary(anvLonDB)

tukLonDB <- glht(LonDB, linfct = mcp(VegCover = "Tukey"))
summary(tukLonDB)

# Biomass 
BiomDB <- glm(log(Biomass+1) ~ VegCover, 
              family = gaussian(), data = GLMDB)
print(summary(BiomDB))
par(mfcol = c(2, 2))
plot(BiomDB)

shapiro.test(BiomDB$residuals)
leveneTest(log(Biomass+1) ~ VegCover, data = GLMDB)

anvBiomDB <- anova(BiomDB, test = "Chisq")
print(anvBiomDB)
summary(anvBiomDB)

tukBiomDB <- glht(BiomDB, linfct = mcp(VegCover = "Tukey"))
summary(tukBiomDB)

# Bees -----
GLMBE$VegCover <- ordered(GLMBE$VegCover, levels = c("Secondary forest", 
                                                     "Cacao agroforestry", 
                                                     "Pastures"))

## Analysis
par(mfrow=c(2,4))
print(summary(GLMBE))
plot(GLMBE)

## Abundance 
IndBE <- glm(Abundance ~ VegCover, family = gaussian(), data = GLMBE)
print(summary(IndBE))
par(mfcol = c(2, 2))
plot(IndBE)

anvIndBE <- anova(IndBE, test = "Chisq")
print(anvIndBE)
summary(anvIndBE)

tukIndBE <- glht(IndBE, linfct = mcp(VegCover = "Tukey"))
summary(tukIndBE)

## Species 
EspBE <- glm(Species ~ VegCover, family = gaussian(), data = GLMBE)
print(summary(EspBE))
par(mfcol = c(2, 2))
plot(EspBE)

anvEspBE <- anova(EspBE, test = "Chisq")
print(anvEspBE)
summary(anvEspBE)

tukEspBE <- glht(EspBE, linfct = mcp(VegCover = "Tukey"))
summary(tukEspBE)

# Butterflies -----
GLMBU$VegCover <- ordered(GLMBU$VegCover, levels = c("Secondary forest", 
                                                     "Cacao agroforestry", 
                                                     "Pastures"))

## Analysis
par(mfrow=c(2,4))
print(summary(GLMBU))
plot(GLMBU)

## Abundance 
IndBU <- glm(Abundance ~ VegCover, family = gaussian(), data = GLMBU)
print(summary(IndBU))
par(mfcol = c(2, 2))
plot(IndBU)

anvIndBU <- anova(RichDB, test = "Chisq")
print(anvIndBU)
summary(anvIndBU)

tukIndBU <- glht(IndBU, linfct = mcp(VegCover = "Tukey"))
summary(tukIndBU)

## Species 
EspBU <- glm(Species ~ VegCover, family = gaussian(), data = GLMBU)
print(summary(EspBU))
par(mfcol = c(2, 2))
plot(EspBU)

anvEspBU <- anova(EspBU, test = "Chisq")
print(anvEspBU)
summary(anvEspBU)

tukEspBU <- glht(EspBU, linfct = mcp(VegCover = "Tukey"))
summary(tukEspBU)
