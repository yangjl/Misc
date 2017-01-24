### Jinliang Yang
### 01-22-2017
### study NAM and 282 phenotypes


pheno <- read.csv("largedata/pheno_long.csv")

tb <- data.frame(table(pheno$trait))

library(g3tools)
library(lme4)

trait <- subset(pheno, trait %in% "TasselLength")
loc <- data.frame(table(trait$location))
loc <- subset(loc, Freq > 0)
fit <- get_BLUP(data = trait, model = value ~ (1 | genotype) + (1 | location) + (1| genotype:location), 
         which.factor = "Line", outfile = "largedata/BLUPs/TasselLength.csv")

