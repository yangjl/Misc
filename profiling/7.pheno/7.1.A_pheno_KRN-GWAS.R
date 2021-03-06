### Jinliang Yang
### calculate the mixed linear model and interaction terms

namlong <- read.csv("data/pheno_NAMRIL_KRN_051612.csv")
namlong$KRN <- as.numeric(namlong$KRN)
namlong <- namlong[!is.na(namlong$KRN),]

library("lme4")
m1 <- lmer(KRN ~ Genotype + (1 |block) + (1 |plot) + (1 |Year) + (1 |Location), data=namlong)
m2 <- lmer(KRN ~ Genotype + (1 |block) + (1 |plot) + (1 |Year) + (1 |Location) + (1 |Genotype:Year) + (1 |Genotype:Location), data=namlong)

m3 <- lmer(KRN ~ Genotype + (1 |Year) + (1 | Location), data=namlong)
m4 <- lmer(KRN ~ Genotype + (1 |Year) + (1 | Location) + (1 |Genotype:Year) + (1 |Genotype:Location), data=namlong)

save(file="cache/krn_lmer.RData", list=c("m1", "m2", "m3", "m4"))

ob <- load("cache/krn_lmer.RData")

res1 <- anova(m1, m2)

