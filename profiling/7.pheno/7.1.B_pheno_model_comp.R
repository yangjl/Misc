### Jinliang Yang
### calculate the mixed linear model and interaction terms


ob <- load("cache/krn_lmer.RData")

res1 <- anova(m1, m2)

save(file="cache/krn_models.RData", list="res1")

res2 <- anova(m3, m4)

save(file="cache/krn_models.RData", list=c("res1", "res2"))

     