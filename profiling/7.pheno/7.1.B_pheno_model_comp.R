### Jinliang Yang
### calculate the mixed linear model and interaction terms


ob <- load("cache/krn_lmer.RData")

res1 <- anova(m1, m2)

save(file="cache/krn_models.RData", list="res1")

res2 <- anova(m3, m4)

save(file="cache/krn_models.RData", list=c("res1", "res2"))

ob2 <- load("cache/krn_models.RData")


Genotype:Location (Intercept) 12.643  
Genotype:Year     (Intercept)  9.503  
plot              (Intercept) 11.521  
block             (Intercept)  0.000  
Location          (Intercept)  2.939  
Year              (Intercept) 10.477  
Residual                      32.055  

Groups   Name        Std.Dev.
plot     (Intercept) 16.611  
block    (Intercept)  1.929  
Location (Intercept)  1.595  
Year     (Intercept) 10.147  
Residual             32.121

