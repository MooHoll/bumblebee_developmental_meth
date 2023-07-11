setwd("~/Dropbox/Leicester_postdoc/Projects/Ben_Developmental_BB/weighted_meth_new_annotation")

library(lme4) # needed for lmer
library(emmeans)
library(car) # needed for Anova


methdev <- read.delim("gene_meth_by_stage_with_replicates.txt")
methdev$dev_stage<-as.factor(methdev$dev_stage) # insuring its treated as a factor in models


mixed.lmer <- lmer(weightedMeth ~ dev_stage + (1|replicate) + (1|gene_id), data = methdev) # treating replicate and gene_id as random effects in a mixed effect model

summary(mixed.lmer)

Anova(mixed.lmer) # to get a p value for the effect of dev_stage

emm1 = emmeans(mixed.lmer, specs = pairwise ~ dev_stage, lmer.df = "asymptotic") #https://cran.r-project.org/web/packages/emmeans/vignettes/sophisticated.html
emm1$contrasts
output<-as.data.frame(emm1$contrasts)# If you want to make a nice table














