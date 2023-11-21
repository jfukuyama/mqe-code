library(phyloseq)
library(adaptiveGPCA)
library(viridis)
library(RColorBrewer)
library(magrittr)
library(tidyverse)
tps = readRDS("Sprockett_Tsimane_ps.rds")
sum_nonzero = colSums(otu_table(tps) > 0)
tps_filt = prune_taxa(sum_nonzero > 5, tps)
tps_filt = prune_samples(sample_data(tps_filt)$Cohort == "Tsimane", tps_filt)
tps_filt = prune_samples(sample_data(tps_filt)$Age_Class != "Child", tps_filt)
otu_table(tps_filt) = log(1 + otu_table(tps_filt))
pp = processPhyloseq(tps_filt)
out.agpca = adaptivegpca(pp$X, pp$Q)
out.agpca$r

t_proposed <- function(x) {
    .5 - .5 * cos(x * pi)
}
rvec = t_proposed(seq(0,1,length.out = 30))

Qeig = eigen(pp$Q)
rvec_smaller = t_proposed(seq(0,1,length.out = 10))

out.ff = gpcaFullFamily(pp$X, pp$Q, k = 2, rvec = rvec, returnLong = TRUE, sampledata = sample_data(tps_filt))

out.ff$locations = out.ff$locations %>%
    mutate(Age_Saturated = pmin(1.6, Age_Years),
           Age_And_Type_Sat = ifelse(Sample_Type == "Feces", Age_Saturated, -Age_Saturated),
           r_label = paste("r =", round(1 - r, digits = 2)),
           Type = ifelse(Sample_Type == "Feces", "Gut", "Oral"))



## Figures for the paper
ggplot(subset(out.ff$locations, r %in% rvec[c(1, 5, 20, 30)]),
       aes(x = Axis1, y = Axis2, color = Age_And_Type_Sat, shape = Type)) +
    geom_point(size = .6) +
    facet_grid(. ~ fct_rev(r_label)) +
    scale_x_continuous(breaks = c(0, .1)) +
    scale_color_distiller(type = "div", palette = 5, guide = "none")
