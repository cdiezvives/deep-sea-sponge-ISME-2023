#' Data analysis scripts for manuscript XXX
#' 
#' 16S comparison among deep-sea sponge species and locations
#' 

library(vegan)
library(pheatmap)
library(RColorBrewer)
library (tidyverse)
library(MASS)
library(ggrepel)
library(ggplot2)
library(effects)
library(scales)

options(digits = 2)

# Import data from mothur ----
counts_115 <-
  readr::read_delim(
    'final_sponge_asv.shared',
    delim = '\t',
    col_types = readr::cols('Group' = 'c', .default = 'i')
  )

taxa_115 = read.table(
  'sponge.trim.contigs.pcr.good.unique.good.filter.unique.precluster.abund.pick.pick.asv.ASV.cons.taxonomy',
  header = T,
  row.names = 1,
  stringsAsFactors = F
)

Infotable_115 <-
  readr::read_csv('data-analysis-ready/16s-deep-shallow-species/Infotable115.csv')

# Relative abundance and presence absence transformations
counts_RA_115 <-
  t(t(counts_115) / colSums(counts_115) * 100)

counts_PA_115 <- counts_115
counts_PA_115[counts_PA_115 > 0 ] <- 1

# Rarefaction ----
set.seed(1010)

counts_rare_115 <- rrarefy(t(counts_115), 20000)
counts_rare_115 <- as.data.frame(t(counts_rare_115))

# Remove empty ASV after rarefaction
counts_rare_115 <- counts_rare_115[ rowSums(counts_rare_115) > 0, ]
counts_rare_115_RA <- as.data.frame(t(t (counts_rare_115) / colSums(counts_rare_115) *100))

taxa_rare_115 <- taxa_115[ rownames(counts_rare_115), ] 

# Dendrogram & Heatmap ----
# Dendrogram at Genus level
levelall = apply(counts_rare_115_RA, 2, function(x) tapply(x, taxa_rare_115$alllevelstaxa, sum))
colSums(levelall)

levelall.imp  <- levelall
levelall.imp <- levelall.imp[rowSums(levelall.imp) > 1, ]

otherall <- levelall[rowSums(levelall) <= 1, , drop=F]
colSums(otherall)
levelall.imp.all <- rbind(levelall.imp, "Bacteria;Xothers" = colSums(otherall))

counts_bc <- vegdist(log2(t(levelall)+1), method="bray")  
hc <- hclust(counts_bc, method = "complete")
plot(rev(hc), hang = -1)

# save the order for the barplot
# hc$labels[hc$order] -> manually rearranged nodes
orden115 <- hc$labels[hc$order]

# barplot at level 2 to simplify colors
level3 = apply(counts_rare_115_RA, 2, function(x) tapply(x, taxa_rare_115$level3taxa, sum)) 
level3.imp  <- level3 
level3.imp <- level3.imp[rowSums(level3) > 5, ]
other3 <- level3[rowSums(level3) <= 5, , drop=F]
colSums(other3)
level3.imp.all <- rbind(level3.imp, "Bacteria; Others < 5%" = colSums(other3)  )
colSums(level3)

# order by dendrogram
Infotable_plot <- as.data.frame(Infotable_115)
rownames(Infotable_plot) <- Infotable_115$Species_ID
level.plot <- level3.imp[,orden115]
Infotable_plot <- Infotable_plot[ orden115,]
identical (colnames(level.plot ), rownames(Infotable_plot))

RdBu_ramp <-
  colorRampPalette(colors = RColorBrewer::brewer.pal(11, 'RdBu'))

pheatmap(log10(level.plot + 1),
         fontsize_row = 5, fontsize_col = 5,
         clustering_method = "average", 
         cluster_cols = F,
         cluster_rows = T,
         color = c(rev(RdBu_ramp(99))),
         annotation_col = Infotable_plot[ ,c("Depth_code_2lev", "Order", "Location_collected"), drop=F]) 

# Distance to centroid for order-depth groups by microbiome at genus level
Infotable_115 <- 
  Infotable_115 %>% 
  mutate(Order_depth= paste(Order, Depth_code_2lev, sep = "_"))

mod <- betadisper(counts_bc , Infotable_115$Order_depth, 
                  type = c("median"), bias.adjust = FALSE,
                  sqrt.dist = FALSE, add = FALSE)
boxplot(mod)
anova(mod)

# Distance to centroid for deep vs. shallow groups by microbiome at genus level
mod <- betadisper(counts_bc , Infotable_115$Depth_code_2lev, 
                  type = c("median"), bias.adjust = FALSE,
                  sqrt.dist = FALSE, add = FALSE)
boxplot(mod)
anova(mod)

# Alpha diversity deep vs. shallow
shannonH <- diversity(t(levelall), index = "shannon")

ggplot(Infotable_115, aes(x = Depth_code_2lev, y = shannonH)) + 
  geom_boxplot() +   
  theme_light() +
  theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

shannonH_model <- lm(shannonH ~ Depth_code_2lev , data = Infotable_115)
anova(shannonH_model)
plot(allEffects(shannonH_model))


# Ordination of 115 samples by Order and depth based on microbiome at genus level
mod <- betadisper(counts_bc , Infotable_115$Depth_code_2lev, 
                  type = c("median"), bias.adjust = FALSE,
                  sqrt.dist = FALSE, add = FALSE)
boxplot(mod)
anova(mod)

mod_pcoa_pcts <- mod$eig/sum(mod$eig) * 100
mod_pcoa_pcts <- round(mod_pcoa_pcts, 2)

site_scores <-
  scores(mod, choices = c(1, 2))$sites %>%
  as_tibble(rownames = 'site_id') %>% 
  left_join(
    Infotable_115 %>% 
      dplyr::select(Species_ID, Depth_code_2lev, Order, Location_collected),
    by = c(site_id = 'Species_ID')) %>% 
  mutate(Order = factor(Order),
         Depth_code_2lev = factor(Depth_code_2lev))

par(pty = 's')

plot(
  mod,
  hull = FALSE,
  ellipse = TRUE,
  label = FALSE,
  conf = 0.70,
  axes = c(1, 2),       ### select axes
  pch = NA , # not plot points
  xlab = sprintf('PCoA_1 [%s %%]', mod_pcoa_pcts[1]),
  ylab = sprintf('PCoA_2 [%s %%]', mod_pcoa_pcts[2]) ## select axes
)

points(site_scores$PCoA1, site_scores$PCoA2, ## select axes
       pch = c(16,17)[site_scores$Depth_code_2lev],
       col = site_scores$Order)

legend(x = 'topright',
       legend = levels(site_scores$Depth_code_2lev),
       pch = c(16,17))

legend(x = 'bottomright',
       legend = levels(site_scores$Order),
       pch = 16,
       col = palette()[1:4])

# Ordination of 115 samples by Species based on microbiome at genus level
mod <- betadisper(counts_bc , Infotable_115$Species, 
                  type = c("median"), bias.adjust = FALSE,
                  sqrt.dist = FALSE, add = FALSE)

mod_pcoa_pcts <- mod$eig/sum(mod$eig) * 100
mod_pcoa_pcts <- round(mod_pcoa_pcts, 2)

site_scores <-
  scores(mod, choices = c(1, 2))$sites %>%
  as_tibble(rownames = 'site_id') %>% 
  left_join(
    Infotable_115 %>% 
      dplyr::select(Species_ID, Depth_code_2lev, Order, Species),
    by = c(site_id = 'Species_ID')) %>% 
  mutate(Order = factor(Order),
         Depth_code_2lev = factor(Depth_code_2lev),
         Species = factor(Species))

par(pty = 's')

plot(
  mod,
  hull = FALSE,
  ellipse = TRUE,
  label = FALSE,
  conf = 0.70,
  axes = c(1, 2),       ### select axes
  pch = NA , # not plot points
  xlab = sprintf('PCoA_1 [%s %%]', mod_pcoa_pcts[1]),
  ylab = sprintf('PCoA_2 [%s %%]', mod_pcoa_pcts[2]) ## select axes
)

points(site_scores$PCoA1, site_scores$PCoA2, ## select axes
       pch = c(16,17)[site_scores$Depth_code_2lev],
       col = site_scores$Species)

legend(x = 'topright',
       legend = levels(site_scores$Depth_code_2lev),
       pch = c(16,17))

legend(x = 'bottomright',
       legend = levels(site_scores$Species),
       pch = 16,
       col = palette()[1:4])
